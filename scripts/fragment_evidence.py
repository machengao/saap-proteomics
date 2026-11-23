"""Calculate site evidence and Prosit spectral angles on the same observed MS2."""
import csv
import hashlib
import json
import math
import re
from pathlib import Path
from collections import defaultdict

import numpy as np
from pyteomics import mass

PROTON = mass.nist_mass["H+"][0][0]
WATER = mass.calculate_mass(formula="H2O")
# UNIMOD identifiers and monoisotopic deltas used in the frozen project competition.
MODS = {1:42.010564684, 4:57.021463721, 5:43.005813656, 7:0.984015583, 21:79.96633052,
        27:-18.010564684,28:-17.026549101,34:14.015650064,35:15.99491462,36:28.031300128,
        37:42.046950192,122:27.99491462,425:31.989829239,737:229.162932}
TOKEN = re.compile(r"([A-Z])(?:\[UNIMOD:(\d+)\])?")
SITE_MODS = {7:"NQ", 35:"M", 21:"STY", 34:"KR", 36:"KR", 37:"K", 122:"K", 425:"M"}
MQ_MODS = {"Oxidation (M)":35, "Deamidation (NQ)":7, "Acetyl (Protein N-term)":1,
           "Phospho (STY)":21, "Carbamidomethyl (C)":4}


def read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def encode(seq, mods=None, chemistry=None):
    mods = dict(mods or {})
    cam = chemistry is None or chemistry.cam
    if chemistry and chemistry.tmt:
        mods.setdefault(0,737)
        for i,a in enumerate(seq,1):
            if a=="K":
                mods.setdefault(i,737)
    return ("[UNIMOD:%s]-" % mods[0] if 0 in mods else "") + "".join(
        a + (f"[UNIMOD:{mods[i]}]" if i in mods else "[UNIMOD:4]" if a == "C" and cam else "")
        for i,a in enumerate(seq,1))


def mq_form(seq, modified):
    """Recover MQ variable modification positions, including nested parenthesized names."""
    s = modified.strip("_")
    cursor, plain, mods = 0, "", {}
    while cursor < len(s):
        if s[cursor].isupper():
            plain += s[cursor]
            cursor += 1
        elif s[cursor] == "(":
            end, depth = cursor+1, 1
            while end < len(s) and depth:
                depth += (s[end] == "(") - (s[end] == ")")
                end += 1
            name = s[cursor+1:end-1]
            if depth or name not in MQ_MODS or len(plain) in mods:
                raise ValueError("unsupported MQ modification: " + name)
            mods[len(plain)] = MQ_MODS[name]
            cursor = end
        else:
            raise ValueError("invalid modified sequence")
    if plain != seq:
        raise ValueError("base PSM sequence mismatch")
    return mods


def decode(form):
    nt = 0.
    if form.startswith("["):
        prefix = re.match(r"^\[UNIMOD:(\d+)\]-", form)
        if not prefix:
            raise ValueError("Unsupported terminal modification")
        nt, form = MODS[int(prefix[1])], form[prefix.end():]
    parts = list(TOKEN.finditer(form))
    if not parts or "".join(m[0] for m in parts) != form:
        raise ValueError("Unsupported ProForma: " + form)
    vals = np.array([mass.std_aa_mass[m[1]] + (MODS[int(m[2])] if m[2] else 0) for m in parts])
    return nt, vals


def neutral(form):
    nt, vals = decode(form)
    return float(vals.sum() + nt + WATER)


def ions(form,z):
    nt, vals = decode(form)
    peaks, labels = [], []
    for n, (suffix, prefix) in enumerate(zip(np.cumsum(vals[::-1])[:-1]+WATER,
                                             np.cumsum(vals)[:-1]+nt),1):
        for series, value in (("y",suffix),("b",prefix)):
            for charge in range(1, min(z,3)+1):
                peaks.append(value/charge + PROTON)
                labels.append(f"{series}{n}+{charge}")
    return np.array(peaks), np.array(labels)


def match(pm, obs, ppm):
    left = np.searchsorted(obs, pm*(1-ppm*1e-6))
    right = np.searchsorted(obs, pm*(1+ppm*1e-6))
    edges = sorted((abs(obs[j]-x)/x,i,j) for i,x in enumerate(pm)
                   for j in range(left[i],right[i]))
    used, paired = set(), {}
    for _,i,j in edges:
        if i not in paired and j not in used:
            paired[i] = j
            used.add(j)
    return paired


def angle(prediction, spectrum, ppm):
    pm, pi = map(np.asarray, prediction)
    obs,it,lo,hi = spectrum
    mask = (pm >= lo) & (pm <= hi) & (pi >= 0)
    pm, pi = pm[mask], pi[mask]
    if not len(pi) or not np.any(pi > 0):
        return None
    paired = match(pm,obs,ppm)
    v = np.zeros(len(pm))
    for i,j in paired.items():
        v[i] = it[j]
    a,b = np.sqrt(pi),np.sqrt(v)
    den = np.linalg.norm(a)*np.linalg.norm(b)
    return float(1-2/math.pi*math.acos(np.clip(float(np.dot(a,b)/den) if den else 0.,0,1)))


def fragment_score(form, z, spectrum, ppm):
    obs,it,lo,hi = spectrum
    pm, labels = ions(form,z)
    mask = (pm >= lo) & (pm <= hi)
    pm, labels = pm[mask],labels[mask]
    paired = match(pm,obs,ppm)
    observed = {j:str(labels[i]) for i,j in paired.items()}
    return {"matched_ions":len(paired),
            "explained_fraction":float(sum(it[j] for j in observed)/it.sum()), "observed":observed}


def spectrum(data):
    obs, it = np.asarray(data["mz"],float), np.asarray(data["intensity"],float)
    lo,hi = float(data["lower"]),float(data["upper"])
    if (obs.ndim != 1 or obs.shape != it.shape or not len(obs) or not np.all(np.isfinite(obs))
            or not np.all(np.isfinite(it)) or not np.all(it>=0) or not it.sum() or not 0<lo<hi):
        raise ValueError("Invalid centroid spectrum or acquisition bounds")
    order = np.argsort(obs)
    return obs[order],it[order],lo,hi


def mzml_peaks(path, wanted):
    from pyteomics import mzml
    result = {}
    with mzml.MzML(str(path), use_index=False) as spectra:
        for entry in spectra:
            number = re.search(r"(?:^|\s)scan=(\d+)(?:\s|$)",entry["id"])
            if not number or number[1] not in wanted or entry.get("ms level") != 2:
                continue
            n = number[1]
            if n in result:
                raise ValueError("Duplicate native scan: " + n)
            if "profile spectrum" in entry:
                raise ValueError("Use centroided mzML for peak matching")
            scan = entry["scanList"]["scan"][0]
            win = scan["scanWindowList"]["scanWindow"][0]
            ion = entry["precursorList"]["precursor"][0]["selectedIonList"]["selectedIon"][0]
            result[n] = {"mz":entry["m/z array"],"intensity":entry["intensity array"],
                         "charge":ion.get("charge state"),"lower":win["scan window lower limit"],
                         "upper":win["scan window upper limit"]}
    return result


def native_cache(directory, wanted, model, energy):
    """Load only relevant batches and require the matching model version, sequence, charge and CE."""
    result = {}
    for path in sorted(directory.glob(model+"_*.rows.json")):
        mapping = json.loads(path.read_text(encoding="utf-8"))
        selected = [(i,(r["peptide_sequences"],int(r["precursor_charges"]))) for i,r in enumerate(mapping)
                    if float(r["effective_CE"])==energy and (r["peptide_sequences"],int(r["precursor_charges"])) in wanted]
        if not selected:
            continue
        response_file = path.with_name(path.name.replace(".rows.json",".response.json"))
        if not response_file.exists():
            continue
        body = json.loads(response_file.read_text(encoding="utf-8"))
        if body["model_name"]!=model or str(body["model_version"])!="1":
            raise ValueError(f"Unexpected model/version in {response_file}")
        out = {x["name"]:np.asarray(x["data"]).reshape(x["shape"]) for x in body["outputs"] if x["name"] in ("mz","intensities")}
        if out["mz"].shape != out["intensities"].shape or len(out["mz"]) != len(mapping):
            raise ValueError("Prediction batch size mismatch")
        for i,key in selected:
            pair = (out["mz"][i],out["intensities"][i])
            if not all(np.all(np.isfinite(arr)) for arr in pair):
                raise ValueError("Nonfinite model prediction")
            if key in result and not all(np.array_equal(a,b) for a,b in zip(result[key],pair)):
                raise ValueError("Conflicting cached prediction")
            result[key] = pair
    return result

