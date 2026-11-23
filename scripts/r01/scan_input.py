"""Read direct MaxQuant DP observations, retaining physical scan identities."""
import csv
import hashlib
import re
import xml.etree.ElementTree as ET
from collections import Counter
from pathlib import PureWindowsPath


FIELDS = {
    "Base peptide sequence": "DP base sequence", "Modification": "DP modification",
    "Cluster index": "DP cluster index", "Cluster mass": "DP cluster mass",
    "Cluster mass SD": "DP cluster mass SD", "Cluster size total": "DP cluster size total",
    "Cluster size forward": "DP cluster size forward", "Cluster size reverse": "DP cluster size reverse",
    "Localization probability": "DP probabilities", "Max Probability": "DP positional probability",
    "Posterior error probability": "DP PEP", "Score": "DP score", "Mass difference": "DP mass difference",
    "Peptide Length Difference": "DP peptide length difference", "Protein IDs": "DP proteins",
    "Gene names": "DP gene names", "Protein names": "DP protein names",
    "Modification scan numbers": "DP mod scan number", "Aminoacid localizations": "DP AA",
}


def search_profile(xml):
    root = ET.parse(xml).getroot()
    groups = root.findall("./parameterGroups/parameterGroup")
    if not (root.findtext("maxQuantVersion") or "").startswith("2.8."):
        raise ValueError("Discovery supports the MaxQuant 2.8.x direct DP schema.")
    if root.findtext("dependentPeptides", "").lower() != "true" or not groups:
        raise ValueError("A dependent peptide search and parameter groups are required.")
    for group in groups:
        mods = [e.text for e in group.findall("fixedModifications/string")]
        if (mods != ["Carbamidomethyl (C)"] or group.findtext("multiplicity") != "1"
                or group.findtext("lcmsRunType") != "Standard" or group.findall("isobaricLabels/*")
                or any((e.text or "").strip() for e in group.findall("labelMods/string"))):
            raise ValueError("This discovery profile requires unlabelled Standard DDA and fixed CAM(C) only.")
    names = [PureWindowsPath(e.text).stem for e in root.findall("filePaths/string")]
    if not names or len(names) != len(set(names)):
        raise ValueError("Missing or ambiguous RAW names in mqpar.xml.")
    return names


def prepare(scans, names, output):
    index = {name: i for i, name in enumerate(names)}
    seen = set()
    counts = Counter()
    with (scans.open(encoding="utf-8-sig", newline="") as src,
          (output / "dp_scan_input.tsv").open("w", encoding="utf-8-sig", newline="") as dst,
          (output / "dp_review.tsv").open("w", encoding="utf-8", newline="") as log):
        reader = csv.DictReader(src, delimiter="\t")
        required = set(FIELDS.values()) | {"Raw file", "Scan number", "DP base raw file", "DP decoy"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing msmsScans columns: {sorted(missing)}")
        fields = list(FIELDS)
        fields[8:8] = ["Base peptide raw file", "Matches in Raw files"]
        writer = csv.DictWriter(dst, fields + ["Found by"], delimiter="\t")
        writer.writeheader()
        review = csv.writer(log, delimiter="\t")
        review.writerow(["raw_file", "scan_number", "reason"])
        for r in reader:
            name, scan = r["Raw file"], int(r["Scan number"])
            key = (name, scan)
            if key in seen or name not in index:
                raise ValueError(f"Duplicate or unknown scan key: {key}")
            seen.add(key)
            base = r["DP base sequence"]
            if not base:
                continue
            counts["dp_scans"] += 1
            reasons = []
            if any(";" in r[k] for k in ("DP mass difference", "DP PEP", "DP positional probability",
                                         "DP mod scan number", "DP peptide length difference")):
                reasons.append("non_scalar_DP_observation")
            if r["DP decoy"] == "+":
                reasons.append("DP_decoy")
            loc = r["DP probabilities"].replace("_", "")
            if re.sub(r"\([0-9.eE+\-]+\)", "", loc) != base:
                reasons.append("localization_sequence_mismatch")
            if PureWindowsPath(r["DP base raw file"]).stem != name:
                reasons.append("unexpected_DP_base_RAW")
            if r["DP mod scan number"] != str(scan):
                reasons.append("DP_mod_scan_mismatch")
            if reasons:
                counts["excluded_dp_scans"] += 1
                review.writerow([name, scan, ";".join(reasons)])
                continue
            normalized = {dest: r[source] for dest, source in FIELDS.items()}
            normalized.update({"Base peptide raw file": index[name], "Matches in Raw files": index[name],
                               "Found by": "direct_msmsScans_DP"})
            writer.writerow(normalized)
    return {"input_ms2_rows": len(seen), "dp_scans": counts["dp_scans"],
            "excluded_dp_scans": counts["excluded_dp_scans"], "raw_order": names}


def label_strict(output, names):
    spectra, events = set(), set()
    with ((output / "saap_candidates_strict.tsv").open(encoding="utf-8", newline="") as src,
          (output / "candidates.tsv").open("w", encoding="utf-8", newline="") as dst):
        reader = csv.DictReader(src, delimiter="\t")
        writer = csv.DictWriter(dst, reader.fieldnames + ["raw_file", "scan_number", "event_id"], delimiter="\t")
        writer.writeheader()
        for r in reader:
            name = names[int(r["matching_raw_file_index"])]
            scan = int(r["modification_scan_number"])
            event = "|".join([r["base_peptide_sequence"].replace("L", "I"), r["peptide_position_1based"],
                              r["candidate_peptide_sequence"].replace("L", "I")])
            event_id = "EV_" + hashlib.sha256(event.encode()).hexdigest()[:16]
            spectra.add((name, scan))
            events.add(event_id)
            writer.writerow(dict(r, raw_file=name, scan_number=scan, event_id=event_id))
    return {"strict_physical_scans": len(spectra), "strict_IL_events": len(events)}
