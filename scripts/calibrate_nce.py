"""Select a model-specific collision energy using ordinary PSMs and sequence-disjoint holdout."""
import argparse
import hashlib
import json
from pathlib import Path
from statistics import median

from fragment_evidence import angle, encode, mq_form, mzml_peaks, native_cache, read, spectrum
from prosit_batch import eligible
from saap_table import tsv
from search_chemistry import Chemistry


def prepare(args):
    chemistry=Chemistry.from_xml(args.mqpar) if args.mqpar else Chemistry()
    excluded=set()
    for r in read(args.relaxed):
        excluded.update(r[k].replace("I","L") for k in ("base_peptide_sequence","candidate_peptide_sequence"))
    controls={}
    for r in read(args.psms):
        seq=r["Sequence"]
        if (seq.replace("I","L") in excluded or r["Decoy"]=="+" or r["Potential contaminant"]=="+"
                or not r["PEP"] or not float(r["PEP"])<=.01 or r["Fragmentation"]!="HCD"):
            continue
        try:
            form=encode(seq,mq_form(seq,r["Modified sequence"]),chemistry)
            z=int(r["Charge"])
        except (ValueError,KeyError):
            continue
        if not eligible(form,z,args.model):
            continue
        token=hashlib.sha256(seq.replace("I","L").encode()).hexdigest()
        split="holdout" if int(token[:8],16)%3==0 else "calibration"
        record=dict(raw_file=r["Raw file"],scan_number=r["Scan number"],peptide_sequences=form,
                    precursor_charges=z,split=split,PEP=float(r["PEP"]))
        # One spectrum per sequence per raw limits influence from repeated measurements.
        key=(token,r["Raw file"])
        if key not in controls or record["PEP"]<controls[key]["PEP"]:
            controls[key]=record
    if args.per_split<1:
        raise ValueError("Invalid panel size")
    panel=[]
    for split in ("calibration","holdout"):
        pool=[r for k,r in sorted(controls.items()) if r["split"]==split]
        if not pool:
            raise ValueError("No eligible ordinary PSMs in "+split)
        panel.extend(pool[:args.per_split])
    tsv(args.out,list(panel[0]),panel)


def select(args):
    config=json.loads(args.config.read_text(encoding="utf-8-sig"))
    rows=list(read(config["controls"]))
    keys={r["raw_file"]+":scan="+r["scan_number"] for r in rows}
    data={}
    if config.get("peaks"):
        for k,v in json.loads(Path(config["peaks"]).read_text(encoding="utf-8-sig")).items():
            if k in keys:data[k]=v
    for raw,path in config.get("mzml",{}).items():
        wanted={r["scan_number"] for r in rows if r["raw_file"]==raw}
        data.update({raw+":scan="+k:v for k,v in mzml_peaks(Path(path),wanted).items()})
    splits={}
    for r in rows:
        import re
        seq=re.sub(r"\[UNIMOD:\d+\]-?","",r["peptide_sequences"]).replace("I","L")
        if splits.setdefault(seq,r["split"])!=r["split"]:
            raise ValueError("Sequence leakage between calibration and holdout")
    pairs={(r["peptide_sequences"],int(r["precursor_charges"])) for r in rows}
    energies=sorted(set(map(float,config["grid"])))
    if len(energies)<2 or not all(0<e<=100 for e in energies):
        raise ValueError("Compare at least two distinct, valid NCE values")
    cache={e:native_cache(Path(config["cache"]),pairs,config["model"],e) for e in energies}
    scores={e:{"calibration":[],"holdout":[]} for e in energies}
    for r in rows:
        z=int(r["precursor_charges"])
        key=(r["peptide_sequences"],z)
        d=data.get(r["raw_file"]+":scan="+r["scan_number"])
        if not d or d.get("charge")!=z or not all(key in c for c in cache.values()):
            continue
        values={e:angle(cache[e][key],spectrum(d),20) for e in energies}
        if all(v is not None for v in values.values()):
            for e,v in values.items():scores[e][r["split"]].append(v)
    if not all(scores[energies[0]].values()):
        raise ValueError("Both sequence partitions need spectra scored across the same NCE grid")
    best=max(energies,key=lambda e:median(scores[e]["calibration"]))
    result={"model":config["model"],"calibrated_nce":best,"common_controls":{s:len(v) for s,v in scores[best].items()},
            "median_calibration_SA":{str(e):median(scores[e]["calibration"]) for e in energies},
            "median_holdout_SA_at_selected_NCE":median(scores[best]["holdout"])}
    with args.out.open("x",encoding="utf-8") as f:json.dump(result,f,indent=2)


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    sub=ap.add_subparsers(required=True)
    p=sub.add_parser("prepare")
    p.add_argument("--relaxed",required=True,type=Path)
    p.add_argument("--psms",required=True,type=Path)
    p.add_argument("--mqpar",type=Path)
    p.add_argument("--model",required=True)
    p.add_argument("--per-split",type=int,default=500)
    p.set_defaults(fn=prepare)
    q=sub.add_parser("select")
    q.add_argument("config",type=Path)
    q.set_defaults(fn=select)
    for command in (p,q):command.add_argument("--out",required=True,type=Path)
    args=ap.parse_args()
    args.fn(args)

