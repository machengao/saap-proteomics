"""Inventory local method evidence for design resolution before reading peak arrays."""
import argparse
import csv
import json
import re
import xml.etree.ElementTree as ET
from itertools import islice
from pathlib import Path

METHOD = re.compile(r"acquisition|SPS.?MS3|DDA|DIA|HCD|CID|ETD|Orbitrap|Exactive|TMT|iTRAQ|SILAC|"
                    r"label.free|MaxLFQ|trypsin|Lys.?C|enzyme|alkylat|iodoacetamide|carbamido|oxidation|"
                    r"deamidation|contaminant|定量|酶切|漏切|修饰|烷基|质谱|采集|通道|污染",re.I)


def inspect(root):
    report={"raw_files":[],"methods":[],"tables":[],"search_parameters":[],"mzml":[]}
    for p in sorted(root.rglob("*")):
        ext=p.suffix.lower()
        if ext in (".raw",".d"):
            report["raw_files"].append(str(p))
        elif p.name.lower()=="mqpar.xml":
            doc=ET.parse(p).getroot()
            groups=[]
            for g in doc.findall("parameterGroups/parameterGroup"):
                groups.append({key:[e.text for e in g.findall(key+"/string")] if g.findall(key+"/string")
                               else g.findtext(key) for key in ("msInstrument","lcmsRunType","lfqMode","multiplicity",
                                    "enzymes","maxMissedCleavages","fixedModifications","variableModifications")})
            report["search_parameters"].append({"path":str(p),"version":doc.findtext("maxQuantVersion"),
                                                "contaminants":doc.findtext("includeContaminants"),"groups":groups})
        elif p.name=="msmsScans.txt":
            with p.open(encoding="utf-8-sig",newline="") as f:
                rows=list(islice(csv.DictReader(f,delimiter="\t"),128))
            report["tables"].append({"path":str(p),"sampled":len(rows),"signals":{
                k:sorted({r.get(k,"") for r in rows}) for k in ("Raw file","Fragmentation","Mass analyzer","Type")}})
        elif ext==".mzml":
            report["mzml"].append(str(p))
        elif ext in (".pdf",".docx",".xlsx") or (ext in (".tsv",".csv") and re.search(r"design|metadata|sample|channel|样本|分组|通道",p.name,re.I)):
            report["methods"].append({"path":str(p),"read_if_relevant":True})
        elif p.is_file() and ext in (".txt",".md") and p.name not in ("msms.txt","peptides.txt","evidence.txt","dependentPeptides.txt") and p.stat().st_size<=2_000_000:
            lines=p.read_text(encoding="utf-8-sig",errors="replace").splitlines()
            matches=[{"line":i,"text":s[:600]} for i,s in enumerate(lines,1) if METHOD.search(s)]
            if matches:report["methods"].append({"path":str(p),"matches":matches[:24],"total_matches":len(matches)})
    return report


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("folder",type=Path)
    args=ap.parse_args()
    print(json.dumps(inspect(args.folder.resolve(strict=True)),ensure_ascii=False,indent=2))

