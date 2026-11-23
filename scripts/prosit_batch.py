"""Submit eligible requests to a selected Koina-compatible model or export the payloads locally."""
import argparse
import hashlib
import json
import os
import tempfile
import urllib.request
from pathlib import Path

import numpy as np

from fragment_evidence import TOKEN, decode, read

CLASSIC="Prosit_2020_intensity_HCD"
PTM="Prosit_2024_intensity_PTMs_gl"
# Intersection of the two pinned PTM preprocessors, limited to masses handled by fragment_evidence.
VOCAB=set("C_34 C_35 C_4 D_34 E_27 E_34 H_21 H_34 H_35 I_34 K_1 K_34 K_35 K_36 K_37 K_4 L_34 M_35 N_34 N_7 P_21 P_35 Q_28 Q_34 Q_7 R_34 R_36 R_7 S_21 T_21 W_35 Y_21 _1".split())


def eligible(seq,z,model):
    try:
        _,vals=decode(seq)
    except (KeyError,ValueError):
        return False
    if not 7<=len(vals)<=30 or not 1<=z<=6:
        return False
    nt=""
    if seq.startswith("["):
        nt,seq=seq.split("]-",1)
        nt="_"+nt.split(":")[1]
    pairs=TOKEN.findall(seq)
    vocab=VOCAB|{"_737","K_737"} if model==PTM else {"C_4","M_35"}
    return ((not nt or nt in vocab) and all((a!="C" or u=="4") and
            (not u or a+"_"+u in vocab) for a,u in pairs))


def save(path,obj):
    data=json.dumps(obj,allow_nan=False)
    with tempfile.NamedTemporaryFile(mode="w",dir=path.parent,encoding="utf-8",delete=False) as f:
        f.write(data)
    os.replace(f.name,path)


def validate(body,request,model,n):
    if body.get("id")!=request["id"] or body.get("model_name")!=model or str(body.get("model_version"))!="1":
        raise ValueError("Model version or request ID mismatch")
    out={x["name"]:x for x in body["outputs"]}
    shapes=[]
    for name in ("mz","intensities","annotation"):
        x=out[name]
        if x["shape"][0]!=n or len(x["data"])!=int(np.prod(x["shape"])):
            raise ValueError("Unexpected tensor: "+name)
        if name!="annotation" and not np.all(np.isfinite(x["data"])):
            raise ValueError("Nonfinite prediction")
        shapes.append(x["shape"])
    if shapes.count(shapes[0])!=len(shapes):
        raise ValueError("Output tensor shapes differ")


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("requests",type=Path)
    ap.add_argument("--model",choices=[CLASSIC,PTM],required=True)
    ap.add_argument("--nce",type=float,required=True)
    ap.add_argument("--cache",type=Path,required=True)
    ap.add_argument("--endpoint",help="Koina server root, including scheme. Omit to export requests only.")
    args=ap.parse_args()
    if not 0<args.nce<=100:
        raise ValueError("Invalid NCE")
    args.cache.mkdir(parents=True,exist_ok=True)
    pairs={(r["peptide_sequences"],int(r["precursor_charges"])) for r in read(args.requests)}
    selected=sorted((f,z) for f,z in pairs if eligible(f,z,args.model))
    for start in range(0,len(selected),100):
        subset=selected[start:start+100]
        n=len(subset)
        values=[("peptide_sequences","BYTES",[f for f,z in subset]),("precursor_charges","INT32",[z for f,z in subset]),
                ("collision_energies","FP32",[args.nce]*n)]
        if args.model==PTM:
            values.append(("fragmentation_types","BYTES",["HCD"]*n))
        inputs=[dict(name=k,shape=[n,1],datatype=d,data=v) for k,d,v in values]
        ident=hashlib.sha256(json.dumps(inputs,sort_keys=True).encode()).hexdigest()[:24]
        request={"id":ident,"inputs":inputs}
        dest=args.cache/(args.model+"_"+ident+".response.json")
        mapping=[dict(peptide_sequences=f,precursor_charges=z,effective_CE=args.nce) for f,z in subset]
        save(dest.with_name(dest.name.replace(".response.",".request.")),request)
        save(dest.with_name(dest.name.replace(".response.",".rows.")),mapping)
        if dest.exists():
            validate(json.loads(dest.read_text(encoding="utf-8")),request,args.model,n)
        elif args.endpoint:
            url=args.endpoint.rstrip("/")+"/v2/models/"+args.model+"/versions/1/infer"
            req=urllib.request.Request(url,data=json.dumps(request).encode(),headers={"Content-Type":"application/json"})
            with urllib.request.urlopen(req,timeout=120) as response:
                body=json.load(response)
            validate(body,request,args.model,n)
            save(dest,body)
    print(json.dumps({"requested":len(pairs),"model_eligible":len(selected),"unsupported":len(pairs)-len(selected),
                      "mode":"inference" if args.endpoint else "export_only"}))


if __name__=="__main__":
    main()

