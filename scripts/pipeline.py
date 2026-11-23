"""Run the two-search human DDA workflow, also accepting existing MaxQuant exports."""
import argparse
import json
import os
import subprocess
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

from search_chemistry import Chemistry
from saap import sha

HERE=Path(__file__).resolve().parent


def tool(name,*arguments):
    subprocess.run([sys.executable,"-B",str(HERE/name),*map(str,arguments)],check=True)


def freeze(stage, inputs, out, job):
    files={k:Path(v).resolve() for k,v in inputs.items()}
    mounts=[]
    try:
        root=Path(os.path.commonpath([str(p.parent) for p in files.values()]))
    except ValueError:
        root=job.parent
        for k,p in files.items():mounts.extend(["--input",k+"="+str(p)])
    switches=[]
    for k,p in files.items():switches.extend(["--"+k,Path(k)/p.name if mounts else p.relative_to(root)])
    tool("saap.py","freeze","--stage",stage,"--root",root,"--job",job,*switches,*mounts)
    tool("saap.py","run","--root",root,"--job",job,"--out",out,*mounts)
    tool("saap.py","verify","--job",job,"--out",out)


def search(phase,config,design,params,existing,execute):
    if existing:
        txt=existing.resolve()
    else:
        if not execute:
            raise ValueError(f"{phase} search needed: provide existing exports or request search execution.")
        paths=[Path(f[phase+"_path"]) for f in design["files"]]
        if any(not p.is_absolute() for p in paths) or len({p.parent for p in paths})!=1:
            raise ValueError("Each search needs absolute RAW paths in its own run directory.")
        run=paths[0].parent
        if (run/"combined").exists():
            raise FileExistsError("Existing MaxQuant combined output: resume from its verified txt or use a fresh run.")
        for f,path in zip(design["files"],paths):
            source=Path(f["raw_source"]).resolve(strict=True)
            if source.name!=path.name:
                raise ValueError("Retain the original RAW filename")
            path.parent.mkdir(parents=True,exist_ok=True)
            if path.exists():
                if not path.samefile(source):
                    raise ValueError("Existing staged path differs from the declared original: "+str(path))
            else:
                try:os.link(source,path)
                except OSError:path.symlink_to(source)
        command=config["maxquant_command"]
        if not isinstance(command,list) or not command:
            raise ValueError("Expected a command argument list, e.g. dotnet and MaxQuantCmd.dll")
        if not Path(design[phase+"_fasta"]).is_file():
            raise FileNotFoundError(design[phase+"_fasta"])
        with (params.parent/(phase+"_MaxQuant.log")).open("x",encoding="utf-8") as log:
            subprocess.run(command+[str(params.resolve())],stdout=log,stderr=subprocess.STDOUT,check=True)
        txt=run/"combined/txt"
    for name in (("mqpar.xml","msms.txt","msmsScans.txt") if phase=="dp" else ("mqpar.xml","msms.txt")):
        if not (txt/name).is_file():raise FileNotFoundError(txt/name)
    actual=ET.parse(txt/"mqpar.xml").getroot()
    if (actual.findtext("dependentPeptides","").lower()=="true")!=(phase=="dp"):
        raise ValueError("Existing search DP setting conflicts with the selected phase")
    return txt


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config",type=Path)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--dp-txt",type=Path)
    ap.add_argument("--closed-txt",type=Path)
    ap.add_argument("--execute-searches",action="store_true")
    ap.add_argument("--parameters-only",action="store_true")
    ap.add_argument("--predict",action="store_true",help="Run requested inference against the explicitly configured Koina endpoint")
    args=ap.parse_args()
    cfg=json.loads(args.config.read_text(encoding="utf-8-sig"))
    design=json.loads(Path(cfg["design"]).read_text(encoding="utf-8-sig"))
    out=args.out.resolve()
    out.mkdir(parents=True,exist_ok=False)
    design["closed_fasta"]=str(out/"database/competition.fasta")
    design_file=out/"design.json"
    design_file.write_text(json.dumps(design,ensure_ascii=False,indent=2),encoding="utf-8")
    params=out/"parameters"
    tool("mq_design.py",design_file,"--out",params)
    if args.parameters_only:return
    if not args.dp_txt and sha(Path(cfg["reference"]))!=sha(Path(design["dp_fasta"])):
        raise ValueError("Original search reference differs from the reference supplied to competition.")
    planned=Chemistry.from_xml(params/"mqpar_dp.xml")
    dp=search("dp",cfg,design,params/"mqpar_dp.xml",args.dp_txt,args.execute_searches)
    dp_xml=ET.parse(dp/"mqpar.xml").getroot()
    if Chemistry.from_xml(dp/"mqpar.xml")!=planned:
        raise ValueError("Actual DP search chemistry differs from the resolved design")
    discovery=out/"discovery"
    freeze("discover",dict(mqpar=dp/"mqpar.xml",scans=dp/"msmsScans.txt"),discovery,out/"discover.job.json")
    db=out/"database"
    refs={k:cfg[k] for k in ("reference","isoforms","contaminants")}
    freeze("competition",dict(candidates=discovery/"candidates.tsv",mqpar=dp/"mqpar.xml",**refs),db,out/"competition.job.json")
    closed=None
    source=""
    db_match="no_strict_candidates"
    if args.closed_txt or json.loads((discovery/"completion.json").read_text())["result"]["strict_hypotheses"]:
        closed=search("closed",cfg,design,params/"mqpar_closed.xml",args.closed_txt,args.execute_searches)
        cr=ET.parse(closed/"mqpar.xml").getroot()
        if Chemistry.from_xml(closed/"mqpar.xml")!=planned:
            raise ValueError("Closed search chemistry differs from DP search")
        dp_raw=[Path(x.text.replace("\\","/")).name for x in dp_xml.findall("filePaths/string")]
        closed_raw=[Path(x.text.replace("\\","/")).name for x in cr.findall("filePaths/string")]
        if sorted(dp_raw)!=sorted(closed_raw):
            raise ValueError("The two exports have different RAW files")
        source=cr.findtext("fastaFiles/FastaFileInfo/fastaFilePath")
        db_match=Path(source).is_file() and sha(Path(source))==sha(db/"competition.fasta")
    options=dict(cfg.get("evidence",{}))
    options.update(discovery=str(discovery),scans=str(dp/"msmsScans.txt"),base_psms=str(dp/"msms.txt"),
                   closed_psms=str(closed/"msms.txt") if closed else None,context_mapping=str(db/"candidate_context_mapping.tsv"),
                   digest_references=refs,pif_calculated=dp_xml.findtext("calcPeakProperties","").lower()=="true",
                   closed_export_fasta=source,closed_export_fasta_matches_generated=db_match)
    evidence=out/"evidence_inputs.json"
    if args.predict:
        endpoint=options["endpoint"]
        if not endpoint:
            raise ValueError("--predict requires the chosen Koina endpoint")
        if not options.get("mzml") and not options.get("peaks"):
            raise ValueError("Provide observed spectra for model calibration and competition")
        model=options.get("model","Prosit_2024_intensity_PTMs_gl")
        cache=Path(options.get("predictions",out/"predictions")).resolve()
        options["predictions"]=str(cache)
        if not options.get("calibrated_nce"):
            controls=out/"ordinary_controls.tsv"
            tool("calibrate_nce.py","prepare","--relaxed",discovery/"saap_candidates_relaxed.tsv",
                 "--psms",dp/"msms.txt","--mqpar",dp/"mqpar.xml","--model",model,"--out",controls)
            grid=options.get("nce_grid",[24,28,32,36])
            for e in grid:
                tool("prosit_batch.py",controls,"--model",model,"--nce",e,"--cache",cache,"--endpoint",endpoint)
            calibration=out/"calibration_inputs.json"
            calibration.write_text(json.dumps(dict(controls=str(controls),grid=grid,model=model,cache=str(cache),
                                                  peaks=options.get("peaks"),mzml=options.get("mzml",{}))),encoding="utf-8")
            selected=out/"calibration.json"
            tool("calibrate_nce.py","select",calibration,"--out",selected)
            options.update(calibrated_nce=json.loads(selected.read_text())["calibrated_nce"],calibration_source=str(selected))
        evidence.write_text(json.dumps(options,ensure_ascii=False,indent=2),encoding="utf-8")
        request_dir=out/"prediction-inputs"
        tool("saap_table.py",evidence,"--out",request_dir,"--requests-only")
        tool("prosit_batch.py",request_dir/"prediction_requests.tsv","--model",model,
             "--nce",options["calibrated_nce"],"--cache",cache,"--endpoint",endpoint)
    evidence.write_text(json.dumps(options,ensure_ascii=False,indent=2),encoding="utf-8")
    tool("saap_table.py",evidence,"--out",out/"results")
    print(out/"results/SAAP_evidence.tsv")


if __name__=="__main__":
    main()

