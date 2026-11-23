"""Freeze or execute a portable, checksum-verified SAAP stage."""
import argparse
import csv
import hashlib
import json
import math
import platform
import sys
import subprocess
from dataclasses import asdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
STAGES = {"discover": ("mqpar", "scans"),
          "competition": ("candidates", "reference", "isoforms", "contaminants")}
CODE = ("saap.py", "scan_input.py", "discovery.py", "competition.py", "search_chemistry.py","../assets/enzyme-pairs-2.8.1.json")


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def code_hashes():
    return {name: sha(HERE / name) for name in CODE}


def dump(path, data):
    with path.open("x", encoding="utf-8") as f:
        json.dump(data, f, ensure_ascii=False, indent=2, allow_nan=False)
        f.write("\n")


def validate_params(params):
    if set(params) != {"mass_tolerance_da", "pep_max", "localization_min"}:
        raise ValueError("Unknown discovery settings")
    if (not all(math.isfinite(x) for x in params.values()) or params["mass_tolerance_da"] <= 0
            or not 0 <= params["pep_max"] <= 1 or not 0 <= params["localization_min"] <= 1):
        raise ValueError("Invalid discovery thresholds")


def bind(root, relative):
    path = (root / relative).resolve()
    if Path(relative).is_absolute() or not path.is_relative_to(root):
        raise ValueError("Input must be relative to the data root: " + relative)
    return path


def mounted_paths(args):
    entries = {}
    for item in args.input:
        key, separator, value = item.partition("=")
        if not separator or key in entries or not value:
            raise ValueError("Use --input role=path once per mounted input")
        entries[key] = Path(value).resolve(strict=True)
    return entries


def freeze(args):
    root = args.root.resolve()
    inputs = {}
    mounted = mounted_paths(args)
    keys=list(STAGES[args.stage])
    if args.stage=="competition" and args.mqpar:
        keys.append("mqpar")
    for key in keys:
        relative = getattr(args, key)
        if not relative:
            raise ValueError("Missing --" + key)
        logical = bind(root, relative)
        path = mounted.get(key, logical)
        inputs[key] = {"path":logical.relative_to(root).as_posix(),"size":path.stat().st_size,"sha256":sha(path)}
    if mounted.keys()-inputs.keys():
        raise ValueError("Unrecognized input role")
    params = {}
    if args.stage == "discover":
        params = {k: getattr(args, k) for k in ("mass_tolerance_da", "pep_max", "localization_min")}
        validate_params(params)
    dump(args.job, {"schema": 1, "stage": args.stage, "inputs": inputs, "parameters": params, "code": code_hashes()})


def execute(args):
    from scan_input import label_strict, prepare, search_profile
    from discovery import discover
    from competition import build
    import numpy
    import pandas
    job = json.loads(args.job.read_text(encoding="utf-8"))
    stage = job["stage"]
    r01=HERE/"r01"
    if job["code"] == {name:sha(r01/name) for name in CODE[:4]}:
        subprocess.run([sys.executable,"-B",str(r01/"saap.py"),"run","--job",str(args.job),"--root",str(args.root),"--out",str(args.out)],check=True)
        return
    keys=set(STAGES[stage])
    if stage=="competition" and "mqpar" in job["inputs"]:
        keys.add("mqpar")
    if job["schema"] != 1 or job["code"] != code_hashes() or set(job["inputs"]) != keys:
        raise ValueError("Job schema, stage inputs, or deployed script versions differ.")
    if stage == "discover":
        validate_params(job["parameters"])
    elif job["parameters"]:
        raise ValueError("Competition uses the documented fixed search profile.")
    paths = {}
    mounted = mounted_paths(args)
    if mounted.keys()-job["inputs"].keys():
        raise ValueError("Unrecognized input role")
    for key, item in job["inputs"].items():
        path = mounted.get(key, bind(args.root.resolve(), item["path"]))
        if path.stat().st_size != item["size"] or sha(path) != item["sha256"]:
            raise ValueError(f"Input checksum mismatch: {key}")
        paths[key] = path
    out = args.out
    out.mkdir(parents=True, exist_ok=False)
    csv.field_size_limit(100_000_000)
    dump(out / "job.json", job)
    try:
        if stage == "discover":
            from search_chemistry import Chemistry
            chemistry=Chemistry.from_xml(paths["mqpar"])
            names = search_profile(paths["mqpar"])
            result = prepare(paths["scans"], names, out)
            result["chemistry"]=asdict(chemistry)
            result.update(discover(out / "dp_scan_input.tsv", out, masses=chemistry.masses, **job["parameters"]))
            result.update(label_strict(out, names))
        else:
            result = build(**paths, output=out)
        report = {"status": "complete", "stage": stage, "result": result,
                  "python": platform.python_version(), "pandas": pandas.__version__, "numpy": numpy.__version__,
                  "artifacts": {p.name: sha(p) for p in sorted(out.iterdir()) if p.is_file()}}
        dump(out / "completion.json", report)
        print(json.dumps(result, ensure_ascii=False))
    except Exception as exc:
        dump(out / "failure.json", {"stage": stage, "error": str(exc), "type": type(exc).__name__})
        raise


def verify(args):
    report = json.loads((args.out / "completion.json").read_text(encoding="utf-8"))
    frozen = json.loads(args.job.read_text(encoding="utf-8"))
    returned = json.loads((args.out / "job.json").read_text(encoding="utf-8"))
    if frozen != returned or report["status"] != "complete" or report["stage"] != frozen["stage"]:
        raise ValueError("Returned run differs from the frozen job or is incomplete.")
    expected = {"job.json"}
    if report["stage"] == "discover":
        expected.update(("saap_candidates_strict.tsv", "saap_candidates_relaxed.tsv", "candidates.tsv",
                         "dp_scan_input.tsv", "dp_review.tsv", "aa_substitution_mass_table.tsv",
                         "saap_candidates_strict.fasta"))
    elif report["stage"] == "competition":
        expected.update(("competition.fasta", "candidate_context_mapping.tsv", "reference_aliases.tsv",
                         "normal_competitors.tsv"))
    else:
        raise ValueError("Unknown returned stage")
    if expected != set(report["artifacts"]):
        raise ValueError("Returned artifact inventory is incomplete or unexpected.")
    for name in expected:
        if sha(args.out / name) != report["artifacts"][name]:
            raise ValueError("Returned artifact checksum mismatch: " + name)
    print("Verified", len(expected), "files for", report["stage"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)
    make = sub.add_parser("freeze", help="Create a portable job using paths relative to --root")
    make.set_defaults(action=freeze)
    make.add_argument("--stage", choices=STAGES, required=True)
    for key in dict.fromkeys(k for stage in STAGES.values() for k in stage):
        make.add_argument("--" + key)
    make.add_argument("--mass-tolerance-da", type=float, default=0.005)
    make.add_argument("--pep-max", type=float, default=0.01)
    make.add_argument("--localization-min", type=float, default=0.90)
    run = sub.add_parser("run", help="Verify inputs and run one stage in a new output directory")
    run.set_defaults(action=execute)
    run.add_argument("--out", required=True, type=Path)
    check = sub.add_parser("verify", help="Check a returned stage against the frozen job and output checksums")
    check.set_defaults(action=verify)
    check.add_argument("--out", required=True, type=Path)
    check.add_argument("--job", required=True, type=Path)
    for command in (make, run):
        command.add_argument("--input",action="append",default=[],help="Bind role=path on an independently mounted volume")
        command.add_argument("--root", required=True, type=Path)
        command.add_argument("--job", required=True, type=Path)
    args = parser.parse_args()
    args.action(args)


if __name__ == "__main__":
    main()
