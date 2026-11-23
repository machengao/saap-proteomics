"""Run the 18-scan exported-table regression on a new host. No RAW inputs."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

def rows(path):
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))

def key(row):
    return tuple(row[k] for k in ("raw_file", "scan_number", "base_peptide_sequence",
                                  "candidate_peptide_sequence", "peptide_position_1based"))

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--skill", required=True, type=Path)
    parser.add_argument("--out", required=True, type=Path)
    args = parser.parse_args()
    source = Path(__file__).resolve().parent
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    def cli(*parts):
        subprocess.run([sys.executable, "-B", str(args.skill / "scripts/saap.py"), *map(str, parts)], check=True)
    d = out / "discovery"
    job = source / "discovery-r2.json"
    cli("run", "--root", source, "--job", job, "--out", d)
    cli("verify", "--job", job, "--out", d)
    current = rows(d / "candidates.tsv")
    expected = {key(r): r for r in rows(source / "expected.tsv")}
    if len(current) != len(expected) or {key(r) for r in current} != set(expected):
        raise ValueError("Different strict candidate set")
    for r in current:
        record = expected[key(r)]
        if r["event_id"] != record["event_id_IL_collapsed"] or not all(
            math.isclose(float(r[k]), float(record[k]), abs_tol=1e-12)
            for k in ("posterior_error_probability", "mass_error_da")
        ):
            raise ValueError(f"Different evidence for {key(r)}")
    archived = json.loads((source / "validation.json").read_text(encoding="utf-8"))
    reference_job = archived["jobs"][1]
    for field in ("reference", "isoforms", "contaminants"):
        name = reference_job["inputs"][field]["path"]
        shutil.copyfile(source / name, out / name)
    comp = out / "competition"
    job = out / "competition.json"
    cli("freeze", "--stage", "competition", "--root", out, "--job", job,
        "--candidates", "discovery/candidates.tsv", "--reference", "ref0.fasta",
        "--isoforms", "ref1.fasta", "--contaminants", "ref2.fasta")
    frozen = json.loads(job.read_text(encoding="utf-8"))
    for field in ("reference", "isoforms", "contaminants"):
        if frozen["inputs"][field] != reference_job["inputs"][field]:
            raise ValueError("Reference extract differs: " + field)
    cli("run", "--root", out, "--job", job, "--out", comp)
    cli("verify", "--job", job, "--out", comp)
    counts = json.loads((comp / "completion.json").read_text(encoding="utf-8"))["result"]
    if (counts["normal_proteins"], counts["variant_proteins"], counts["testable_sequences"],
            len(rows(comp / "candidate_context_mapping.tsv"))) != (79, 6, 8, 72):
        raise ValueError("Different competition coverage")
    print("PASS: 14 strict hypotheses and competition coverage match the frozen subset.")

if __name__ == "__main__":
    main()
