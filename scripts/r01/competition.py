"""Build a Trypsin/P competition database from direct DP candidates and protein references."""
import bisect
import csv
import gzip
from collections import defaultdict
from pathlib import Path


def records(path: Path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig") as handle:
        header, parts = None, []
        for line in handle:
            line = line.strip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(parts)
                header, parts = line[1:], []
            elif line:
                parts.append(line)
        if header is not None:
            yield header, "".join(parts)


def tryptic(parent, start, peptide):
    end = start + len(peptide)
    return ((start == 0 or parent[start - 1] in "KR")
            and (end == len(parent) or peptide[-1] in "KR")
            and sum(a in "KR" for a in peptide[:-1]) <= 2)


def positions(text, fragment):
    pos = text.find(fragment)
    while pos >= 0:
        yield pos
        pos = text.find(fragment, pos + 1)


def emit(handle, header, sequence):
    handle.write(">" + header + "\n")
    for i in range(0, len(sequence), 60):
        handle.write(sequence[i:i+60] + "\n")


def build(candidates, reference, isoforms, contaminants, output):
    # Mass and cleavage profile of the frozen unlabelled DDA competition search.
    from discovery import AA_MASS
    proteins, contexts, normal_id = [], defaultdict(list), {}
    first_source, contaminant_proteins = [], []
    sources = [("original_search", reference), ("human_reference_with_isoforms", isoforms),
               ("contaminant_library", contaminants)]
    with (output / "reference_aliases.tsv").open("w", encoding="utf-8", newline="") as f:
        aliases = csv.writer(f, delimiter="\t")
        aliases.writerow(["source", "accession", "fasta_id", "original_header"])
        for source, path in sources:
            seen = set()
            for head, seq in records(path):
                token = head.split()[0]
                acc = token.split("|")[1] if token.startswith(("sp|", "tr|")) else token
                if acc in seen or not seq or not set(seq) <= set(AA_MASS) | {"U", "X"}:
                    raise ValueError(f"Invalid or duplicate protein: {source}/{acc}")
                seen.add(acc)
                if source == "contaminant_library":
                    ident = acc
                    contaminant_proteins.append((acc, seq))
                else:
                    if seq not in normal_id:
                        normal_id[seq] = f"SFFREF{len(proteins)+1:07d}"
                        proteins.append(seq)
                        first_source.append((source, acc))
                    ident = normal_id[seq]
                contexts[acc].append((source, seq))
                aliases.writerow([source, acc, ident, head])
            if not seen:
                raise ValueError(f"Empty reference: {path}")
    starts, offset = [], 0
    for seq in proteins:
        starts.append(offset)
        offset += len(seq) + 1
    folded = "*".join(proteins).replace("L", "I")
    with candidates.open(encoding="utf-8-sig", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        required = {"candidate_id", "candidate_peptide_sequence", "base_peptide_sequence", "protein_ids",
                    "raw_file", "scan_number", "event_id", "peptide_position_1based", "origin_aa", "destination_aa"}
        if required.difference(reader.fieldnames or []):
            raise ValueError("Expected the discovery candidates.tsv schema.")
        rows = list(reader)
    hits, variants, routes = {}, {}, {}
    unique_ids, tested = set(), set()
    with ((output / "competition.fasta").open("w", encoding="utf-8", newline="\n") as fasta,
          (output / "candidate_context_mapping.tsv").open("w", encoding="utf-8", newline="") as f,
          (output / "normal_competitors.tsv").open("w", encoding="utf-8", newline="") as h):
        mapping = csv.writer(f, delimiter="\t")
        mapping.writerow(["candidate_id", "raw_file", "scan_number", "event_id", "candidate_peptide",
                          "route", "source", "source_accession", "base_start_1based", "site_1based",
                          "variant_fasta_id", "normal_digestible_ids", "variant_context_digestible", "coverage"])
        normal = csv.writer(h, delimiter="\t")
        normal.writerow(["peptide", "protein_id", "start_1based", "match_type", "digestible", "source"])
        for seq, (source, acc) in zip(proteins, first_source):
            emit(fasta, f"sff|{normal_id[seq]}|reference source={source} accession={acc}", seq)
        for r in rows:
            peptide, base = r["candidate_peptide_sequence"], r["base_peptide_sequence"]
            site = int(r["peptide_position_1based"]) - 1
            if (r["candidate_id"] in unique_ids or not 0 <= site < len(base) or len(peptide) != len(base)
                    or sum(a != b for a, b in zip(base, peptide)) != 1 or base[site] != r["origin_aa"]
                    or peptide[site] != r["destination_aa"] or not set(peptide) <= set(AA_MASS)):
                raise ValueError(f"Invalid substitution or duplicate ID: {r['candidate_id']}")
            unique_ids.add(r["candidate_id"])
            if peptide not in hits:
                ids, types = set(), set()
                for acc, sequence in contaminant_proteins:
                    for start in positions(sequence.replace("L", "I"), peptide.replace("L", "I")):
                        match = "exact" if sequence[start:start+len(peptide)] == peptide else "IL_equivalent"
                        normal.writerow([peptide, acc, start+1, match, tryptic(sequence, start, peptide),
                                         "contaminant_library"])
                for pos in positions(folded, peptide.replace("L", "I")):
                    ix = bisect.bisect_right(starts, pos) - 1
                    parent, start = proteins[ix], pos - starts[ix]
                    exact = parent[start:start+len(peptide)] == peptide
                    match = "exact" if exact else "IL_equivalent"
                    digestible = tryptic(parent, start, peptide)
                    normal.writerow([peptide, normal_id[parent], start+1, match, digestible, "reference"])
                    if digestible:
                        ids.add(normal_id[parent])
                        types.add(match)
                hits[peptide] = ";".join(sorted(ids))
                routes[peptide] = ("normal_exact_digestible" if "exact" in types else
                                   "normal_IL_equivalent_digestible" if ids else "candidate_protein_context_required")
            prefix = [r["candidate_id"], r["raw_file"], r["scan_number"], r["event_id"], peptide, routes[peptide]]
            in_range = len(peptide) >= 7 and 18.010565 + sum(AA_MASS[a] for a in peptide) <= 5000
            if hits[peptide] and in_range:
                tested.add(peptide)
            found = False
            for accession in sorted(set(r["protein_ids"].split(";"))):
                acc = accession.removeprefix("CON__").removeprefix("REV__")
                for source, parent in sorted(contexts.get(acc, [])):
                    for start in positions(parent, base):
                        found = True
                        pos = start + site
                        mutant = parent[:pos] + peptide[site] + parent[pos+1:]
                        digestible = tryptic(mutant, start, peptide)
                        ident = ""
                        if not hits[peptide]:
                            ident = normal_id.get(mutant) or variants.get(mutant)
                            if not ident:
                                ident = f"SFFVAR{len(variants)+1:06d}"
                                variants[mutant] = ident
                                emit(fasta, f"sff|{ident}|SAAP_HYPOTHESIS source={source} parent={acc} "
                                     f"substitution={base[site]}{pos+1}{peptide[site]}", mutant)
                        if digestible and in_range:
                            tested.add(peptide)
                        mapping.writerow(prefix + [source, acc, start+1, pos+1, ident, hits[peptide], digestible, "mapped"])
            if not found:
                mapping.writerow(prefix + ["", "", "", "", "", hits[peptide], "", "unresolved_parent"])
    return {"strict_hypotheses": len(rows), "unique_candidates": len(hits), "normal_proteins": len(proteins),
            "variant_proteins": len(variants), "testable_sequences": len(tested),
            "unresolved_or_untestable_sequences": sorted(hits.keys() - tested)}
