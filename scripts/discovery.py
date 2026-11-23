#!/usr/bin/env python3
"""Scan-level SAAP mass enumeration adapted from the project MaxQuant discovery script."""
from __future__ import annotations

import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd


# Monoisotopic residue masses (unmodified residues within a peptide), Da.
AA_MASS = {
    "A": 71.037114, "R": 156.101111, "N": 114.042927, "D": 115.026943,
    # This workflow's MaxQuant search used fixed Carbamidomethyl (C), so the
    # effective cysteine residue mass includes +57.021464 Da, as in decode.
    "C": 160.030649, "E": 129.042593, "Q": 128.058578, "G": 57.021464,
    "H": 137.058912, "I": 113.084064, "L": 113.084064, "K": 128.094963,
    "M": 131.040485, "F": 147.068414, "P": 97.052764, "S": 87.032028,
    "T": 101.047679, "W": 186.079313, "Y": 163.063329, "V": 99.068414,
}

REQUIRED_COLUMNS = {
    "Base peptide sequence", "Modification", "Cluster index", "Cluster mass",
    "Cluster size total", "Cluster size forward", "Cluster size reverse",
    "Base peptide raw file", "Matches in Raw files", "Localization probability",
    "Max Probability", "Posterior error probability", "Score", "Mass difference",
    "Peptide Length Difference", "Protein IDs", "Gene names", "Protein names",
    "Modification scan numbers",
}

KNOWN_OR_LIKELY_NON_SUBSTITUTION = re.compile(
    r"oxidation|deamidation|loss water|loss ammonia|cation:|13c-12c|"
    r"carbamid|acetyl|methyl|phospho|ubiquitin|glyc|formyl|carbamyl|"
    r"pyro|adduct|isotope",
    flags=re.IGNORECASE,
)


def split_semicolon(value: object) -> list[str]:
    if pd.isna(value):
        return []
    return [x.strip() for x in str(value).split(";")]


def safe_float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def aligned_value(value: object, index: int) -> str:
    values = split_semicolon(value)
    if not values:
        return ""
    if index < len(values):
        return values[index]
    if len(values) == 1:
        return values[0]
    return ""


def parse_localization(text: object, base_sequence: str) -> list[tuple[int, str, float]]:
    """Return (0-based position, residue, probability) from MQ localization text."""
    if pd.isna(text):
        return []
    text = str(text).replace("_", "")
    tokens = re.findall(r"([A-Z])(?:\(([0-9.eE+\-]+)\))?", text)
    parsed: list[tuple[int, str, float]] = []
    for pos, (aa, prob) in enumerate(tokens):
        if pos >= len(base_sequence):
            break
        # Ignore malformed localization strings whose residue order differs.
        if base_sequence[pos] != aa:
            continue
        parsed.append((pos, aa, safe_float(prob) if prob else 0.0))
    return parsed


def substitution_table(masses=None) -> pd.DataFrame:
    masses = AA_MASS if masses is None else masses
    rows = []
    for old, old_mass in masses.items():
        for new, new_mass in masses.items():
            if old == new:
                continue
            if {old, new} == {"I", "L"}:
                continue  # indistinguishable by ordinary MS mass
            rows.append({
                "origin_aa": old,
                "destination_aa": new,
                "substitution": f"{old}>{new}",
                "theoretical_delta_da": new_mass - old_mass,
                "destination_is_IL_ambiguous": new in {"I", "L"},
            })
    return pd.DataFrame(rows).sort_values(
        ["origin_aa", "theoretical_delta_da", "destination_aa"]
    ).reset_index(drop=True)


def modification_class(value: object) -> str:
    if pd.isna(value) or not str(value).strip():
        return "unnamed"
    value = str(value).strip()
    if value.lower().startswith("diffuse"):
        return "diffuse_localization"
    if value.lower() == "unmodified":
        return "unmodified_annotation"
    if value.lower().startswith("comp.:"):
        return "composition_only"
    if KNOWN_OR_LIKELY_NON_SUBSTITUTION.search(value):
        return "known_or_likely_ptm_artifact"
    return "other_named_change"


def unique_fasta(candidates: pd.DataFrame, path: Path) -> int:
    unique = candidates.drop_duplicates("candidate_peptide_sequence").copy()
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        for i, row in enumerate(unique.itertuples(index=False), start=1):
            seq = row.candidate_peptide_sequence
            header = (
                f">SAAP_CAND_{i:06d}|{row.substitution}|pep_pos_{row.peptide_position_1based}"
                f"|cluster_{row.cluster_index}|source_DP_row_{row.dp_row_number}"
            )
            handle.write(header + "\n" + seq + "\n")
    return len(unique)


def discover(input_path: Path, output_dir: Path, mass_tolerance_da: float, pep_max: float, localization_min: float, masses=None):
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.exists():
        raise FileNotFoundError(input_path)
    dp = pd.read_csv(input_path, sep="\t", low_memory=False, keep_default_na=False, na_values=[""], dtype={"Base peptide sequence": str})
    missing = sorted(REQUIRED_COLUMNS - set(dp.columns))
    if missing:
        raise ValueError("Missing required columns: " + ", ".join(missing))

    subs = substitution_table(masses)
    subs.to_csv(output_dir / "aa_substitution_mass_table.tsv", sep="\t", index=False)
    subs_by_origin = {
        aa: list(subs[subs["origin_aa"] == aa].itertuples(index=False))
        for aa in AA_MASS
    }

    candidate_rows: list[dict[str, object]] = []
    skipped = Counter()

    for table_index, row in dp.iterrows():
        base = str(row["Base peptide sequence"]).strip()
        if not base or base == "nan" or any(aa not in AA_MASS for aa in base):
            skipped["invalid_base_sequence"] += 1
            continue

        localizations = parse_localization(row["Localization probability"], base)
        if not localizations:
            skipped["no_parseable_localization"] += 1
            continue

        mass_values = split_semicolon(row["Mass difference"])
        if not mass_values:
            skipped["no_mass_difference"] += 1
            continue

        for obs_index, mass_text in enumerate(mass_values):
            observed = safe_float(mass_text)
            if not np.isfinite(observed):
                skipped["non_numeric_mass_difference"] += 1
                continue

            length_difference = safe_float(
                aligned_value(row["Peptide Length Difference"], obs_index)
            )
            scan_number = aligned_value(row["Modification scan numbers"], obs_index)
            feature_id = aligned_value(row.get("Feature Ids", ""), obs_index)
            raw_match_index = aligned_value(row["Matches in Raw files"], obs_index)
            mod_intensity = safe_float(
                aligned_value(row.get("Modification intensities", ""), obs_index)
            )

            for pos, old_aa, loc_prob in localizations:
                matches = (
                    match for match in subs_by_origin[old_aa]
                    if abs(match.theoretical_delta_da - observed)
                    <= mass_tolerance_da
                )
                for match in matches:
                    new_aa = match.destination_aa
                    candidate_sequence = base[:pos] + new_aa + base[pos + 1:]
                    mod_class = modification_class(row["Modification"])
                    mass_error_da = observed - match.theoretical_delta_da
                    candidate_rows.append({
                        "candidate_id": "",  # populated after sorting
                        "dp_row_number": int(table_index) + 2,
                        "observation_index": obs_index,
                        "base_peptide_sequence": base,
                        "candidate_peptide_sequence": candidate_sequence,
                        "substitution": match.substitution,
                        "origin_aa": old_aa,
                        "destination_aa": new_aa,
                        "peptide_position_1based": pos + 1,
                        "localization_probability_at_site": loc_prob,
                        "max_probability": safe_float(row["Max Probability"]),
                        "posterior_error_probability": safe_float(
                            row["Posterior error probability"]
                        ),
                        "dp_score": safe_float(row["Score"]),
                        "observed_delta_da": observed,
                        "theoretical_delta_da": match.theoretical_delta_da,
                        "mass_error_da": mass_error_da,
                        "mass_error_mda": mass_error_da * 1000.0,
                        "mass_tolerance_da_used": mass_tolerance_da,
                        "peptide_length_difference": length_difference,
                        "modification_annotation": row["Modification"],
                        "modification_class": mod_class,
                        "ptm_or_artifact_warning": mod_class
                        in {"known_or_likely_ptm_artifact", "diffuse_localization"},
                        "cluster_index": row["Cluster index"],
                        "cluster_mass": row["Cluster mass"],
                        "cluster_mass_sd": row.get("Cluster mass SD", np.nan),
                        "cluster_size_total": row["Cluster size total"],
                        "cluster_size_forward": row["Cluster size forward"],
                        "cluster_size_reverse": row["Cluster size reverse"],
                        "base_peptide_raw_file_indices": row["Base peptide raw file"],
                        "matching_raw_file_index": raw_match_index,
                        "feature_id": feature_id,
                        "modification_scan_number": scan_number,
                        "found_by": row.get("Found by", ""),
                        "modification_intensity_observation": mod_intensity,
                        "protein_ids": row["Protein IDs"],
                        "gene_names": row["Gene names"],
                        "protein_names": row["Protein names"],
                        "base_intensity_KO": row.get("Base intensity Exp: KO", np.nan),
                        "base_intensity_WT": row.get("Base intensity Exp: WT", np.nan),
                        "candidate_intensity_KO": row.get(
                            "Modification intensity Exp: KO", np.nan
                        ),
                        "candidate_intensity_WT": row.get(
                            "Modification intensity Exp: WT", np.nan
                        ),
                        "ratio_candidate_base_KO": row.get(
                            "Ratio mod/base intensity Exp: KO", np.nan
                        ),
                        "ratio_candidate_base_WT": row.get(
                            "Ratio mod/base intensity Exp: WT", np.nan
                        ),
                        "localization_probability_string": row[
                            "Localization probability"
                        ],
                        "aminoacid_localizations": row.get(
                            "Aminoacid localizations", ""
                        ),
                    })

    if not candidate_rows:
        fields = (
            'candidate_id dp_row_number observation_index base_peptide_sequence '
            'candidate_peptide_sequence substitution origin_aa destination_aa '
            'peptide_position_1based localization_probability_at_site max_probability posterior_error_probability '
            'dp_score observed_delta_da theoretical_delta_da mass_error_da '
            'mass_error_mda mass_tolerance_da_used peptide_length_difference modification_annotation '
            'modification_class ptm_or_artifact_warning cluster_index cluster_mass '
            'cluster_mass_sd cluster_size_total cluster_size_forward cluster_size_reverse '
            'base_peptide_raw_file_indices matching_raw_file_index feature_id modification_scan_number '
            'found_by modification_intensity_observation protein_ids gene_names '
            'protein_names base_intensity_KO base_intensity_WT candidate_intensity_KO '
            'candidate_intensity_WT ratio_candidate_base_KO ratio_candidate_base_WT localization_probability_string '
            'aminoacid_localizations passes_strict_discovery_filters '
        ).split()
        for filename in ("saap_candidates_relaxed.tsv", "saap_candidates_strict.tsv"):
            pd.DataFrame(columns=fields).to_csv(output_dir / filename, sep="\t", index=False)
        (output_dir / "saap_candidates_strict.fasta").write_bytes(b"")
        return {"direct_target_dp": len(dp), "relaxed_hypotheses": 0,
                "strict_hypotheses": 0, "strict_sequences": 0, "skipped": dict(skipped)}

    candidates = pd.DataFrame(candidate_rows).sort_values(
        ["posterior_error_probability", "mass_error_da", "dp_row_number"],
        key=lambda s: s.abs() if s.name == "mass_error_da" else s,
    ).reset_index(drop=True)
    candidates["candidate_id"] = [f"SAAP_{i:07d}" for i in range(1, len(candidates) + 1)]

    candidates["passes_strict_discovery_filters"] = (
        (candidates["posterior_error_probability"] <= pep_max)
        & (candidates["max_probability"] >= localization_min)
        & (candidates["localization_probability_at_site"] >= localization_min)
        & (candidates["peptide_length_difference"] == 0)
        & (candidates["cluster_size_forward"] > candidates["cluster_size_reverse"])
        & (candidates["mass_error_da"].abs() <= mass_tolerance_da)
    )
    candidates["passes_strict_discovery_filters"] &= ~candidates["ptm_or_artifact_warning"]

    strict = candidates[candidates["passes_strict_discovery_filters"]].copy()
    strict = strict.drop_duplicates(
        ["base_peptide_sequence", "candidate_peptide_sequence", "cluster_index",
         "matching_raw_file_index", "modification_scan_number"]
    )

    relaxed_path = output_dir / "saap_candidates_relaxed.tsv"
    strict_path = output_dir / "saap_candidates_strict.tsv"
    fasta_path = output_dir / "saap_candidates_strict.fasta"
    candidates.to_csv(relaxed_path, sep="\t", index=False)
    strict.to_csv(strict_path, sep="\t", index=False)
    fasta_entries = unique_fasta(strict, fasta_path)

    return {"direct_target_dp": len(dp), "relaxed_hypotheses": len(candidates),
            "strict_hypotheses": len(strict), "strict_sequences": fasta_entries,
            "skipped": dict(skipped)}
