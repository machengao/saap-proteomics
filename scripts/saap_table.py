"""Build one substitution table, retaining every mass-compatible direct DP hypothesis."""
import argparse
import bisect
import csv
import hashlib
import json
import math
from itertools import combinations_with_replacement, product
from collections import defaultdict
from pathlib import Path

import numpy as np
from competition import records
from discovery import AA_MASS
from fragment_evidence import (MODS, SITE_MODS, angle, encode, fragment_score, ions,
                               mq_form, mzml_peaks, native_cache, neutral, read, spectrum)
from saap import sha
from search_chemistry import Chemistry


def packed(value):
    return json.dumps(value, ensure_ascii=False, separators=(",",":"), allow_nan=False)


def identity(row):
    return (row["base_peptide_sequence"], int(row["peptide_position_1based"]), row["candidate_peptide_sequence"])


def tsv(path, fields, rows):
    with path.open("x", encoding="utf-8-sig", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def digest_index(refs, target_masses, chemistry):
    """Index unmodified, fixed-C peptides only near requested precursor masses."""
    found = defaultdict(set)
    targets = sorted(set(target_masses))
    if not targets:
        return found
    for category,path in refs.items():
        for head,seq in records(Path(path)):
            cuts = [0]+[i for i in range(1,len(seq)) if seq[i-1:i+1] in chemistry.pairs]
            if cuts[-1]!=len(seq):
                cuts.append(len(seq))
            for i,start in enumerate(cuts[:-1]):
                for end in cuts[i+1:i+2+chemistry.missed]:
                    peptide = seq[start:end]
                    if len(peptide)<chemistry.min_length or not set(peptide)<=AA_MASS.keys():
                        continue
                    m = chemistry.mass(peptide)
                    if m>chemistry.max_mass:
                        continue
                    tol = max(.005,m*10.001e-6)
                    left,right = bisect.bisect_left(targets,m-tol),bisect.bisect_right(targets,m+tol)
                    for target in targets[left:right]:
                        found[target].add((encode(peptide,chemistry=chemistry),category))
    return found


def prepared(config, wanted):
    """Join by physical scan, base, position and candidate; historical ranking IDs are run-local."""
    units, links = {}, defaultdict(lambda: defaultdict(set))
    if not config.get("units"):
        return units, links
    for r in read(config["units"]):
        key = ((r["raw_file"],r["scan_number"]), (r["base_sequence"],int(r["position_1based"]),r["candidate_sequence"]))
        if key in wanted and r["unit_id"].startswith("DP_"):
            units[r["unit_id"]] = key
    for r in read(config["links"]):
        if r["unit_id"] not in units:
            continue
        key = units[r["unit_id"]]
        if r["prediction_role"]=="candidate":
            links[key][r["unit_id"]].add((r["proforma_sequence"],"candidate"))
        elif all(r[f]=="True" for f in ("mass_window_compatible","precursor_compatible","context_eligible")):
            links[key][r["unit_id"]].add((r["proforma_sequence"],r["relationship"]))
    return units,links


def ptm_forms(sequence, fixed_state, target, chemistry):
    """Bounded same-base PTM competition; test mass before enumerating site combinations."""
    bf=encode(sequence,fixed_state,chemistry)
    tolerance=max(.005,target*1e-5)
    delta=target-neutral(bf)
    allowed={code:[i for i,a in enumerate(sequence,1) if a in residues and i not in fixed_state
                   and not (a=="K" and chemistry.tmt)] for code,residues in SITE_MODS.items()}
    allowed={code:sites for code,sites in allowed.items() if sites}
    for size in (1,2):
        for codes in combinations_with_replacement(allowed,size):
            if abs(sum(MODS[c] for c in codes)-delta)>tolerance:
                continue
            for positions in product(*(allowed[c] for c in codes)):
                if len(set(positions))==size:
                    yield encode(sequence,{**fixed_state,**dict(zip(positions,codes))},chemistry)


def score_unit(state, spectra, pred, ppm):
    z, key, form, alternatives = state["charge"],state["key"],state["form"],state["alternatives"]
    base,site,peptide = state["identity"]
    result = {"candidate_form":form,"eligible_forms":len(alternatives),"coverage":state["coverage"],
              "comparison_unit":state.get("unit",""),"precursor_error_ppm":state.get("precursor_error_ppm")}
    if not form or key not in spectra:
        return dict(result,status="chemical_state_unresolved" if not form else "missing_spectrum")
    if state["fragmentation"] not in ("HCD","CID"):
        return dict(result,status="fragmentation_outside_by_model")
    data = spectra[key]
    if not data.get("charge") or int(data["charge"])!=z:
        return dict(result,status="charge_missing_or_conflicting")
    spec = spectrum(data)
    candidate = fragment_score(form,z,spec,ppm)
    matches = candidate.pop("observed")
    result.update(candidate)
    _,observed_intensity,_,_=spec
    cuts=set()
    for j,label in matches.items():
        if observed_intensity[j]>=observed_intensity.max()*.01:
            n=int(label.split("+")[0][1:])
            cuts.add(n if label.startswith("b") else len(peptide)-n)
    adjacent={n for n in (site-1,site) if 0<n<len(peptide)}
    result.update(adjacent_observed_cuts=sorted(cuts & adjacent),adjacent_required_cuts=sorted(adjacent),
                  both_available_adjacent_cuts_supported=bool(adjacent) and adjacent<=cuts)
    sa = angle(pred[(form,z)],spec,ppm) if state["fragmentation"]=="HCD" and (form,z) in pred else None
    result.update(prosit_angle=sa, status="scored", identical_chemical_competitors=[],
                  best_ms2=None,best_prosit=None,unscored_prosit_forms=0,unsupported_fragments=0)
    candidate_folded = form.replace("I","L")
    for other, roles in sorted(alternatives.items()):
        if other.replace("I","L")==candidate_folded:
            result["identical_chemical_competitors"].extend(sorted(roles))
            continue
        try:
            competitor = fragment_score(other,z,spec,ppm)
        except (ValueError,KeyError):
            result["unsupported_fragments"]+=1
            continue
        detail = {"form":other,"roles":sorted(roles),"score":competitor["explained_fraction"]}
        if result["best_ms2"] is None or detail["score"]>result["best_ms2"]["score"]:
            result["best_ms2"]=detail
        if state["fragmentation"]!="HCD" or (other,z) not in pred:
            result["unscored_prosit_forms"]+=1
            continue
        other_sa = angle(pred[(other,z)],spec,ppm)
        if other_sa is None:
            result["unscored_prosit_forms"]+=1
        elif result["best_prosit"] is None or other_sa>result["best_prosit"]["score"]:
            result["best_prosit"]={"form":other,"roles":sorted(roles),"score":other_sa}
    # Discriminating ions must have no theoretical counterpart within tolerance in the strongest MS2 alternative.
    if result["best_ms2"]:
        other_mz,_ = ions(result["best_ms2"]["form"],z)
        obs, it, _, _ = spec
        diagnostic=[]
        for j,label in matches.items():
            cut = int(label.split("+")[0][1:])
            cut = cut if label[0]=="b" else len(peptide)-cut
            contains_site = site<=cut if label[0]=="b" else site>cut
            if contains_site and it[j]>=it.max()*.01 and not np.any(abs(obs[j]-other_mz)<=other_mz*ppm*1e-6):
                diagnostic.append({"ion":label,"mz":float(obs[j]),"cut":cut})
        result["diagnostic_vs_best_ms2"]=diagnostic
        result["ms2_fraction_gap"]=candidate["explained_fraction"]-result["best_ms2"]["score"]
    if sa is not None and result["best_prosit"]:
        result["prosit_gap"]=sa-result["best_prosit"]["score"]
    return result


def build(cfg,out,requests_only=False):
    discovery = Path(cfg["discovery"])
    report=json.loads((discovery/"completion.json").read_text(encoding="utf-8"))
    if report["stage"]!="discover" or report["status"]!="complete":
        raise ValueError("A completed discovery run is required")
    for name in ("saap_candidates_relaxed.tsv","job.json"):
        if sha(discovery/name)!=report["artifacts"][name]:
            raise ValueError("Changed discovery artifact: "+name)
    job = json.loads((discovery/"job.json").read_text(encoding="utf-8"))
    if sha(Path(cfg["scans"]))!=job["inputs"]["scans"]["sha256"]:
        raise ValueError("The supplied scan table differs from the discovery input")
    names=report["result"]["raw_order"]
    chemistry=Chemistry(**report["result"].get("chemistry",{}))
    rows=list(read(discovery/"saap_candidates_relaxed.tsv"))
    if len(rows)!=report["result"]["relaxed_hypotheses"]:
        raise ValueError("Incomplete relaxed hypothesis table")
    groups=defaultdict(list)
    by_scan=defaultdict(list)
    for r in rows:
        key=(names[int(r["matching_raw_file_index"])],r["modification_scan_number"])
        ident=identity(r)
        groups[ident].append((key,r))
        by_scan[key].append(ident)
    scans={(r["Raw file"],r["Scan number"]):r for r in read(cfg["scans"]) if (r["Raw file"],r["Scan number"]) in by_scan}
    if set(scans)!=set(by_scan):
        raise ValueError("Candidate scan missing")
    wanted={(k,i) for k,ids in by_scan.items() for i in ids}
    _,frozen=prepared(cfg,wanted)
    base_keys={(Path(r["DP base raw file"].replace("\\","/")).stem,r["DP base scan number"]) for r in scans.values()}
    base_forms,closed=defaultdict(list),defaultdict(list)
    for field,keys,store in (("base_psms",base_keys,base_forms),("closed_psms",by_scan,closed)):
        if cfg.get(field):
            for r in read(cfg[field]):
                k=(r["Raw file"],r["Scan number"])
                if k in keys and r["Decoy"]!="+":
                    store[k].append(r)
    states=[]
    for ident,rs in groups.items():
        base,site,pep=ident
        for key,r in rs:
            sc=scans[key]
            bkey=(Path(sc["DP base raw file"].replace("\\","/")).stem,sc["DP base scan number"])
            base_psms=[p for p in base_forms[bkey] if p["Sequence"]==base]
            catalogs=frozen.get((key,ident),{})
            forms=set()
            target=None
            alternatives=defaultdict(set)
            for uid,catalog in catalogs.items():
                alts=defaultdict(set)
                for f,role in catalog:
                    if role!="candidate" and f:
                        alts[f].add(role)
                for form,role in catalog:
                    if role=="candidate":
                        states.append(dict(key=key,identity=ident,row=r,form=form,charge=int(sc["Charge"]),
                                           alternatives=alts,target=None,coverage="frozen_comparison_units",unit=uid))
            if catalogs and any(role=="candidate" for catalog in catalogs.values() for _,role in catalog):
                continue
            note="same_scan_mass_window; inherited_base_mods; up_to_two_base_PTM; alternate_sites; closed_PSM"
            if not forms:
                for bp in base_psms:
                    try:
                        mods=mq_form(base,bp["Modified sequence"])
                    except ValueError:
                        continue
                    bf=encode(base,mods,chemistry)
                    target=neutral(bf)+float(r["observed_delta_da"])
                    tol=max(.005,target*1e-5)
                    # A lost modification can move only to an allowed residue or be lost.
                    candidate_states=[mods]
                    if site in mods:
                        retained={pos:code for pos,code in mods.items() if pos!=site}
                        code=mods[site]
                        candidate_states=[retained]+[{**retained,i:code} for i,a in enumerate(pep,1)
                            if i not in retained and a in SITE_MODS.get(code,"")]
                    compatible={encode(pep,m,chemistry) for m in candidate_states}
                    compatible={f for f in compatible if abs(neutral(f)-target)<=tol}
                    if not compatible:
                        continue
                    forms.update(compatible)
                    for other_base,other_site,other in by_scan[key]:
                        if other_base==base and other_site not in mods:
                            f=encode(other,mods,chemistry)
                            if f not in compatible and abs(neutral(f)-target)<=tol:
                                alternatives[f].add("alternative_substitution")
                    for f in ptm_forms(base,mods,target,chemistry):
                        alternatives[f].add("base_PTM")
                    for psm in closed[key]:
                        if psm["Charge"]!=sc["Charge"]:
                            continue
                        try:
                            f=encode(psm["Sequence"],mq_form(psm["Sequence"],psm["Modified sequence"]),chemistry)
                        except ValueError:
                            continue
                        if abs(neutral(f)-target)<=tol:
                            ids=psm.get("Proteins","").split(";")
                            normal=any(i.startswith(("SFFREF","CON__")) for i in ids)
                            if f not in compatible or normal or psm.get("Potential contaminant")=="+":
                                alternatives[f].add("closed_reference_or_contaminant_PSM" if normal or
                                    psm.get("Potential contaminant")=="+" else "closed_search_PSM")
            for form in sorted(forms) or [""]:
                states.append(dict(key=key,identity=ident,row=r,form=form,charge=int(sc["Charge"]),
                                   alternatives={k:set(v) for k,v in alternatives.items()},
                                   target=neutral(form)+float(r["mass_error_da"]) if form else None,
                                   coverage=note))
    for s in states:
        s["fragmentation"]=scans[s["key"]]["Fragmentation"]
        if s["form"]:
            m = neutral(s["form"])
            measured = scans[s["key"]]["Mass"]
            if measured and math.isfinite(float(measured)):
                s["precursor_error_ppm"]=(float(measured)-m)/m*1e6
            if s["coverage"]!="frozen_comparison_units":
                for f in list(s["alternatives"]):
                    if abs(neutral(f)-s["target"])>max(.005,s["target"]*1e-5):
                        del s["alternatives"][f]
    targets={s["target"] for s in states if s["target"] is not None and s["coverage"]!="frozen_comparison_units"}
    hits=digest_index(cfg.get("digest_references",{}),targets,chemistry)
    for s in states:
        if s["coverage"]!="frozen_comparison_units" and cfg.get("digest_references"):
            s["coverage"]+="; reference_digest:"+",".join(chemistry.enzymes)
            for f,role in hits[s["target"]]:
                s["alternatives"].setdefault(f,set()).add(role)
    spectra={}
    if cfg.get("peaks"):
        data=json.loads(Path(cfg["peaks"]).read_text(encoding="utf-8-sig"))
        spectra={k:data[k[0]+":scan="+k[1]] for k in by_scan if k[0]+":scan="+k[1] in data}
    for raw,path in cfg.get("mzml",{}).items():
        missing={k[1] for k in by_scan if k[0]==raw and k not in spectra}
        if missing:
            spectra.update({(raw,n):v for n,v in mzml_peaks(Path(path),missing).items()})
    wanted_predictions={(f,s["charge"]) for s in states for f in [s["form"],*s["alternatives"]]
                        if f and s["fragmentation"]=="HCD"}
    model=cfg.get("model","Prosit_2024_intensity_PTMs_gl")
    energy=cfg.get("calibrated_nce")
    pred={}
    if cfg.get("predictions"):
        if not energy or not cfg.get("calibration_source"):
            raise ValueError("Provide calibrated NCE and source for cached predictions")
        pred=native_cache(Path(cfg["predictions"]),wanted_predictions,model,float(energy))
    requests=[{"peptide_sequences":f,"precursor_charges":z} for f,z in sorted(wanted_predictions-pred.keys())]
    tsv(out/"prediction_requests.tsv",["peptide_sequences","precursor_charges"],requests)
    if requests_only:
        return {"candidate_observations":len(rows),"missing_prediction_forms":len(requests),"stage":"prediction_inputs"}
    mapping=defaultdict(set)
    if cfg.get("context_mapping"):
        for r in read(cfg["context_mapping"]):
            mapping[(r["raw_file"],r["scan_number"],r["candidate_id"],r["candidate_peptide"])].add(
                "|".join(r[k] for k in ("source","source_accession","site_1based","route","coverage")))
    qs={}
    if cfg.get("class_q"):
        for r in read(cfg["class_q"]):
            if r["sequence"] in qs:
                raise ValueError("Duplicate sequence in scoped class-q table")
            qs[r["sequence"]]={k:r[k] for k in ("main_peptide_q","prosit_peptide_q","secondary_peptide_q","FDR_result")}
    extra=defaultdict(list)
    if cfg.get("prior_unit_evidence"):
        for r in read(cfg["prior_unit_evidence"]):
            ident=(r["base_sequence"],int(r["position_1based"]),r["candidate_sequence"])
            if ident in groups:
                extra[ident].append({k:r[k] for k in ("unit_id","scan_key","model","proforma_sequence","comparison_status",
                    "PIF","RT_residual","adjacent_cuts_supported","diagnostic_cuts_vs_base","spectral_angle",
                    "best_competitor_proforma","candidate_minus_competitor_SA","identical_chemical_normal_links")})
    abundance=defaultdict(list)
    if cfg.get("abundance"):
        events={r["event"]:r for r in read(cfg["paired_events"])}
        for a in read(cfg["abundance"]):
            e=events[a["event"]]
            for label in e["candidate_labels"].split(";"):
                ident=(e["base"],int(e["position"]),label)
                if ident in groups:
                    abundance[ident].append({k:a[k] for k in ("raw_file","modification_stratum","matched_charges",
                                              "candidate_to_base_signal_ratio","paired_xic_qc_pass","status")})
    evidence=defaultdict(list)
    for s in states:
        key,r=s["key"],s["row"]
        result=score_unit(s,spectra,pred,float(cfg.get("fragment_ppm",20)))
        result.update(raw=key[0],scan=key[1],id=r["candidate_id"],DP_PEP=r["posterior_error_probability"],
                      localization=r["localization_probability_at_site"],mass_error_Da=r["mass_error_da"],
                      DP_annotation=r["modification_annotation"],strict=r["passes_strict_discovery_filters"],
                      max_localization=r["max_probability"],DP_PTM_warning=r["ptm_or_artifact_warning"],
                      cluster_forward=r["cluster_size_forward"],cluster_reverse=r["cluster_size_reverse"],
                      peptide_length_difference=r["peptide_length_difference"],source_genes=r["gene_names"],
                      source_proteins=r["protein_ids"],contexts=sorted(mapping[(key[0],key[1],r["candidate_id"],s["identity"][2])]),
                      DP_score=r["dp_score"],observed_delta_Da=r["observed_delta_da"],cluster=r["cluster_index"],
                      PIF=(scans[key].get("Parent intensity fraction") if cfg.get("pif_calculated") else None),
                      closed_PSMs=[{k:v for k,v in p.items() if k in ("Sequence","Modified sequence","Proteins","PEP","Score","PIF","Potential contaminant","Charge","id") or k.startswith("Reporter intensity")} for p in closed[key]])
        evidence[s["identity"]].append(result)
    fields=["base_peptide","position_1based","substitution","candidate_peptide","IL_event","DP_hypotheses",
            "physical_scans","strict_hypotheses","min_DP_PEP","max_site_probability","MS2_scored_forms",
            "Prosit_scored_forms","best_MS2_competitor","MS2_fraction_gap","best_Prosit_competitor",
            "Prosit_angle_gap","representative_raw_scan","evidence_status","prosit_model","calibrated_NCE",
            "calibration_source","closed_library_matches_current_build","class_q_json","all_scan_evidence_json","prior_unit_evidence_json","paired_abundance_json",
            "representative_MS2_explained_fraction","competitor_MS2_explained_fraction","representative_Prosit_angle",
            "competitor_Prosit_angle","representative_diagnostic_ions","证据注释"]
    final=[]
    for ident,rs in sorted(groups.items()):
        base,site,pep=ident
        ev=evidence[ident]
        def finite(values,func):
            vals=[float(x) for x in values if x!="" and math.isfinite(float(x))]
            return func(vals) if vals else ""
        # Fixed display rule; keep all scans/forms, including contradictory evidence, in the same row.
        display=sorted(ev,key=lambda x:(x["status"]!="scored",x.get("prosit_gap") is None,
                                          x["raw"],int(x["scan"]),x.get("candidate_form","")))[0]
        status=sorted({e["status"] for e in ev})
        if any(e.get("identical_chemical_competitors") for e in ev):status.append("same_chemical_identity")
        if any(e.get("prosit_gap",0)<0 or e.get("ms2_fraction_gap",0)<0 for e in ev):status.append("stronger_competitor")
        if any(e.get("prosit_angle") is None or e.get("unscored_prosit_forms") for e in ev):status.append("prediction_coverage_incomplete")
        if not any(r["passes_strict_discovery_filters"]=="True" for k,r in rs):status.append("relaxed_only")
        if any(e.get("precursor_error_ppm") is None or abs(e["precursor_error_ppm"])>10 for e in ev):status.append("precursor_review")
        if not any(e.get("best_ms2") for e in ev):status.append("no_distinct_scored_competitor")
        strict_count=sum(r["passes_strict_discovery_filters"]=="True" for _,r in rs)
        scan_count=len({k for k,_ in rs})
        observed_count=sum(e["status"]=="scored" for e in ev)
        predicted_count=sum(e.get("prosit_gap") is not None for e in ev)
        note=[f"{scan_count} 张发现谱，{len(rs)} 条质量相容假设，{strict_count} 条通过严格初筛。",
              f"{observed_count} 个扫描与化学形式组合有实测碎片评分，{predicted_count} 个可比较双方预测谱。"]
        if "same_chemical_identity" in status:note.append("有同化学身份或 I/L 等价解释，谱图无法区分其来源。")
        if "stronger_competitor" in status:note.append("至少一个比较单元有得分更高的竞争解释。")
        if "precursor_review" in status:note.append("部分前体质量误差缺失或超过 10 ppm，须核对。")
        if "chemical_state_unresolved" in status:note.append("部分基底修饰或候选化学形式未解析。")
        if "missing_spectrum" in status:note.append("部分扫描缺少实测峰。")
        if "prediction_coverage_incomplete" in status:note.append("模型或预测缓存未覆盖部分比较。")
        event=base.replace("L","I")+"|"+str(site)+"|"+pep.replace("L","I")
        final.append(dict(zip(fields,[base,site,base[site-1]+">"+pep[site-1],pep,"EV_"+hashlib.sha256(event.encode()).hexdigest()[:16],
            len(rs),len({k for k,r in rs}),sum(r["passes_strict_discovery_filters"]=="True" for k,r in rs),
            finite((r["posterior_error_probability"] for k,r in rs),min),finite((r["localization_probability_at_site"] for k,r in rs),max),
            sum(e["status"]=="scored" for e in ev),sum(e.get("prosit_angle") is not None for e in ev),
            (display.get("best_ms2") or {}).get("form",""),display.get("ms2_fraction_gap",""),
            (display.get("best_prosit") or {}).get("form",""),display.get("prosit_gap",""),display["raw"]+":scan="+display["scan"],
            ";".join(status),model,energy or "",cfg.get("calibration_source",""),
            cfg.get("closed_export_fasta_matches_generated","not_checked"),packed(qs.get(pep,{"status":"not_estimated"})),
            packed(ev),packed(extra[ident]),packed(abundance[ident]),display.get("explained_fraction",""),
            (display.get("best_ms2") or {}).get("score",""),display.get("prosit_angle",""),
            (display.get("best_prosit") or {}).get("score",""),packed(display.get("diagnostic_vs_best_ms2",[])),
            "".join(note)])))
    if sum(r["DP_hypotheses"] for r in final)!=len(rows):
        raise ValueError("Lost DP observations during aggregation")
    tsv(out/"SAAP_evidence.tsv",fields,final)
    return {"candidate_observations":len(rows),"substitution_rows":len(final),"states":len(states),
            "observed_spectra":len(spectra),"matched_cached_predictions":len(pred),"missing_prediction_forms":len(requests),
            "model":model,"nce":energy,"table_sha256":sha(out/"SAAP_evidence.tsv")}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("config",type=Path)
    ap.add_argument("--out",required=True,type=Path)
    ap.add_argument("--requests-only",action="store_true")
    args=ap.parse_args()
    cfg=json.loads(args.config.read_text(encoding="utf-8-sig"))
    args.out.mkdir(parents=True,exist_ok=False)
    (args.out/"inputs.json").write_text(packed(cfg),encoding="utf-8")
    result=build(cfg,args.out,args.requests_only)
    (args.out/"completion.json").write_text(packed(result),encoding="utf-8")
    print(packed(result))


if __name__=="__main__":
    main()

