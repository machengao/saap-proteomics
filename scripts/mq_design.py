"""Render two version-matched MaxQuant jobs from the experiment design resolved by the skill."""
import argparse
import copy
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path, PureWindowsPath

INSTRUMENTS = ["Orbitrap","Bruker Q-TOF","Sciex Q-TOF","Agilent Q-TOF",
               "Bruker TIMS","Astral","Waters Q-TOF","Shimadzu Q-TOF"]
ENZYME_MODES = ["Specific","SemispecificFreeNterm","SemispecificFreeCterm","Semispecific","Unspecific","None"]
DEFAULT_TEMPLATE = Path(__file__).resolve().parent.parent / "assets" / "orbitrap-dda-2.8.1.xml"


def set_value(node, key, value):
    leaf = node.find(key)
    if leaf is None:
        raise ValueError("Missing template field: " + key)
    leaf.text = str(value)


def array(node, key, tag, values):
    element = node.find(key)
    if element is None:
        raise ValueError("Missing template array: " + key)
    element.clear()
    for value in values:
        ET.SubElement(element, tag).text = str(value)


def render(spec, phase):
    design = spec["groups"]
    if not design or not spec["files"] or not isinstance(spec["include_contaminants"], bool):
        raise ValueError("Expected groups, acquisition files and an explicit contaminant decision.")
    template = ET.parse(design[0].get("template",DEFAULT_TEMPLATE))
    root = template.getroot()
    if root.findtext("maxQuantVersion") != "2.8.1.0":
        raise ValueError("Use the audited 2.8.1.0 template or adapt the serializer to the installed version.")
    conf = Path(spec["maxquant_conf"])
    known = {name: {n.get("title") for n in ET.parse(conf / (name + ".xml")).getroot()}
             for name in ("enzymes", "modifications")}
    local = conf / "modifications.local.xml"
    if local.is_file():
        known["modifications"].update(n.get("title") for n in ET.parse(local).getroot())
    groups = root.find("parameterGroups")
    groups.clear()
    for item in design:
        src = ET.parse(item.get("template",DEFAULT_TEMPLATE)).getroot()
        if src.findtext("maxQuantVersion") != root.findtext("maxQuantVersion"):
            raise ValueError("Group template versions differ")
        group = copy.deepcopy(src.findall("./parameterGroups/parameterGroup")[item.get("template_group", 0)])
        if item["acquisition"] != "DDA":
            raise ValueError("Direct DP requires a supported DDA acquisition. Route DIA separately.")
        if not item.get("instrument") or not item.get("fragmentation") or not item.get("basis"):
            raise ValueError("Document instrument, fragmentation and parameter evidence per group.")
        family = item["instrument_family"]
        if INSTRUMENTS.index(family) != int(group.findtext("msInstrument")):
            raise ValueError("Select a native template for the identified instrument family: "+family)
        set_value(group,"enzymeMode",ENZYME_MODES.index(item.get("enzyme_mode","Specific")))
        # Global detector/fragment profiles apply to all groups, so require compatible templates.
        for tag in ("msmsParamsArray", "fragmentationParamsArray"):
            if ET.tostring(root.find(tag)) != ET.tostring(src.find(tag)):
                raise ValueError("Group templates have different global detector or fragment settings: " + tag)
        for field, catalog, values in (("enzymes","enzymes",item["enzymes"]),
                    ("fixedModifications","modifications",item["fixed_modifications"]),
                    ("variableModifications","modifications",item["variable_modifications"])):
            missing = set(values) - known[catalog]
            if missing or len(values) != len(set(values)):
                raise ValueError(f"Unknown or duplicate {field}: {missing or values}")
            array(group, field, "string", values)
        if not item["enzymes"] or set(item["fixed_modifications"]) & set(item["variable_modifications"]):
            raise ValueError("Check enzyme and overlapping fixed/variable modifications.")
        missed = item["missed_cleavages"]
        if type(missed) is not int or missed < 0:
            raise ValueError("Invalid missed cleavage count")
        set_value(group, "maxMissedCleavages", missed)
        mode = item["quantification"]
        if mode in ("unlabelled", "MaxLFQ", "TMT6-MS2", "TMT6-MS3"):
            if group.findtext("lcmsRunType") not in ("Standard", "Reporter MS2", "Reporter MS3"):
                raise ValueError("Select a matching native template for this instrument/acquisition.")
            labelled = mode.startswith("TMT")
            set_value(group, "lcmsRunType", "Reporter " + mode[-3:] if labelled else "Standard")
            set_value(group, "lfqMode", 1 if mode == "MaxLFQ" else 0)
            set_value(group, "multiplicity", 1)
            if mode=="MaxLFQ":
                set_value(group,"lfqMinRatioCount",item.get("lfq_min_ratio_count",2))
            array(group, "labelMods", "string", [""])
            labels = group.find("isobaricLabels")
            labels.clear()
            if labelled:
                if item.get("reporter_tolerance") is not None:
                    if not math.isfinite(item["reporter_tolerance"]) or item["reporter_tolerance"]<=0:
                        raise ValueError("Invalid reporter ion tolerance")
                    set_value(group,"reporterMassTolerance",item["reporter_tolerance"])
                if any("TMT" in m.upper() for m in item["fixed_modifications"] + item["variable_modifications"]):
                    raise ValueError("TMT is already encoded by reporter labels.")
                corrections = item["channel_corrections"]
                for channel in range(126, 132):
                    coeff = corrections[str(channel)]
                    if len(coeff) != 4 or not all(math.isfinite(x) and 0 <= x <= 100 for x in coeff):
                        raise ValueError(f"Invalid isotope corrections for {channel}")
                    row = ET.SubElement(labels, "IsobaricLabelInfo")
                    for tag, prefix in (("internalLabel", "Lys"), ("terminalLabel", "Nter")):
                        mod = f"TMT6plex-{prefix}{channel}"
                        if mod not in known["modifications"]:
                            raise ValueError("Label absent from local library: " + mod)
                        ET.SubElement(row, tag).text = mod
                    for tag, val in zip(("correctionFactorM2","correctionFactorM1","correctionFactorP1","correctionFactorP2"), coeff):
                        ET.SubElement(row, tag).text = str(val)
                    ET.SubElement(row, "tmtLike").text = "True"
                if mode == "TMT6-MS3":
                    set_value(root, "writeMs3ScansTable", True)
        elif mode != "native-template":
            raise ValueError("Unsupported quantification setting: " + mode)
        groups.append(group)
    files = spec["files"]
    paths = [f[phase + "_path"] for f in files]
    names = [PureWindowsPath(p).stem for p in paths]
    if (len(names) != len(set(names)) or any(not f["experiment"] for f in files)
            or any(type(f.get("fraction",32767)) != int or not 0 < f.get("fraction",32767) <= 32767 for f in files)
            or any(not 0 <= f["group"] < len(design) for f in files)):
        raise ValueError("Duplicate RAW or invalid experiment/group assignment.")
    for f in files:
        a, b = (PureWindowsPath(f[k]) for k in ("dp_path", "closed_path"))
        if a.name != b.name or a.parent == b.parent:
            raise ValueError("Use independent search directories and the same original RAW filename.")
    for key, tag, values in (("filePaths","string",paths),("experiments","string",[f["experiment"] for f in files]),
            ("fractions","short",[f.get("fraction",32767) for f in files]), ("ptms","boolean",["False"]*len(files)),
            ("paramGroupIndices","int",[f["group"] for f in files]),("referenceChannel","string",[f.get("reference_channel","") for f in files])):
        array(root, key, tag, values)
    fasta = root.findall("fastaFiles/FastaFileInfo")
    if len(fasta) != 1 or any(len(root.find(k)) for k in ("fastaFilesProteogenomics","fastaFilesFirstSearch")):
        raise ValueError("Use a template with a single reviewed FASTA parser.")
    set_value(fasta[0], "fastaFilePath", spec[phase+"_fasta"])
    if phase == "closed":
        set_value(fasta[0], "identifierParseRule", r">[^|]*\|(.*?)\|")
    for key, val in {"dependentPeptides":phase=="dp","massDifferenceSearch":False,"calcPeakProperties":True,
                     "includeContaminants":spec["include_contaminants"],"writeMsmsScansTable":True,
                     "writeMsScansTable":True,"matchBetweenRuns":False,"numThreads":spec["threads"],
                     "fixedSearchFolder":"","customTxtFolder":""}.items():
        set_value(root,key,val)
    return template


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("design", type=Path)
    ap.add_argument("--out", required=True, type=Path)
    args = ap.parse_args()
    spec = json.loads(args.design.read_text(encoding="utf-8-sig"))
    if not 1 <= spec["threads"] or spec.get("unresolved") or not spec["evidence_sources"]:
        raise ValueError("Resolve required design fields and document their sources before rendering.")
    trees = {phase: render(spec, phase) for phase in ("dp","closed")}
    args.out.mkdir(parents=True, exist_ok=False)
    for phase, tree in trees.items():
        tree.write(args.out / f"mqpar_{phase}.xml", encoding="utf-8", xml_declaration=True, short_empty_elements=False)
    (args.out / "parameter_decisions.json").write_text(json.dumps(spec, ensure_ascii=False, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()

