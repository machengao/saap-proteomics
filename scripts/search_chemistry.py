"""Recover a homogeneous cleavage and residue-mass profile from the actual MaxQuant search."""
import json
import xml.etree.ElementTree as ET
from pathlib import Path
from functools import cached_property
from dataclasses import dataclass

from discovery import AA_MASS

TMT6=229.162932
PAIRS=Path(__file__).resolve().parent.parent/"assets/enzyme-pairs-2.8.1.json"


@dataclass(frozen=True)
class Chemistry:
    cam: bool=True
    tmt: bool=False
    enzymes: tuple=("Trypsin/P",)
    missed: int=2
    min_length: int=7
    max_mass: float=5000.

    @cached_property
    def masses(self):
        values=dict(AA_MASS)
        if not self.cam:values["C"]=103.009185
        if self.tmt:values["K"]+=TMT6
        return values

    @cached_property
    def pairs(self):
        catalog=json.loads(PAIRS.read_text(encoding="utf-8"))
        if not self.enzymes or any(e not in catalog for e in self.enzymes):
            raise ValueError(f"No verified cleavage table for {self.enzymes}")
        return set().union(*(catalog[e] for e in self.enzymes))

    def digestible(self,parent,start,peptide):
        end=start+len(peptide)
        return ((start==0 or parent[start-1:start+1] in self.pairs)
                and (end==len(parent) or parent[end-1:end+1] in self.pairs)
                and sum(peptide[i:i+2] in self.pairs for i in range(len(peptide)-1))<=self.missed)

    def mass(self,seq):
        return 18.010565+(TMT6 if self.tmt else 0.)+sum(self.masses[a] for a in seq)

    @classmethod
    def from_xml(cls,path):
        root=ET.parse(path).getroot()
        found=set()
        for group in root.findall("parameterGroups/parameterGroup"):
            fixed={e.text for e in group.findall("fixedModifications/string")}
            mode=group.findtext("lcmsRunType")
            labels=group.findall("isobaricLabels/IsobaricLabelInfo")
            if (mode not in ("Standard","Reporter MS2","Reporter MS3")
                    or group.findtext("multiplicity")!="1" or any((e.text or "").strip() for e in group.findall("labelMods/string"))
                    or fixed-{"Carbamidomethyl (C)"} or group.findtext("enzymeMode")!="0"
                    or root.findtext("independentEnzymes","false").lower()=="true"):
                raise ValueError("Unsupported fixed chemistry, labelling or digestion profile; review actual mqpar.")
            tmt=mode.startswith("Reporter")
            if tmt:
                expected={(f"TMT6plex-Lys{n}",f"TMT6plex-Nter{n}") for n in range(126,132)}
                if len(labels)!=6 or {(l.findtext("internalLabel"),l.findtext("terminalLabel")) for l in labels}!=expected:
                    raise ValueError("Reporter chemistry currently supports a complete TMT6 channel set.")
            elif labels:
                raise ValueError("Reporter labels conflict with Standard acquisition")
            found.add(cls("Carbamidomethyl (C)" in fixed,tmt,tuple(e.text for e in group.findall("enzymes/string")),
                          int(group.findtext("maxMissedCleavages")),int(root.findtext("minPeptideLength")),
                          float(root.findtext("maxPeptideMass"))))
        if len(found)!=1:
            raise ValueError("Run distinct chemical/digestion profiles as separate cohorts before combining results.")
        profile=found.pop()
        if not profile.pairs or profile.missed<0 or profile.min_length<1 or profile.max_mass<=0:
            raise ValueError("Invalid chemistry or digestion parameters")
        return profile

