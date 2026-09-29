#!/usr/bin/env python3
"""Reproducible 30-compound developability demonstration, Python >=3.10.

Published CNS-MPO transforms and Delaney ESOL coefficients are implemented
explicitly. RDKit descriptors replace the original proprietary descriptors.
All pKa, logD, hERG, Caco-2, HIA, and conformational estimates are exploratory;
none is a validated assay prediction. No network access is used at runtime.
Run ``python run_task1_mpo_admet_developability.py --help`` for options.
SPDX-License-Identifier: MIT
"""
from __future__ import annotations

import argparse
from collections import Counter
from datetime import datetime, timezone
import hashlib
import json
import logging
import math
import multiprocessing as mp
from pathlib import Path
import platform
import sys
import tempfile
from typing import Any

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import rdkit
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, ChemicalFeatures, Crippen, Descriptors, Lipinski
from rdkit.Chem import rdMolDescriptors
import seaborn as sns

VERSION = "1.0.0"
LOG = logging.getLogger("task1")
GROUPS = ("Oral Drugs", "Toxic Dropouts", "bRo5 Modalities")
PALETTE = ("#187F88", "#D2664B", "#735CB0")
WARNING = ("Illustrative benchmark, not external validation. pKa/logD and ADMET "
           "surrogates are uncalibrated; hERG is not overall safety. 3D geometry "
           "does not establish solvent exposure or absolute chameleonicity.")

# Frozen PubChem parent structures, with stereochemistry exactly as deposited.
# Embedded to make this script independently runnable without companion data.
PANEL_JSON = r'''[
  {
    "name": "Aspirin",
    "archetype": "Oral Drugs",
    "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O",
    "pubchem_cid": 2244,
    "molecular_formula": "C9H8O4",
    "inchikey": "BSYNRYMUTXBXSQ-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Aspirin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/2244",
    "retrieved_utc": "2026-09-10T05:32:50.211318+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/2244"
  },
  {
    "name": "Diazepam",
    "archetype": "Oral Drugs",
    "smiles": "CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3",
    "pubchem_cid": 3016,
    "molecular_formula": "C16H13ClN2O",
    "inchikey": "AAOVKJBEBIDNHE-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Diazepam/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/3016",
    "retrieved_utc": "2026-09-10T05:32:51.850789+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/3016"
  },
  {
    "name": "Sildenafil",
    "archetype": "Oral Drugs",
    "smiles": "CCCC1=NN(C2=C1N=C(NC2=O)C3=C(C=CC(=C3)S(=O)(=O)N4CCN(CC4)C)OCC)C",
    "pubchem_cid": 135398744,
    "molecular_formula": "C22H30N6O4S",
    "inchikey": "BNRNXUUZRGQAQC-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Sildenafil/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/135398744",
    "retrieved_utc": "2026-09-10T05:32:53.792600+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/135398744"
  },
  {
    "name": "Gefitinib",
    "archetype": "Oral Drugs",
    "smiles": "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4",
    "pubchem_cid": 123631,
    "molecular_formula": "C22H24ClFN4O3",
    "inchikey": "XGALLCVXEZPNRQ-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Gefitinib/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/123631",
    "retrieved_utc": "2026-09-10T05:32:55.754835+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/123631"
  },
  {
    "name": "Omeprazole",
    "archetype": "Oral Drugs",
    "smiles": "CC1=CN=C(C(=C1OC)C)CS(=O)C2=NC3=C(N2)C=C(C=C3)OC",
    "pubchem_cid": 4594,
    "molecular_formula": "C17H19N3O3S",
    "inchikey": "SUBDBMMJDZJVOS-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Omeprazole/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/4594",
    "retrieved_utc": "2026-09-10T05:32:57.737976+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/4594"
  },
  {
    "name": "Metoprolol",
    "archetype": "Oral Drugs",
    "smiles": "CC(C)NCC(COC1=CC=C(C=C1)CCOC)O",
    "pubchem_cid": 4171,
    "molecular_formula": "C15H25NO3",
    "inchikey": "IUBSYMUCCVWXPE-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Metoprolol/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/4171",
    "retrieved_utc": "2026-09-10T05:32:59.705667+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/4171"
  },
  {
    "name": "Losartan",
    "archetype": "Oral Drugs",
    "smiles": "CCCCC1=NC(=C(N1CC2=CC=C(C=C2)C3=CC=CC=C3C4=NNN=N4)CO)Cl",
    "pubchem_cid": 3961,
    "molecular_formula": "C22H23ClN6O",
    "inchikey": "PSIFNNKUMBGKDQ-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Losartan/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/3961",
    "retrieved_utc": "2026-09-10T05:33:01.575261+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/3961"
  },
  {
    "name": "Warfarin",
    "archetype": "Oral Drugs",
    "smiles": "CC(=O)CC(C1=CC=CC=C1)C2=C(C3=CC=CC=C3OC2=O)O",
    "pubchem_cid": 54678486,
    "molecular_formula": "C19H16O4",
    "inchikey": "PJVWKTKQMONHTI-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Warfarin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/54678486",
    "retrieved_utc": "2026-09-10T05:33:03.754480+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/54678486"
  },
  {
    "name": "Fluconazole",
    "archetype": "Oral Drugs",
    "smiles": "C1=CC(=C(C=C1F)F)C(CN2C=NC=N2)(CN3C=NC=N3)O",
    "pubchem_cid": 3365,
    "molecular_formula": "C13H12F2N6O",
    "inchikey": "RFHAOTPXVQNOHP-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Fluconazole/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/3365",
    "retrieved_utc": "2026-09-10T05:33:05.660445+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/3365"
  },
  {
    "name": "Captopril",
    "archetype": "Oral Drugs",
    "smiles": "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O",
    "pubchem_cid": 44093,
    "molecular_formula": "C9H15NO3S",
    "inchikey": "FAKRSMQSSFJEIM-RQJHMYQMSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Captopril/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/44093",
    "retrieved_utc": "2026-09-10T05:33:07.690954+00:00",
    "clinical_note": "Established oral clinical reference; parent active moiety, not a commercial salt.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/44093"
  },
  {
    "name": "Terfenadine",
    "archetype": "Toxic Dropouts",
    "smiles": "CC(C)(C)C1=CC=C(C=C1)C(CCCN2CCC(CC2)C(C3=CC=CC=C3)(C4=CC=CC=C4)O)O",
    "pubchem_cid": 5405,
    "molecular_formula": "C32H41NO2",
    "inchikey": "GUGOEEXESWIERI-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Terfenadine/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5405",
    "retrieved_utc": "2026-09-10T05:33:09.331402+00:00",
    "clinical_note": "Historical withdrawal; hERG/QT liability and exposure-dependent interaction risk.",
    "clinical_source": "https://pubmed.ncbi.nlm.nih.gov/8772706/"
  },
  {
    "name": "Astemizole",
    "archetype": "Toxic Dropouts",
    "smiles": "COC1=CC=C(C=C1)CCN2CCC(CC2)NC3=NC4=CC=CC=C4N3CC5=CC=C(C=C5)F",
    "pubchem_cid": 2247,
    "molecular_formula": "C28H31FN4O",
    "inchikey": "GXDALQBWZGODGZ-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Astemizole/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/2247",
    "retrieved_utc": "2026-09-10T05:33:11.386416+00:00",
    "clinical_note": "Historical withdrawal; QT prolongation and ventricular arrhythmia risk.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/2247"
  },
  {
    "name": "Cisapride",
    "archetype": "Toxic Dropouts",
    "smiles": "CO[C@H]1CN(CC[C@H]1NC(=O)C2=CC(=C(C=C2OC)N)Cl)CCCOC3=CC=C(C=C3)F",
    "pubchem_cid": 6917698,
    "molecular_formula": "C23H29ClFN3O4",
    "inchikey": "DCSUBABJRXZOMT-IRLDBZIGSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Cisapride/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/6917698",
    "retrieved_utc": "2026-09-10T05:33:14.033900+00:00",
    "clinical_note": "Withdrawn or restricted in multiple markets; QT/arrhythmia liability.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/6917698"
  },
  {
    "name": "Cerivastatin",
    "archetype": "Toxic Dropouts",
    "smiles": "CC(C)C1=C(C(=C(C(=N1)C(C)C)COC)C2=CC=C(C=C2)F)/C=C/[C@H](C[C@H](CC(=O)O)O)O",
    "pubchem_cid": 446156,
    "molecular_formula": "C26H34FNO5",
    "inchikey": "SEERZIQQUAZTOL-ANMDKAQQSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Cerivastatin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/446156",
    "retrieved_utc": "2026-09-10T05:33:16.374539+00:00",
    "clinical_note": "Historical withdrawal for rhabdomyolysis; not a hERG-positive label.",
    "clinical_source": "https://www.accessdata.fda.gov/drugsatfda_docs/nda/2003/21-366_Crestor_Admindocs_P1.pdf"
  },
  {
    "name": "Troglitazone",
    "archetype": "Toxic Dropouts",
    "smiles": "CC1=C(C2=C(CCC(O2)(C)COC3=CC=C(C=C3)CC4C(=O)NC(=O)S4)C(=C1O)C)C",
    "pubchem_cid": 5591,
    "molecular_formula": "C24H27NO5S",
    "inchikey": "GXPHKUHSUJUWKP-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Troglitazone/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5591",
    "retrieved_utc": "2026-09-10T05:33:18.758954+00:00",
    "clinical_note": "Historical withdrawal for hepatotoxicity; not a hERG-positive label.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/5591"
  },
  {
    "name": "Rofecoxib",
    "archetype": "Toxic Dropouts",
    "smiles": "CS(=O)(=O)C1=CC=C(C=C1)C2=C(C(=O)OC2)C3=CC=CC=C3",
    "pubchem_cid": 5090,
    "molecular_formula": "C17H14O4S",
    "inchikey": "RZJQGNCSTQAWON-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Rofecoxib/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5090",
    "retrieved_utc": "2026-09-10T05:33:21.195382+00:00",
    "clinical_note": "Historical withdrawal for cardiovascular thrombotic risk; not a hERG-positive label.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/5090"
  },
  {
    "name": "Nefazodone",
    "archetype": "Toxic Dropouts",
    "smiles": "CCC1=NN(C(=O)N1CCOC2=CC=CC=C2)CCCN3CCN(CC3)C4=CC(=CC=C4)Cl",
    "pubchem_cid": 4449,
    "molecular_formula": "C25H32ClN5O2",
    "inchikey": "VRBKIVRKKCLPHA-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Nefazodone/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/4449",
    "retrieved_utc": "2026-09-10T05:33:23.476696+00:00",
    "clinical_note": "Hepatotoxicity; withdrawn in some markets, not a universal withdrawal.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/4449"
  },
  {
    "name": "Benoxaprofen",
    "archetype": "Toxic Dropouts",
    "smiles": "CC(C1=CC2=C(C=C1)OC(=N2)C3=CC=C(C=C3)Cl)C(=O)O",
    "pubchem_cid": 39941,
    "molecular_formula": "C16H12ClNO3",
    "inchikey": "MITFXPHMIHQXPI-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Benoxaprofen/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/39941",
    "retrieved_utc": "2026-09-10T05:33:25.960031+00:00",
    "clinical_note": "Historical withdrawal associated with hepatic toxicity.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/39941"
  },
  {
    "name": "Ximelagatran",
    "archetype": "Toxic Dropouts",
    "smiles": "CCOC(=O)CN[C@H](C1CCCCC1)C(=O)N2CC[C@H]2C(=O)NCC3=CC=C(C=C3)/C(=N/O)/N",
    "pubchem_cid": 9574101,
    "molecular_formula": "C24H35N5O5",
    "inchikey": "ZXIBCJHYVWYIKI-PZJWPPBQSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Ximelagatran/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/9574101",
    "retrieved_utc": "2026-09-10T05:33:29.831130+00:00",
    "clinical_note": "Withdrawn after liver injury reports; not a hERG-positive label.",
    "clinical_source": "https://www.ema.europa.eu/en/medicines/human/EPAR/ximelagatran-36-mg-film-coated-tablets"
  },
  {
    "name": "Fialuridine",
    "archetype": "Toxic Dropouts",
    "smiles": "C1=C(C(=O)NC(=O)N1[C@H]2[C@H]([C@@H]([C@H](O2)CO)O)F)I",
    "pubchem_cid": 50313,
    "molecular_formula": "C9H10FIN2O5",
    "inchikey": "IPVFGAYTKQKGBM-BYPJNBLXSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Fialuridine/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/50313",
    "retrieved_utc": "2026-09-10T05:33:33.343630+00:00",
    "clinical_note": "Investigational clinical failure with severe mitochondrial/hepatic toxicity.",
    "clinical_source": "https://pubmed.ncbi.nlm.nih.gov/7565947/"
  },
  {
    "name": "Cyclosporine A",
    "archetype": "bRo5 Modalities",
    "smiles": "CC[C@H]1C(=O)N(CC(=O)N([C@H](C(=O)N[C@H](C(=O)N([C@H](C(=O)N[C@H](C(=O)N[C@@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N([C@H](C(=O)N1)[C@@H]([C@H](C)C/C=C/C)O)C)C(C)C)C)CC(C)C)C)CC(C)C)C)C)C)CC(C)C)C)C(C)C)CC(C)C)C)C",
    "pubchem_cid": 5284373,
    "molecular_formula": "C62H111N11O12",
    "inchikey": "PMATZTZNYRCHOR-CGLBZJNRSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/cyclosporine/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5284373",
    "retrieved_utc": "2026-09-10T05:33:35.225066+00:00",
    "clinical_note": "Oral cyclic peptide; formulation-dependent exposure.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/5284373"
  },
  {
    "name": "Venetoclax",
    "archetype": "bRo5 Modalities",
    "smiles": "CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C",
    "pubchem_cid": 49846579,
    "molecular_formula": "C45H50ClN7O7S",
    "inchikey": "LQBVNQSMGBZMKD-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Venetoclax/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/49846579",
    "retrieved_utc": "2026-09-10T05:33:37.635967+00:00",
    "clinical_note": "Oral nonmacrocyclic BCL-2 inhibitor; amorphous dispersion and food effects.",
    "clinical_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9338003/"
  },
  {
    "name": "ARV-110",
    "archetype": "bRo5 Modalities",
    "smiles": "C1CC(CCC1NC(=O)C2=NN=C(C=C2)N3CCC(CC3)CN4CCN(CC4)C5=C(C=C6C(=C5)C(=O)N(C6=O)C7CCC(=O)NC7=O)F)OC8=CC(=C(C=C8)C#N)Cl",
    "pubchem_cid": 134414307,
    "molecular_formula": "C41H43ClFN9O6",
    "inchikey": "CLCTZVRHDOAUGJ-UHFFFAOYSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/bavdegalutamide/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/134414307",
    "retrieved_utc": "2026-09-10T05:33:39.632677+00:00",
    "clinical_note": "Bavdegalutamide; investigational oral PROTAC reference, not an approved-drug claim. PubChem record leaves stereochemistry unspecified.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/summary/summary.cgi?sid=441604884"
  },
  {
    "name": "Rifampicin",
    "archetype": "bRo5 Modalities",
    "smiles": "C[C@H]1/C=C/C=C(\\C(=O)NC2=C(C(=C3C(=C2O)C(=C(C4=C3C(=O)[C@](O4)(O/C=C/[C@@H]([C@H]([C@H]([C@@H]([C@@H]([C@@H]([C@H]1O)C)O)C)OC(=O)C)C)OC)C)C)O)O)/C=N/N5CCN(CC5)C)/C",
    "pubchem_cid": 135398735,
    "molecular_formula": "C43H58N4O12",
    "inchikey": "JQXXHWHPUNPDRT-WLSIYKJHSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Rifampicin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/135398735",
    "retrieved_utc": "2026-09-10T05:33:41.430105+00:00",
    "clinical_note": "Oral ansamycin macrocycle.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/135398735"
  },
  {
    "name": "Paclitaxel",
    "archetype": "bRo5 Modalities",
    "smiles": "CC1=C2[C@H](C(=O)[C@@]3([C@H](C[C@@H]4[C@]([C@H]3[C@@H]([C@@](C2(C)C)(C[C@@H]1OC(=O)[C@@H]([C@H](C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C",
    "pubchem_cid": 36314,
    "molecular_formula": "C47H51NO14",
    "inchikey": "RCINICONZNJXQF-MZXODVADSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Paclitaxel/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/36314",
    "retrieved_utc": "2026-09-10T05:33:43.525737+00:00",
    "clinical_note": "Taxane; conventional systemic products are intravenous. Included as a challenging bRo5 counterexample, not an oral success.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/36314"
  },
  {
    "name": "Tacrolimus",
    "archetype": "bRo5 Modalities",
    "smiles": "C[C@@H]1C[C@@H]([C@@H]2[C@H](C[C@H]([C@@](O2)(C(=O)C(=O)N3CCCC[C@H]3C(=O)O[C@@H]([C@@H]([C@H](CC(=O)[C@@H](/C=C(/C1)\\C)CC=C)O)C)/C(=C/[C@@H]4CC[C@H]([C@@H](C4)OC)O)/C)O)C)OC)OC",
    "pubchem_cid": 445643,
    "molecular_formula": "C44H69NO12",
    "inchikey": "QJJXYPPXXYFBGM-LFZNUXCKSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Tacrolimus/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/445643",
    "retrieved_utc": "2026-09-10T05:33:45.471094+00:00",
    "clinical_note": "Oral macrocyclic immunosuppressant.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/445643"
  },
  {
    "name": "Sirolimus",
    "archetype": "bRo5 Modalities",
    "smiles": "C[C@@H]1CC[C@H]2C[C@@H](/C(=C/C=C/C=C/[C@H](C[C@H](C(=O)[C@@H]([C@@H](/C(=C/[C@H](C(=O)C[C@H](OC(=O)[C@@H]3CCCCN3C(=O)C(=O)[C@@]1(O2)O)[C@H](C)C[C@@H]4CC[C@H]([C@@H](C4)OC)O)C)/C)O)OC)C)C)/C)OC",
    "pubchem_cid": 5284616,
    "molecular_formula": "C51H79NO13",
    "inchikey": "QFJCIRLUMZQUOT-HPLJOQBZSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Sirolimus/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5284616",
    "retrieved_utc": "2026-09-10T05:33:47.078922+00:00",
    "clinical_note": "Oral macrocyclic immunosuppressant.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/5284616"
  },
  {
    "name": "Erythromycin",
    "archetype": "bRo5 Modalities",
    "smiles": "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O",
    "pubchem_cid": 12560,
    "molecular_formula": "C37H67NO13",
    "inchikey": "ULGZDMOVFRHVEP-RWJQBGPGSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Erythromycin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/12560",
    "retrieved_utc": "2026-09-10T05:33:49.385322+00:00",
    "clinical_note": "Oral macrolide antibacterial.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/12560"
  },
  {
    "name": "Azithromycin",
    "archetype": "bRo5 Modalities",
    "smiles": "CC[C@@H]1[C@@]([C@@H]([C@H](N(C[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)C)O)(C)O",
    "pubchem_cid": 447043,
    "molecular_formula": "C38H72N2O12",
    "inchikey": "MQTOSJVFKKJCRP-BICOPXKESA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Azithromycin/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/447043",
    "retrieved_utc": "2026-09-10T05:33:51.161699+00:00",
    "clinical_note": "Oral macrolide antibacterial.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/447043"
  },
  {
    "name": "Daclatasvir",
    "archetype": "bRo5 Modalities",
    "smiles": "CC(C)[C@@H](C(=O)N1CCC[C@H]1C2=NC=C(N2)C3=CC=C(C=C3)C4=CC=C(C=C4)C5=CN=C(N5)[C@@H]6CCCN6C(=O)[C@H](C(C)C)NC(=O)OC)NC(=O)OC",
    "pubchem_cid": 25154714,
    "molecular_formula": "C40H50N8O6",
    "inchikey": "FKRSSPOQAMALKA-CUPIEXAXSA-N",
    "structure_source": "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/Daclatasvir/property/IsomericSMILES,MolecularFormula,InChIKey/JSON",
    "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/25154714",
    "retrieved_utc": "2026-09-10T05:33:53.173974+00:00",
    "clinical_note": "Oral nonmacrocyclic antiviral; historical clinical reference, not a current availability claim.",
    "clinical_source": "https://pubchem.ncbi.nlm.nih.gov/compound/25154714"
  }
]'''


def low_desirability(x: float, ideal: float, poor: float) -> float:
    """One below ideal, zero above poor, linear interpolation in between."""
    if not (math.isfinite(x) and math.isfinite(ideal) and math.isfinite(poor)):
        raise ValueError("Desirability inputs must be finite")
    if ideal >= poor:
        raise ValueError("ideal must be below poor")
    return float(np.clip((poor - x) / (poor - ideal), 0.0, 1.0))


def high_desirability(x: float, poor: float, ideal: float) -> float:
    return 1.0 - low_desirability(x, poor, ideal)


def window_desirability(x: float, lo0: float, lo1: float,
                        hi1: float, hi0: float) -> float:
    if not lo0 < lo1 <= hi1 < hi0:
        raise ValueError("Invalid desirability window")
    return min(high_desirability(x, lo0, lo1), low_desirability(x, hi1, hi0))


def sigmoid(x: float) -> float:
    if x >= 0:
        return 1.0 / (1.0 + math.exp(-x))
    e = math.exp(x)
    return e / (1.0 + e)


def cns_mpo(clogp: float, logd: float, mw: float, tpsa: float,
            hbd: float, basic_pka: float | None) -> dict[str, float]:
    """Wager 2010 transforms; nonbasic molecules receive a pKa term of 1.

    The transforms are original, but RDKit/estimated inputs make the output
    an approximate-input implementation rather than the Pfizer calculator.
    """
    terms = {
        "mpo_d_clogp": low_desirability(clogp, 3, 5),
        "mpo_d_logd74": low_desirability(logd, 2, 4),
        "mpo_d_mw": low_desirability(mw, 360, 500),
        "mpo_d_tpsa": window_desirability(tpsa, 20, 40, 90, 120),
        "mpo_d_hbd": low_desirability(hbd, 0.5, 3.5),
        "mpo_d_basic_pka": 1.0 if basic_pka is None else low_desirability(basic_pka, 8, 10),
    }
    terms["cns_mpo_approx"] = sum(terms.values())
    # User-requested cLogP optimum 2-4 is a custom score, not Pfizer CNS-MPO.
    terms["custom_mpo_logp_2_4"] = (
        terms["cns_mpo_approx"] - terms["mpo_d_clogp"]
        + window_desirability(clogp, 0, 2, 4, 6)
    )
    return terms


def esol_log_s(clogp: float, mw: float, rotb: int, aromatic_fraction: float) -> float:
    """Delaney 2004 coefficients, log10(mol/L), with RDKit input descriptors.

    The aromatic correction is aromatic HEAVY ATOM fraction, not ring count.
    The rotatable-bond coefficient is 0.066 (not 0.66).
    """
    return 0.16 - 0.63 * clogp - 0.0062 * mw + 0.066 * rotb - 0.74 * aromatic_fraction


def molar_to_ug_ml(log_s: float, mw: float) -> float:
    # mol/L * g/mol * 1e6 ug/g / 1e3 mL/L
    return 10.0 ** log_s * mw * 1000.0


def parse_molecule(record: dict[str, Any]) -> Chem.Mol:
    mol = Chem.MolFromSmiles(record["smiles"])
    if mol is None or mol.GetNumHeavyAtoms() == 0:
        raise ValueError(f"Invalid SMILES: {record.get('name', '<unnamed>')}")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError(f"Use one parent active moiety, not salts/mixtures: {record['name']}")
    if Chem.GetFormalCharge(mol) != 0:
        raise ValueError(f"This benchmark requires neutral parent structures: {record['name']}")
    if any(a.GetFormalCharge() for a in mol.GetAtoms()):
        # Resonance-separated charge (e.g. nitro) is supported; permanent
        # zwitterions are not modelled by the simple neutral-fraction formula.
        charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() > 0]
        if any(a.GetAtomicNum() == 7 and a.GetTotalDegree() == 4 for a in charged):
            raise ValueError(f"Permanent zwitterion is outside ionization model: {record['name']}")
    expected_key = record.get("inchikey")
    if expected_key and Chem.MolToInchiKey(mol) != expected_key:
        raise ValueError(f"InChIKey mismatch: {record['name']}")
    expected_formula = record.get("molecular_formula")
    if expected_formula and rdMolDescriptors.CalcMolFormula(mol) != expected_formula:
        raise ValueError(f"Molecular formula mismatch: {record['name']}")
    return mol


def load_panel(path: Path | None = None) -> list[dict[str, Any]]:
    records = json.loads(path.read_text(encoding="utf-8") if path else PANEL_JSON)
    if not isinstance(records, list) or len(records) != 30:
        raise ValueError("Benchmark panel must contain exactly 30 records")
    names, keys = set(), set()
    for r in records:
        if not isinstance(r, dict) or not all(isinstance(r.get(k), str) and r[k]
                                            for k in ("name", "archetype", "smiles")):
            raise ValueError("Every record needs nonempty name, archetype, smiles strings")
        mol = parse_molecule(r)
        key = Chem.MolToInchiKey(mol)
        if r["name"] in names or key in keys:
            raise ValueError(f"Duplicate compound: {r['name']}")
        names.add(r["name"])
        keys.add(key)
        violations = sum((Descriptors.MolWt(mol) > 500, Crippen.MolLogP(mol) > 5,
                          Lipinski.NumHDonors(mol) > 5, Lipinski.NumHAcceptors(mol) > 10))
        if r["archetype"] == "Oral Drugs" and violations:
            raise ValueError(f"Oral reference is not strictly Ro5 compliant: {r['name']}")
        if r["archetype"] == "bRo5 Modalities" and not violations:
            raise ValueError(f"bRo5 reference has no Ro5 violation: {r['name']}")
    if Counter(r["archetype"] for r in records) != Counter({g: 10 for g in GROUPS}):
        raise ValueError("Panel requires exactly 10 compounds in each supported archetype")
    return records


# Deliberately simple chemical-class assumptions; these are NOT fitted pKa models.
# First match per site wins, preventing duplicate rules at one nitrogen.
BASE_RULES = (
    ("amidine_or_guanidine", "[NX2;H0,H1]=[CX3]([NX3])", 11.5),
    ("aliphatic_amine", "[NX3;!$(N-[C,S,P]=[O,S,N]);!$(N-a);!$(N-O);!$(N-N)]", 9.0),
    ("aniline_like", "[NX3;$(N-a);!$(N-[C,S,P]=[O,S,N])]", 5.0),
    ("aromatic_n", "[nH0;+0]", 4.5),
)
ACID_RULES = (
    ("carboxylic_acid", "[OX2H]-[CX3](=O)", 4.5),
    ("tetrazole", "[nH]1nnnc1", 4.9),
    ("sulfonamide", "[NX3;H1,H2]-S(=O)(=O)", 7.0),
    ("imide", "[NX3H]([CX3]=O)[CX3]=O", 8.5),
    ("phenol", "[OX2H]-[c]", 10.0),
    ("thiol", "[SX2H]", 9.5),
)


def ionization(mol: Chem.Mol, clogp: float,
               override: dict[str, Any] | None = None) -> dict[str, Any]:
    basic_sites: dict[int, tuple[str, float]] = {}
    acidic_sites: dict[int, tuple[str, float]] = {}
    for rules, sites in ((BASE_RULES, basic_sites), (ACID_RULES, acidic_sites)):
        for label, smarts, pka in rules:
            query = Chem.MolFromSmarts(smarts)
            for match in mol.GetSubstructMatches(query):
                sites.setdefault(match[0], (label, pka))
    base = max((p[1] for p in basic_sites.values()), default=None)
    acid = min((p[1] for p in acidic_sites.values()), default=None)
    provenance = "unvalidated_SMARTS_class_assumptions"
    if override:
        provenance = override["source"]
        base = override.get("basic_pka", base)
        acid = override.get("acidic_pka", acid)
    # Dominant independent acid/base approximation. No microstate enumeration,
    # coupled pKas, ion-pair partitioning, or charged-species lipid solubility.
    b = 0.0 if base is None else 10.0 ** (base - 7.4)
    a = 0.0 if acid is None else 10.0 ** (7.4 - acid)
    neutral = 1.0 / ((1 + b) * (1 + a))
    logd = clogp + math.log10(neutral)
    if override and "logd74" in override:
        logd = override["logd74"]
    return dict(basic_pka_assumed=base, acidic_pka_assumed=acid,
                basic_site_indices=list(basic_sites),
                basic_site_count=len(basic_sites), acidic_site_count=len(acidic_sites),
                ionization_rules=";".join(sorted({v[0] for v in [*basic_sites.values(), *acidic_sites.values()]})) or "no_matched_site",
                ionization_source=provenance, logd74_approx=logd,
                neutral_fraction_approx=neutral, cation_fraction_approx=b/(1+b),
                charge_burden_approx=b/(1+b)+a/(1+a))


def validate_overrides(raw: Any, names: set[str]) -> dict[str, Any]:
    if not isinstance(raw, dict):
        raise ValueError("Overrides must be an object keyed by compound name")
    for name, entry in raw.items():
        if name not in names or not isinstance(entry, dict):
            raise ValueError(f"Invalid override compound: {name}")
        if not isinstance(entry.get("source"), str) or not entry["source"].strip():
            raise ValueError(f"Override needs a source: {name}")
        if set(entry) - {"source", "basic_pka", "acidic_pka", "logd74"}:
            raise ValueError(f"Unknown override field for {name}")
        for k in ("basic_pka", "acidic_pka", "logd74"):
            if k not in entry or (entry[k] is None and k != "logd74"):
                continue
            v = entry[k]
            if isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v) or not -20 <= v <= 25:
                raise ValueError(f"Invalid {k} override for {name}")
    return raw


def empty_geometry(status: str) -> dict[str, Any]:
    return dict(geometry_status=status, embedding_method=None, conformers_embedded=0, conformers_converged=0,
                imhb_vectors_mean=None, imhb_vectors_min=None, imhb_vectors_max=None,
                unpaired_donor_vectors_mean=None, unpaired_acceptor_sites_mean=None,
                chi_geometry_range_proxy=None, absolute_chameleonic_hb_index=None,
                volume_3d_a3=None, esp_ring_proxy_e_per_a=None)


def geometry_calculation(smiles: str, nconf: int, seed: int, embed_timeout: int = 12) -> dict[str, Any]:
    """Gas-phase ETKDG/MMFF ensemble: geometric accessibility proxy only.

    A D-H vector is internally paired if D...A <=3.5 A, H...A <=2.6 A,
    angle D-H...A >=120 degrees, and graph distance D-A >=4 bonds.
    Each D-H vector counts at most once. An unpaired vector is NOT proven
    solvent-exposed. CHI proxy is the range of vector counts across conformers,
    NOT an absolute or experimentally validated chameleonicity index.
    """
    mol = Chem.AddHs(Chem.MolFromSmiles(smiles))
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.numThreads = 1
    params.pruneRmsThresh = 0.5
    params.maxIterations = 100
    if hasattr(params, "timeout"):
        params.timeout = embed_timeout
    AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=params)
    # Some RDKit releases return [-1] on timeout, rather than an empty list.
    # Actual attached conformers are the authoritative source of valid IDs.
    ids = [c.GetId() for c in mol.GetConformers()]
    method = "ETKDGv3"
    if not ids:
        params.useRandomCoords = True
        # Cyclic N-methyl peptides may require cis amides; an all-trans
        # embedding constraint can otherwise make valid structures impossible.
        params.forceTransAmides = False
        AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=params)
        ids = [c.GetId() for c in mol.GetConformers()]
        method = "ETKDGv3_random_coordinates_relaxed_amides_retry"
    result = empty_geometry("embedding_failed")
    result["embedding_method"] = method
    result["conformers_embedded"] = len(ids)
    if not ids:
        return result
    if not AllChem.MMFFHasAllMoleculeParams(mol):
        result["geometry_status"] = "MMFF_parameters_unavailable"
        return result
    optimized = AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=1, maxIters=800)
    usable = [cid for cid, (status, energy) in zip(ids, optimized)
              if status == 0 and math.isfinite(energy)]
    result["conformers_converged"] = len(usable)
    if not usable:
        result["geometry_status"] = "MMFF_not_converged"
        return result
    factory = ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir) / "BaseFeatures.fdef"))
    features = factory.GetFeaturesForMol(mol)
    donors = sorted({idx for f in features if f.GetFamily() == "Donor" for idx in f.GetAtomIds()})
    acceptors = sorted({idx for f in features if f.GetFamily() == "Acceptor" for idx in f.GetAtomIds()})
    vectors = [(d, h.GetIdx()) for d in donors for h in mol.GetAtomWithIdx(d).GetNeighbors() if h.GetAtomicNum() == 1]
    graph = Chem.GetDistanceMatrix(mol)
    counts, unpaired_a = [], []
    for cid in usable:
        xyz = mol.GetConformer(cid).GetPositions()
        internal = 0
        paired_acceptors = set()
        for donor, hydrogen in vectors:
            candidates = []
            for acceptor in acceptors:
                if graph[donor, acceptor] < 4:
                    continue
                dh, ah = xyz[donor] - xyz[hydrogen], xyz[acceptor] - xyz[hydrogen]
                dist_ha, dist_da = np.linalg.norm(ah), np.linalg.norm(xyz[donor] - xyz[acceptor])
                denominator = np.linalg.norm(dh) * dist_ha
                if denominator < 1e-10:
                    continue
                angle = math.degrees(math.acos(float(np.clip(np.dot(dh, ah)/denominator, -1, 1))))
                if dist_da <= 3.5 and dist_ha <= 2.6 and angle >= 120:
                    candidates.append((float(dist_ha), acceptor))
            if candidates:
                internal += 1
                paired_acceptors.add(min(candidates)[1])
        counts.append(internal)
        unpaired_a.append(len(acceptors) - len(paired_acceptors))
    result.update(geometry_status="ok" if len(usable) == len(ids) else "partial_convergence",
                  imhb_vectors_mean=float(np.mean(counts)), imhb_vectors_min=min(counts),
                  imhb_vectors_max=max(counts), unpaired_donor_vectors_mean=len(vectors)-float(np.mean(counts)),
                  unpaired_acceptor_sites_mean=float(np.mean(unpaired_a)),
                  chi_geometry_range_proxy=max(counts)-min(counts),
                  volume_3d_a3=float(AllChem.ComputeMolVolume(mol, confId=usable[0])))
    # Neutral-form Gasteiger point-charge potential at aromatic centroids.
    # This is a crude electrostatic descriptor, without receptor or protonation
    # microstates; it is NOT a quantum ESP or hERG binding calculation.
    try:
        AllChem.ComputeGasteigerCharges(mol, throwOnParamFailure=True)
        charges = np.array([float(a.GetProp("_GasteigerCharge")) for a in mol.GetAtoms()])
        xyz = mol.GetConformer(usable[0]).GetPositions()
        rings = [r for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
        if np.isfinite(charges).all():
            potentials = [float(np.sum(charges / np.maximum(np.linalg.norm(xyz - xyz[list(r)].mean(axis=0), axis=1), 1.0))) for r in rings]
            result["esp_ring_proxy_e_per_a"] = max(potentials, default=0.0)
    except (RuntimeError, ValueError):
        pass
    return result


def geometry_worker(connection: Any, smiles: str, nconf: int, seed: int, embed_timeout: int) -> None:
    try:
        connection.send(geometry_calculation(smiles, nconf, seed, embed_timeout))
    except Exception as exc:
        connection.send(empty_geometry(f"geometry_error:{type(exc).__name__}:{exc}"))
    finally:
        connection.close()


def bounded_geometry(smiles: str, nconf: int, seed: int, timeout: float) -> dict[str, Any]:
    if nconf == 0:
        return empty_geometry("disabled")
    context = mp.get_context("spawn")
    receive, send = context.Pipe(duplex=False)
    process = context.Process(target=geometry_worker,
                              args=(send, smiles, nconf, seed, max(1, int(timeout * .30))))
    process.start()
    send.close()
    try:
        if receive.poll(timeout):
            try:
                return receive.recv()
            except EOFError:
                return empty_geometry("worker_exited_without_result")
        return empty_geometry("timeout")
    finally:
        if process.is_alive():
            process.join(timeout=0.2)
        if process.is_alive():
            process.terminate()
        process.join(timeout=5)
        if process.is_alive():
            process.kill()
            process.join(timeout=5)
        receive.close()
        process.close()


def calculate_record(record: dict[str, Any], conformers: int = 6, seed: int = 20260910,
                     timeout: float = 45, override: dict[str, Any] | None = None) -> dict[str, Any]:
    mol = parse_molecule(record)
    mw, clogp = Descriptors.MolWt(mol), Crippen.MolLogP(mol)
    tpsa = rdMolDescriptors.CalcTPSA(mol)
    hbd, hba = Lipinski.NumHDonors(mol), Lipinski.NumHAcceptors(mol)
    rotb = rdMolDescriptors.CalcNumRotatableBonds(mol, rdMolDescriptors.NumRotatableBondsOptions.Strict)
    aromatic_rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    aromatic_fraction = sum(a.GetIsAromatic() for a in mol.GetAtoms()) / mol.GetNumHeavyAtoms()
    fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
    violations = sum((mw > 500, clogp > 5, hbd > 5, hba > 10))
    ion = ionization(mol, clogp, override)
    mpo = cns_mpo(clogp, ion["logd74_approx"], mw, tpsa, hbd, ion["basic_pka_assumed"])
    geometry = bounded_geometry(record["smiles"], conformers, seed, timeout)
    log_s = esol_log_s(clogp, mw, rotb, aromatic_fraction)
    volume = geometry["volume_3d_a3"]
    volume_method = "RDKit_grid_volume_first_converged_conformer"
    if volume is None:
        volume = 1.2 * mw
        volume_method = "unvalidated_1.2_times_MW_fallback"
    # Unfitted coefficients below define transparent ranking scenarios only.
    charge = ion["charge_burden_approx"]
    charge_volume = volume * (1 + 0.5 * charge)
    log_papp = 1.8 - 0.012 * tpsa - 0.35 * charge - 0.45 * math.log10(charge_volume / 300)
    papp = 10 ** log_papp  # numeric unit: 1e-6 cm/s
    # Desolvation-energy-like score in kcal/mol, NOT a computed free energy.
    dg_proxy = 0.025 * tpsa + 1.5 * charge + 0.15 * hbd - 0.30 * clogp
    hia = 100 * sigmoid((3.0 - dg_proxy) / 1.2)
    ring_atoms = [r for r in mol.GetRingInfo().AtomRings() if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    distance = Chem.GetDistanceMatrix(mol)
    separations = [min(float(distance[n, a]) for a in ring) for n in ion["basic_site_indices"] for ring in ring_atoms]
    matches = sum(3 <= d <= 8 for d in separations)
    topology = min(matches / 2.0, 1.0)
    esp = geometry["esp_ring_proxy_e_per_a"]
    esp_term = 0.0 if esp is None else float(np.clip(max(esp, 0) / 0.25, 0, 1))
    h_z = (-3.0 + 0.70 * (clogp - 2) + 1.5 * ion["cation_fraction_approx"]
           + 0.5 * min(aromatic_rings, 5) + 0.8 * topology + 0.25 * esp_term)
    herg = sigmoid(h_z)
    # Fixed anchors avoid cohort-dependent min/max normalization.
    oral_terms = {
        "oral_d_solubility": high_desirability(log_s, -7, -3),
        "oral_d_permeability": high_desirability(log_papp, -1, 1),
        "oral_d_herg_proxy": 1 - herg,
        "oral_d_polarity": low_desirability(tpsa, 90, 200),
        "oral_d_fsp3": high_desirability(fsp3, 0, 0.42),
        "oral_d_aromatic_rings": low_desirability(aromatic_rings, 3, 6),
    }
    oral_score = 100 * sum(oral_terms.values()) / len(oral_terms)
    # pKa +/-1 scenario envelope is sensitivity, NOT confidence/prediction bounds.
    sensitivity = []
    for bshift in (-1, 0, 1):
        for ashift in (-1, 0, 1):
            changed = {"source": "pKa_sensitivity_scenario",
                       "basic_pka": None if ion["basic_pka_assumed"] is None else ion["basic_pka_assumed"] + bshift,
                       "acidic_pka": None if ion["acidic_pka_assumed"] is None else ion["acidic_pka_assumed"] + ashift}
            if override and "logd74" in override:
                changed["logd74"] = override["logd74"]
            varied = ionization(mol, clogp, changed)
            sensitivity.append(cns_mpo(clogp, varied["logd74_approx"], mw, tpsa, hbd, varied["basic_pka_assumed"])["cns_mpo_approx"])
    out = dict(record)
    out.update(mw_da=mw, clogp_rdkit=clogp, tpsa_a2=tpsa, hbd=hbd, hba=hba,
               rotatable_bonds=rotb, aromatic_rings=aromatic_rings,
               aromatic_atom_fraction=aromatic_fraction, fsp3=fsp3,
               aromatic_rings_le3=aromatic_rings <= 3, fsp3_ge042=fsp3 >= 0.42,
               ro5_violations=violations, bro5_by_any_violation=bool(violations),
               largest_ring_atoms=max(map(len, mol.GetRingInfo().AtomRings()), default=0),
               unspecified_stereocenters=rdMolDescriptors.CalcNumUnspecifiedAtomStereoCenters(mol),
               **ion, **mpo, **geometry, **oral_terms,
               esol_logs_mol_l=log_s, solubility_ug_ml_esol=molar_to_ug_ml(log_s, mw),
               volume_used_a3=volume, volume_method=volume_method,
               charge_weighted_volume_proxy_a3=charge_volume,
               log10_caco2_papp_proxy=log_papp, caco2_papp_proxy_1e_minus6_cm_s=papp,
               desolvation_energy_proxy_kcal_mol=dg_proxy, hia_percent_proxy=hia,
               basic_n_aromatic_distance_min_bonds=min(separations, default=None),
               basic_n_aromatic_3to8_bond_pairs=matches, herg_topology_term=topology,
               herg_esp_term=esp_term, herg_esp_available=esp is not None,
               herg_liability_proxy=herg, oral_developability_score_proxy=oral_score,
               mpo_pka_sensitivity_min=min(sensitivity), mpo_pka_sensitivity_max=max(sensitivity),
               esol_domain_flag="extrapolation_caution" if mw > 500 or hba > 10 or hbd > 5 else "descriptor_substitution_unvalidated",
               admet_model_status="unvalidated_exploratory_surrogates",
               absolute_chi_status="not_identifiable_from_this_ensemble",
               clinical_outcome_used_in_scoring=False)
    return out


def plot_figures(df: pd.DataFrame, folder: Path, seed: int) -> None:
    folder.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid", context="notebook", font="DejaVu Sans")
    plt.rcParams.update({"figure.dpi": 120, "savefig.dpi": 300,
                         "axes.spines.top": False, "axes.spines.right": False,
                         "svg.hashsalt": "task1-mpo-v1"})
    def save(fig: Any, stem: str) -> None:
        fig.savefig(folder / (stem + ".png"), dpi=300, facecolor="white")
        fig.savefig(folder / (stem + ".svg"), facecolor="white", metadata={"Date": None})
        plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 6.6))
    sns.violinplot(data=df, x="archetype", y="cns_mpo_approx", order=GROUPS,
                   hue="archetype", hue_order=GROUPS, palette=PALETTE, legend=False,
                   inner=None, cut=0, linewidth=1, ax=ax)
    sns.boxplot(data=df, x="archetype", y="cns_mpo_approx", order=GROUPS,
                width=.16, color="white", showfliers=False,
                boxprops={"zorder": 3}, medianprops={"zorder": 4, "color": "#172B3A"},
                whiskerprops={"zorder": 4}, capprops={"zorder": 4}, ax=ax)
    rng = np.random.default_rng(seed)
    for index, group in enumerate(GROUPS):
        values = df.loc[df.archetype == group, "cns_mpo_approx"]
        ax.scatter(index + rng.uniform(-.11, .11, len(values)), values, s=24,
                   color="#172B3A", alpha=.8, zorder=4)
    ax.set(xlabel="", ylabel="CNS-MPO (0–6; approximate inputs)", ylim=(-.1, 6.35))
    ax.set_xticks(range(3), [f"{g}\n(n=10)" for g in GROUPS])
    ax.set_title("01  |  Developability across clinical archetypes", loc="left", pad=20, weight="bold")
    fig.text(.08, .035, "Wager 2010 transforms • RDKit cLogP + assumed ionization • descriptive panel, not a clinical classifier", fontsize=9)
    fig.subplots_adjust(left=.09, right=.97, bottom=.16, top=.88)
    save(fig, "fig1_mpo_distribution_violin")

    representatives = [n for n in ("Diazepam", "Terfenadine", "Venetoclax", "Cyclosporine A") if n in set(df.name)]
    if len(representatives) < 4:
        representatives = df.name.head(4).tolist()
    fig, axes = plt.subplots(1, 4, figsize=(15, 5.6), subplot_kw={"projection": "polar"})
    labels = ["Solubility", "Permeability*", "Safety\n(hERG only)*", "TPSA", "Fsp³", "CNS-MPO"]
    angles = np.linspace(0, 2 * np.pi, 6, endpoint=False)
    for ax, name in zip(axes, representatives):
        r = df.set_index("name").loc[name]
        values = [r.oral_d_solubility, r.oral_d_permeability, r.oral_d_herg_proxy,
                  r.oral_d_polarity, r.oral_d_fsp3, r.cns_mpo_approx / 6]
        color = PALETTE[GROUPS.index(r.archetype)]
        theta, values = np.r_[angles, angles[0]], np.r_[values, values[0]]
        ax.set_theta_offset(np.pi / 2)
        ax.set_theta_direction(-1)
        ax.plot(theta, values, color=color, linewidth=2)
        ax.fill(theta, values, color=color, alpha=.18)
        ax.set_xticks(angles, labels, fontsize=9)
        ax.tick_params(axis="x", pad=12)
        ax.set_ylim(0, 1)
        ax.set_yticks([.25, .5, .75, 1], ["", "0.5", "", "1.0"], fontsize=8, color="#697781")
        ax.set_title(name, y=1.27, fontsize=13, weight="bold")
    fig.suptitle("02  |  Fixed-anchor developability profiles", x=.035, y=.99, ha="left", weight="bold", fontsize=16)
    fig.text(.035, .035, "All axes: higher desirability is better. *Unvalidated proxies. Safety axis covers hERG liability only.\nSolubility uses log-molar ESOL; normalization anchors are documented in MODEL_CARD.md.", fontsize=10)
    fig.subplots_adjust(left=.06, right=.94, top=.76, bottom=.20, wspace=.7)
    save(fig, "fig2_admet_radar_profiles")

    fig, ax = plt.subplots(figsize=(11, 7.7))
    ymax = max(1350, float(df.mw_da.max()) + 140)
    xmin = min(-1.5, float(df.clogp_rdkit.min()) - .8)
    xmax = max(8.5, float(df.clogp_rdkit.max()) + 1)
    ax.fill_between([xmin, 5], 0, 500, color="#D8EEE4", alpha=.5)
    ax.axhline(500, ls="--", color="#668277", lw=1.2)
    ax.axvline(5, ls="--", color="#668277", lw=1.2)
    for group, marker in zip(GROUPS, ("o", "s", "^")):
        part = df[df.archetype == group]
        points = ax.scatter(part.clogp_rdkit, part.mw_da,
                            s=35 + 14 * part.rotatable_bonds, c=part.oral_developability_score_proxy,
                            cmap="viridis", vmin=0, vmax=100, marker=marker,
                            edgecolors="white", linewidths=.7, alpha=.92)
    offsets = {"Cyclosporine A": (-65, 15), "Venetoclax": (-85, -22), "ARV-110": (10, 14),
               "Paclitaxel": (-78, 16), "Terfenadine": (8, 12), "Diazepam": (8, -18),
               "Rifampicin": (-78, -24)}
    for name, offset in offsets.items():
        rows = df[df.name == name]
        if not rows.empty:
            r = rows.iloc[0]
            ax.annotate(name, (r.clogp_rdkit, r.mw_da), xytext=offset, textcoords="offset points",
                         fontsize=9, arrowprops={"arrowstyle": "-", "color": "#64727C", "lw": .7})
    ax.text(xmin + .15, 90, "Ro5 MW / cLogP region\nHBD and HBA constraints are not shown", fontsize=9, color="#356253")
    ax.set(xlim=(xmin, xmax), ylim=(0, ymax), xlabel="RDKit cLogP", ylabel="Molecular weight (Da)")
    ax.set_title("03  |  Chemical space beyond Rule of Five", loc="left", pad=20, weight="bold")
    handles = [Line2D([], [], marker=m, color="none", markerfacecolor="#6E7D88", markersize=8, label=g)
               for g, m in zip(GROUPS, ("o", "s", "^"))]
    legend = ax.legend(handles=handles, loc="upper left", frameon=True, fontsize=9)
    ax.add_artist(legend)
    size_handles = [ax.scatter([], [], s=35+14*n, facecolor="none", edgecolor="#697781", label=str(n)) for n in (0, 10, 20)]
    ax.legend(handles=size_handles, title="Rotatable bonds", loc="lower right", frameon=True, fontsize=8)
    bar = fig.colorbar(points, ax=ax, pad=.025, shrink=.83)
    bar.set_label("Oral developability score (0–100; exploratory)")
    fig.text(.085, .025, "bRo5 membership does not imply oral success. Paclitaxel is an IV comparator; ARV-110 is investigational.\nScores use fixed desirability functions and do not include clinical outcome labels.", fontsize=9)
    fig.subplots_adjust(left=.09, right=.92, bottom=.14, top=.9)
    save(fig, "fig3_bro5_chemical_space_landscape")


def atomic_text(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", newline="\n",
                                     dir=path.parent, delete=False, suffix=".tmp") as handle:
        temp = Path(handle.name)
        handle.write(content)
    try:
        temp.replace(path)
    finally:
        temp.unlink(missing_ok=True)


def canonical_bytes(value: Any) -> bytes:
    return json.dumps(value, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--output-dir", type=Path, default=Path(__file__).resolve().parent,
                        help="Parent for figures_task1/ and results_task1/ (default: this project's directory)")
    parser.add_argument("--panel", type=Path, help="Optional curated JSON replacement, exactly 10 compounds per archetype")
    parser.add_argument("--ionization-overrides", type=Path, help="JSON keyed by name, source required; see MODEL_CARD.md")
    parser.add_argument("--conformers", type=int, default=6, help="ETKDG conformers per molecule; 0 skips 3D (default 6)")
    parser.add_argument("--conformer-timeout", type=float, default=45, help="Wall seconds per isolated 3D worker (default 45)")
    parser.add_argument("--seed", type=int, default=20260910)
    parser.add_argument("--strict-3d", action="store_true", help="Fail if any compound has no converged conformer")
    args = parser.parse_args(argv)
    if not 0 <= args.conformers <= 100:
        parser.error("--conformers must be between 0 and 100")
    if not math.isfinite(args.conformer_timeout) or not 1 <= args.conformer_timeout <= 3600:
        parser.error("--conformer-timeout must be finite and between 1 and 3600 seconds")
    if not 0 <= args.seed <= 2147483647:
        parser.error("--seed must be an integer in [0, 2147483647]")
    if args.strict_3d and args.conformers == 0:
        parser.error("--strict-3d requires at least one conformer")
    logging.basicConfig(level=logging.INFO, format="%(levelname)s %(message)s")
    records = load_panel(args.panel)
    raw = json.loads(args.ionization_overrides.read_text(encoding="utf-8")) if args.ionization_overrides else {}
    overrides = validate_overrides(raw, {r["name"] for r in records})
    LOG.warning(WARNING)
    results = []
    for index, record in enumerate(records, 1):
        LOG.info("[%02d/30] %s", index, record["name"])
        result = calculate_record(record, args.conformers, args.seed, args.conformer_timeout, overrides.get(record["name"]))
        if result["geometry_status"] not in ("ok", "disabled"):
            LOG.warning("%s: %s", record["name"], result["geometry_status"])
        if args.strict_3d and not result["conformers_converged"]:
            raise RuntimeError(f"No converged 3D ensemble: {record['name']}")
        results.append(result)
    # JSON strictness catches accidental nonfinite results rather than silently
    # writing nonstandard NaN tokens. Unavailable scientific quantities are null.
    serialized = json.dumps(results, ensure_ascii=False, indent=2, allow_nan=False)
    df = pd.DataFrame(results)
    required = ["cns_mpo_approx", "solubility_ug_ml_esol", "caco2_papp_proxy_1e_minus6_cm_s",
                "hia_percent_proxy", "herg_liability_proxy", "oral_developability_score_proxy"]
    if not np.isfinite(df[required].to_numpy(dtype=float)).all():
        raise RuntimeError("Nonfinite required numerical results")
    base = args.output_dir.resolve()
    folder = base / "results_task1"
    folder.mkdir(parents=True, exist_ok=True)
    atomic_text(folder / "developability_results.json", serialized + "\n")
    csv_frame = df.copy()
    csv_frame["basic_site_indices"] = csv_frame.basic_site_indices.map(json.dumps)
    atomic_text(folder / "developability_results.csv", csv_frame.to_csv(index=False, float_format="%.8g", lineterminator="\n"))
    summary = df.groupby("archetype", sort=False).agg(
        n=("name", "count"), cns_mpo_median=("cns_mpo_approx", "median"),
        cns_mpo_mean=("cns_mpo_approx", "mean"), oral_score_median=("oral_developability_score_proxy", "median"),
        herg_proxy_median=("herg_liability_proxy", "median"), ro5_violation_median=("ro5_violations", "median"))
    atomic_text(folder / "group_summary.csv", summary.to_csv(float_format="%.8g", lineterminator="\n"))
    plot_figures(df, base / "figures_task1", args.seed)
    manifest = dict(version=VERSION, created_utc=datetime.now(timezone.utc).isoformat(),
                    python=sys.version, platform=platform.platform(),
                    dependencies={m.__name__: m.__version__ for m in (rdkit, np, pd, matplotlib, sns)},
                    seed=args.seed, conformers_requested=args.conformers,
                    conformer_timeout_seconds=args.conformer_timeout,
                    embedding_attempt_timeout_seconds=max(1, int(args.conformer_timeout * .30)),
                    panel_sha256=hashlib.sha256(canonical_bytes(records)).hexdigest(),
                    ionization_overrides=overrides,
                    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                    geometry_status_counts=dict(Counter(r["geometry_status"] for r in results)),
                    applicability_statement=WARNING,
                    score_semantics="original_CNS_MPO_transforms_with_approximate_inputs",
                    artifacts={})
    artifacts = [folder / "developability_results.json", folder / "developability_results.csv", folder / "group_summary.csv",
                 *sorted((base / "figures_task1").glob("fig[123]_*.png")),
                 *sorted((base / "figures_task1").glob("fig[123]_*.svg"))]
    manifest["artifacts"] = {p.relative_to(base).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest() for p in artifacts}
    atomic_text(folder / "run_manifest.json", json.dumps(manifest, indent=2, ensure_ascii=False, allow_nan=False)+"\n")
    LOG.info("Completed %s compounds. Outputs: %s", len(results), base)
    return 0


if __name__ == "__main__":
    mp.freeze_support()
    try:
        raise SystemExit(main())
    except (ValueError, OSError, RuntimeError) as error:
        LOG.error("Task failed: %s", error)
        raise SystemExit(1) from error
