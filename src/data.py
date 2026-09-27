"""EGFR kinase-domain structures: the osimertinib resistance cascade.

Uses structures bound to the same drug (osimertinib, ligand YY3) that differ only in genotype,
to control for ligand-induced-fit confounding. Persistent homology is invariant to rigid motion,
so binding-pocket topology can be compared without superposition.

Source: RCSB PDB (https://www.rcsb.org). Public structures.
"""
from __future__ import annotations

import os
import urllib.request
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore")
from Bio.PDB import PDBParser  # noqa: E402

DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "pdb"

# Osimertinib (YY3)-bound cascade, in clinical order: WT -> sensitizing -> gatekeeper -> loss of covalent anchor
OSIMERTINIB_SERIES = [
    {"id": "6JXT", "label": "WT", "mutations": [], "role": "wild type (drug-bound)"},
    {"id": "6JWL", "label": "L858R", "mutations": ["L858R"], "role": "sensitizing mutation"},
    {"id": "6JX0", "label": "T790M", "mutations": ["T790M"], "role": "gatekeeper (1st/2nd-gen resistance)"},
    {"id": "6LUD", "label": "L858R/T790M/C797S", "mutations": ["L858R", "T790M", "C797S"],
     "role": "★osimertinib resistance (loss of covalent anchor C797)"},
]
# 1st-generation (gefitinib, IRE) control pair: auxiliary analysis (gatekeeper resistance)
GEFITINIB_SERIES = [
    {"id": "2ITY", "label": "WT", "mutations": [], "role": "wild type (gefitinib)"},
    {"id": "4I22", "label": "L858R/T790M", "mutations": ["L858R", "T790M"], "role": "gatekeeper resistance (gefitinib)"},
]

LIGAND = {"6JXT": "YY3", "6JWL": "YY3", "6JX0": "YY3", "6LUD": "YY3",
          "2ITY": "IRE", "4I22": "IRE"}

# Reference structure / cutoff that define the binding pocket
POCKET_REFERENCE = "6JXT"
POCKET_CUTOFF = 8.0  # Å, from ligand atoms

_SKIP_HET = {"HOH", "NAG", "SO4", "GOL", "EDO", "CL", "NA", "MG", "ACT", "PO4", "MES", "NO3"}
_parser = PDBParser(QUIET=True)


def fetch(pdb_id: str) -> Path:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    path = DATA_DIR / f"{pdb_id}.pdb"
    if not path.exists():
        urllib.request.urlretrieve(
            f"https://files.rcsb.org/download/{pdb_id}.pdb", path
        )
    return path


def _first_polymer_chain(model):
    for ch in model:
        if any(r.id[0] == " " for r in ch):
            return ch
    raise ValueError("no polymer chain")


def load_chain(pdb_id: str):
    model = _parser.get_structure(pdb_id, fetch(pdb_id))[0]
    return _first_polymer_chain(model)


def ligand_atoms(chain, pdb_id: str):
    lig = LIGAND[pdb_id]
    return [a for r in chain if r.resname == lig for a in r if a.element != "H"]


def pocket_residue_numbers() -> list[int]:
    """Residue numbers within 8 Å of the ligand in the WT reference structure (fixed pocket definition)."""
    chain = load_chain(POCKET_REFERENCE)
    lig_xyz = np.array([a.coord for a in ligand_atoms(chain, POCKET_REFERENCE)])
    nums = set()
    for r in chain:
        if r.id[0] != " ":
            continue
        for a in r:
            if a.element == "H":
                continue
            if np.linalg.norm(lig_xyz - a.coord, axis=1).min() < POCKET_CUTOFF:
                nums.add(r.id[1])
                break
    return sorted(nums)


def pocket_atoms_by_residue(pdb_id: str, resnums: list[int]) -> dict[int, np.ndarray]:
    """Residue number -> heavy-atom coordinates (N,3) of that residue. The fixed pocket residue set is applied to every structure."""
    chain = load_chain(pdb_id)
    wanted = set(resnums)
    out: dict[int, np.ndarray] = {}
    for r in chain:
        if r.id[0] == " " and r.id[1] in wanted:
            xyz = [a.coord for a in r if a.element != "H"]
            if xyz:
                out[r.id[1]] = np.array(xyz, dtype=float)
    return out


def pocket_cloud(pdb_id: str, resnums: list[int]) -> np.ndarray:
    d = pocket_atoms_by_residue(pdb_id, resnums)
    return np.vstack([d[n] for n in sorted(d)])
