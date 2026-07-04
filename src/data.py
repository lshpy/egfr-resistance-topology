"""EGFR 키나아제 도메인 구조 — 오시머티닙 내성 캐스케이드.

같은 약물(osimertinib, 리간드 YY3)이 결합한 채 유전형만 다른 구조들을 사용해
'리간드 유도적합' 교란을 통제한다. 지속 호몰로지는 강체운동에 불변이므로
중첩 없이도 결합 포켓의 위상을 비교할 수 있다.

출처: RCSB PDB (https://www.rcsb.org). 공개 구조.
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

# 오시머티닙(YY3) 결합 캐스케이드 — 임상 순서: WT → 민감화 → 게이트키퍼 → 공유앵커 소실
OSIMERTINIB_SERIES = [
    {"id": "6JXT", "label": "WT", "mutations": [], "role": "야생형(약물 결합)"},
    {"id": "6JWL", "label": "L858R", "mutations": ["L858R"], "role": "민감화 돌연변이"},
    {"id": "6JX0", "label": "T790M", "mutations": ["T790M"], "role": "게이트키퍼(1·2세대 내성)"},
    {"id": "6LUD", "label": "L858R/T790M/C797S", "mutations": ["L858R", "T790M", "C797S"],
     "role": "★오시머티닙 내성(공유앵커 C797 소실)"},
]
# 1세대(gefitinib, IRE) 통제쌍 — 보조 분석(게이트키퍼 내성)
GEFITINIB_SERIES = [
    {"id": "2ITY", "label": "WT", "mutations": [], "role": "야생형(gefitinib)"},
    {"id": "4I22", "label": "L858R/T790M", "mutations": ["L858R", "T790M"], "role": "게이트키퍼 내성(gefitinib)"},
]

LIGAND = {"6JXT": "YY3", "6JWL": "YY3", "6JX0": "YY3", "6LUD": "YY3",
          "2ITY": "IRE", "4I22": "IRE"}

# 결합 포켓을 정의하는 기준 구조/컷오프
POCKET_REFERENCE = "6JXT"
POCKET_CUTOFF = 8.0  # Å, 리간드 원자로부터

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
    """WT 기준 구조에서 리간드 8Å 이내 잔기 번호(고정 포켓 정의)."""
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
    """잔기번호 → 그 잔기 중원자 좌표(N,3). 고정 포켓 잔기집합을 모든 구조에 적용."""
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
