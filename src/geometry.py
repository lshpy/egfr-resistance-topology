"""구조 중첩(superposition)과 잔기별 국소 변형 정량화.

전체 포켓 위상은 잡음에 묻히므로, 변화가 임상 돌연변이 잔기에 국소화되는지를
(1) 잔기별 Hausdorff 변형과 (2) 돌연변이 주변 국소 위상으로 보인다.
"""
from __future__ import annotations

import numpy as np
from Bio.PDB import Superimposer
from scipy.spatial.distance import cdist

from . import data
from . import topology as topo


def _residue_map(chain, resnums):
    wanted = set(resnums)
    return {r.id[1]: r for r in chain if r.id[0] == " " and r.id[1] in wanted}


def superpose(mobile_id: str, ref_id: str, resnums):
    """mobile 구조를 ref에 포켓 Cα로 중첩하고 회전·이동을 mobile 전체에 적용."""
    ref_chain = data.load_chain(ref_id)
    mob_chain = data.load_chain(mobile_id)
    ref_res, mob_res = _residue_map(ref_chain, resnums), _residue_map(mob_chain, resnums)
    common = [n for n in sorted(resnums) if n in ref_res and n in mob_res
              and "CA" in ref_res[n] and "CA" in mob_res[n]]
    sup = Superimposer()
    sup.set_atoms([ref_res[n]["CA"] for n in common], [mob_res[n]["CA"] for n in common])
    rot, tran = sup.rotran
    return ref_res, mob_res, common, sup.rms, (rot, tran)


def _heavy(res, rot=None, tran=None):
    xyz = np.array([a.coord for a in res if a.element != "H"], dtype=float)
    if rot is not None:
        xyz = xyz @ rot + tran
    return xyz


def per_residue_hausdorff(mobile_id: str, ref_id: str, resnums):
    """중첩 후 잔기별 WT↔변이 중원자 Hausdorff 거리(Å). 돌연변이·리모델링 잔기에서 커짐."""
    ref_res, mob_res, common, ca_rms, (rot, tran) = superpose(mobile_id, ref_id, resnums)
    out = {}
    for n in common:
        a = _heavy(ref_res[n])
        b = _heavy(mob_res[n], rot, tran)
        d = cdist(a, b)
        out[n] = float(max(d.min(axis=1).max(), d.min(axis=0).max()))  # 대칭 Hausdorff
    return {"ca_rmsd": float(ca_rms), "hausdorff": out}


def local_pocket_cloud(chain_id: str, center_resnum: int, resnums, radius: float = 7.0):
    """중심 잔기 주변 radius Å 안의 포켓 중원자 점구름(국소 위상용)."""
    chain = data.load_chain(chain_id)
    rmap = _residue_map(chain, resnums)
    if center_resnum not in rmap:
        return None
    center = _heavy(rmap[center_resnum])
    pts = []
    for n, r in rmap.items():
        xyz = _heavy(r)
        if cdist(center, xyz).min() < radius:
            pts.extend(xyz.tolist())
    return np.array(pts)


def local_topology_shift(mutant_id: str, ref_id: str, center_resnum: int, resnums,
                         radius: float = 7.0, dims=(1, 2)):
    """돌연변이 자리 주변 국소 포켓의 WT↔변이 위상 거리 + 잡음 기준선 대비."""
    a = local_pocket_cloud(ref_id, center_resnum, resnums, radius)
    b = local_pocket_cloud(mutant_id, center_resnum, resnums, radius)
    if a is None or b is None or len(a) < 5 or len(b) < 5:
        return None
    dist = topo.total_wdist(topo.diagram(a), topo.diagram(b), dims)
    base = topo.noise_baseline(a, sigma=0.3, n=15, dims=dims)
    return {"distance": dist, "noise_p95": base["p95"], "noise_mean": base["mean"],
            "ratio_vs_noise": dist / (base["p95"] + 1e-9), "n_atoms": int(len(a))}
