"""Structural superposition and quantification of per-residue local deformation.

Whole-pocket topology is buried in noise, so we show whether changes localize at the clinical
mutation residues via (1) per-residue Hausdorff deformation and (2) local topology around the mutation.
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
    """Superpose the mobile structure onto ref using pocket Cα atoms and apply the rotation/translation to all of mobile."""
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
    """Per-residue WT<->mutant heavy-atom Hausdorff distance (Å) after superposition. Large at mutated/remodelled residues."""
    ref_res, mob_res, common, ca_rms, (rot, tran) = superpose(mobile_id, ref_id, resnums)
    out = {}
    for n in common:
        a = _heavy(ref_res[n])
        b = _heavy(mob_res[n], rot, tran)
        d = cdist(a, b)
        out[n] = float(max(d.min(axis=1).max(), d.min(axis=0).max()))  # symmetric Hausdorff
    return {"ca_rmsd": float(ca_rms), "hausdorff": out}


def local_pocket_cloud(chain_id: str, center_resnum: int, resnums, radius: float = 7.0):
    """Pocket heavy-atom point cloud within radius Å of the centre residue (for local topology)."""
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
    """WT<->mutant topological distance of the local pocket around the mutation site, relative to a noise baseline."""
    a = local_pocket_cloud(ref_id, center_resnum, resnums, radius)
    b = local_pocket_cloud(mutant_id, center_resnum, resnums, radius)
    if a is None or b is None or len(a) < 5 or len(b) < 5:
        return None
    dist = topo.total_wdist(topo.diagram(a), topo.diagram(b), dims)
    base = topo.noise_baseline(a, sigma=0.3, n=15, dims=dims)
    return {"distance": dist, "noise_p95": base["p95"], "noise_mean": base["mean"],
            "ratio_vs_noise": dist / (base["p95"] + 1e-9), "n_atoms": int(len(a))}
