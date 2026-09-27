"""Persistent homology of the binding pocket and topological distances between variants.

- Computes H0/H1/H2 persistence diagrams of the pocket heavy-atom point cloud with a Vietoris–Rips complex.
- Quantifies topological differences between variant structures with the Wasserstein distance.
- Estimates a noise baseline at the scale of crystallographic coordinate error by bootstrap (honesty check).
- Explains "which residue drives the topological change" via per-residue leave-one-out attribution (XAI).
"""
from __future__ import annotations

import numpy as np
from persim import wasserstein
from ripser import ripser

MAXDIM = 2


def diagram(cloud: np.ndarray, maxdim: int = MAXDIM):
    """Point cloud -> [H0, H1, H2] persistence diagrams. Infinite bars (the inf in H0) are removed."""
    dgms = ripser(cloud, maxdim=maxdim)["dgms"]
    cleaned = []
    for d in dgms:
        finite = d[np.isfinite(d).all(axis=1)]
        cleaned.append(finite)
    return cleaned


def wdist(dgm_a, dgm_b, dim: int) -> float:
    return float(wasserstein(dgm_a[dim], dgm_b[dim]))


def total_wdist(dgm_a, dgm_b, dims=(1, 2)) -> float:
    return float(sum(wdist(dgm_a, dgm_b, d) for d in dims))


def noise_baseline(cloud: np.ndarray, sigma: float = 0.4, n: int = 20,
                   dims=(1, 2), seed: int = 0) -> dict:
    """Self-vs-self Wasserstein distribution obtained by adding sigma (Å) Gaussian noise to the coordinates.
    Baseline for the distance that "no difference" produces at the scale of crystal-structure coordinate uncertainty."""
    rng = np.random.default_rng(seed)
    base = diagram(cloud)
    ds = []
    for _ in range(n):
        jittered = cloud + rng.normal(0, sigma, size=cloud.shape)
        ds.append(total_wdist(base, diagram(jittered), dims))
    ds = np.array(ds)
    return {"mean": float(ds.mean()), "std": float(ds.std()),
            "p95": float(np.percentile(ds, 95)), "samples": ds}


def residue_attribution(atoms_a: dict, atoms_b: dict, dims=(1, 2)) -> dict[int, float]:
    """How much the topological distance drops when each residue is removed from both structures (attribution).
    Larger values mean that residue drives the topological difference between the two structures."""
    resnums = sorted(set(atoms_a) & set(atoms_b))
    full_a = np.vstack([atoms_a[n] for n in resnums])
    full_b = np.vstack([atoms_b[n] for n in resnums])
    base = total_wdist(diagram(full_a), diagram(full_b), dims)
    attr = {}
    for r in resnums:
        keep = [n for n in resnums if n != r]
        ca = np.vstack([atoms_a[n] for n in keep])
        cb = np.vstack([atoms_b[n] for n in keep])
        attr[r] = base - total_wdist(diagram(ca), diagram(cb), dims)
    return attr
