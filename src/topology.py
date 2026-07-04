"""결합 포켓의 지속 호몰로지(persistent homology)와 변이 간 위상 거리.

- Vietoris–Rips 복합체로 포켓 중원자 점구름의 H0/H1/H2 지속 다이어그램 계산.
- 변이 구조 사이의 위상 차이를 Wasserstein 거리로 정량화.
- 결정학적 좌표 오차 규모의 잡음 기준선을 부트스트랩으로 추정(정직성).
- 잔기별 leave-one-out 기여도로 '어느 잔기가 위상 변화를 만드는지' 설명(XAI).
"""
from __future__ import annotations

import numpy as np
from persim import wasserstein
from ripser import ripser

MAXDIM = 2


def diagram(cloud: np.ndarray, maxdim: int = MAXDIM):
    """점구름 → [H0, H1, H2] 지속 다이어그램. 무한 지속(H0의 inf)은 제거."""
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
    """좌표에 sigma(Å) 가우시안 잡음을 넣어 얻는 자기-자기 Wasserstein 분포.
    결정 구조의 좌표 불확실성 규모에서 '차이 없음'이 만드는 거리의 기준선."""
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
    """각 잔기를 양쪽에서 제거했을 때 위상 거리가 얼마나 줄어드는지(기여도).
    값이 클수록 그 잔기가 두 구조의 위상 차이를 주도함."""
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
