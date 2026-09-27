"""Full analysis: structural and topological signals of the EGFR osimertinib resistance cascade.

Key findings:
- The binding-pocket backbone is conserved (Cα RMSD ~1 Å); deformation is localized at the clinical mutation residues.
- Deformation magnitude separates resistance types: L858R (steric) >> T790M (steric) >> C797S (nearly silent).
- The structural "silence" of C797S is itself the finding: osimertinib resistance comes from chemistry
  (loss of the covalent anchor), not shape, so structural/topological descriptors are blind to this type
  -> structural and chemical features must be used together.
- Note: a topological distance computed over the whole pocket is buried in coordinate noise
  (localized, resolution-appropriate descriptors are needed).
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from . import data
from . import geometry as geo
from . import topology as topo

RESULTS = Path(__file__).resolve().parent.parent / "results"
KEY = {790: "T790M\n(gatekeeper)", 797: "C797S\n(covalent anchor)", 858: "L858R\n(activation loop)"}
# Structure that most cleanly carries each mutation (vs WT 6JXT)
SITE_STRUCT = {858: ("6JWL", "L858R"), 790: ("6JX0", "T790M"), 797: ("6LUD", "C797S")}


def run(dims=(1, 2)):
    RESULTS.mkdir(exist_ok=True)
    resnums = data.pocket_residue_numbers()
    series = data.OSIMERTINIB_SERIES
    ids = [v["id"] for v in series]
    labels = {v["id"]: v["label"] for v in series}
    wt = ids[0]

    clouds = {vid: data.pocket_cloud(vid, resnums) for vid in ids}
    dgms = {vid: topo.diagram(clouds[vid]) for vid in ids}

    # --- (A) Whole-pocket topology vs noise (an honest methodological result) ---
    baseline = topo.noise_baseline(clouds[wt], sigma=0.4, n=20, dims=dims)
    global_dist = {vid: topo.total_wdist(dgms[wt], dgms[vid], dims) for vid in ids[1:]}

    # --- (B) Per-residue Hausdorff deformation: localization (WT vs triple mutant 6LUD) ---
    haus = geo.per_residue_hausdorff("6LUD", wt, resnums)
    nonkey = [v for n, v in haus["hausdorff"].items() if n not in KEY]
    struct_background = float(np.median(nonkey))

    # --- (C) Per-mutation site deformation: grading resistance types ---
    site_footprint = {}
    for site, (sid, name) in SITE_STRUCT.items():
        h = geo.per_residue_hausdorff(sid, wt, resnums)
        site_footprint[name] = {"site": site, "struct": sid,
                                "hausdorff": h["hausdorff"].get(site),
                                "ca_rmsd": h["ca_rmsd"]}

    results = {
        "pocket_residues": resnums,
        "n_pocket_atoms": {vid: int(len(clouds[vid])) for vid in ids},
        "global_topology": {
            "noise_baseline": {k: baseline[k] for k in ("mean", "std", "p95")},
            "distance_from_WT": {labels[v]: global_dist[v] for v in ids[1:]},
            "note": "whole-pocket Wasserstein < noise baseline -> localized descriptors needed",
        },
        "localization_pocket_CA_RMSD": haus["ca_rmsd"],
        "structural_background_hausdorff_median": struct_background,
        "site_footprint": site_footprint,
    }
    (RESULTS / "results.json").write_text(json.dumps(results, indent=2, ensure_ascii=False))

    _fig_footprint(site_footprint, struct_background)
    _fig_localization(haus["hausdorff"])
    _fig_global_noise(global_dist, labels, baseline, ids)
    _fig_persistence(dgms, labels, ids)
    return results


def _fig_footprint(site_footprint, background):
    order = ["L858R", "T790M", "C797S"]
    vals = [site_footprint[k]["hausdorff"] for k in order]
    colors = ["#c0392b", "#e67e22", "#95a5a6"]
    fig, ax = plt.subplots(figsize=(5.2, 4))
    ax.bar(order, vals, color=colors)
    ax.axhline(background, color="grey", ls="--", lw=1.2,
               label=f"structural background (median {background:.1f} Å)")
    for i, v in enumerate(vals):
        ax.text(i, v + 0.15, f"{v:.1f} Å", ha="center", fontsize=9, weight="bold")
    ax.set_ylabel("Footprint at mutated residue (Hausdorff, Å)", fontsize=9)
    ax.set_title("Resistance mutations graded by structural footprint\n"
                 "C797S near-silent: its resistance is chemical, not structural", fontsize=10)
    ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(RESULTS / "fig_site_footprint.png", dpi=150); plt.close(fig)


def _fig_localization(hmap):
    items = sorted(hmap.items())
    xs = [n for n, _ in items]; ys = [v for _, v in items]
    colors = ["#c0392b" if n in KEY else "#3b6ea5" for n in xs]
    fig, ax = plt.subplots(figsize=(8, 3.6))
    ax.bar([str(n) for n in xs], ys, color=colors)
    ax.set_ylabel("Hausdorff deformation vs WT (Å)")
    ax.set_title("Per-residue deformation, WT vs L858R/T790M/C797S "
                 "(red = clinical mutation sites)")
    ax.tick_params(axis="x", labelsize=6, rotation=90)
    fig.tight_layout(); fig.savefig(RESULTS / "fig_localization.png", dpi=150); plt.close(fig)


def _fig_global_noise(global_dist, labels, baseline, ids):
    vids = ids[1:]
    vals = [global_dist[v] for v in vids]
    fig, ax = plt.subplots(figsize=(5, 3.8))
    ax.bar([labels[v] for v in vids], vals, color="#3b6ea5")
    ax.axhspan(0, baseline["p95"], color="grey", alpha=0.3,
               label=f"coordinate-noise floor (p95 {baseline['p95']:.0f})")
    ax.set_ylabel("Whole-pocket topological distance (H1+H2)")
    ax.set_title("Naive whole-pocket topology is noise-dominated\n(motivates localized descriptors)")
    ax.legend(fontsize=8); plt.xticks(rotation=15, fontsize=7)
    fig.tight_layout(); fig.savefig(RESULTS / "fig_global_noise.png", dpi=150); plt.close(fig)


def _fig_persistence(dgms, labels, ids):
    fig, axes = plt.subplots(1, len(ids), figsize=(3.1 * len(ids), 3.1), sharex=True, sharey=True)
    for ax, vid in zip(axes, ids):
        h1 = dgms[vid][1]
        if len(h1):
            ax.scatter(h1[:, 0], h1[:, 1], s=9, alpha=0.6, color="#2c7fb8")
        ax.plot([0, 12], [0, 12], color="grey", lw=0.7)
        ax.set_title(labels[vid], fontsize=8); ax.set_xlabel("birth", fontsize=8)
    axes[0].set_ylabel("death (H1)", fontsize=8)
    fig.suptitle("Binding-pocket H1 persistence diagrams", fontsize=10)
    fig.tight_layout(); fig.savefig(RESULTS / "fig_persistence_H1.png", dpi=150); plt.close(fig)


if __name__ == "__main__":
    r = run()
    print("=== Localization: pocket Cα RMSD (WT vs 6LUD) = %.2f Å ===" % r["localization_pocket_CA_RMSD"])
    print("Structural background (median Hausdorff of non-mutated residues) = %.2f Å" % r["structural_background_hausdorff_median"])
    print("=== Per-mutation site deformation (resistance-type grading) ===")
    for name in ["L858R", "T790M", "C797S"]:
        f = r["site_footprint"][name]
        print(f"  {name} (res {f['site']}, {f['struct']}): {f['hausdorff']:.2f} Å")
    print("=== Whole-pocket topology vs noise (honest methodological result) ===")
    b = r["global_topology"]["noise_baseline"]
    print(f"  noise p95 = {b['p95']:.1f}")
    for lab, d in r["global_topology"]["distance_from_WT"].items():
        print(f"  {lab}: {d:.1f}  ({'BELOW' if d < b['p95'] else 'above'} noise)")
