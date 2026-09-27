"""Explainable EGFR resistance structure analysis demo (Gradio).

Pick a mutation to see the structural footprint (Hausdorff) it leaves in the osimertinib binding pocket,
the pocket persistence diagram, and a "steric vs chemical" interpretation.
"""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import gradio as gr

from src import data
from src import geometry as geo
from src import topology as topo

WT = "6JXT"
CHOICES = {
    "L858R (activation loop, sensitizing)": ("6JWL", 858),
    "T790M (gatekeeper, 1st/2nd-gen resistance)": ("6JX0", 790),
    "L858R/T790M/C797S (osimertinib resistance)": ("6LUD", 797),
}
_resnums = data.pocket_residue_numbers()
_INTERP = {
    790: ("steric", "The gatekeeper Thr→Met substitution enlarges the side chain and narrows the drug-binding site. The structural change is clear, so structural/topological descriptors detect it."),
    858: ("steric", "The activation-loop Leu→Arg substitution introduces a large, charged side chain that causes a major structural rearrangement. It has the largest footprint."),
    797: ("chemical", "Cys→Ser is nearly isosteric, so the shape barely changes. Resistance arises because the covalent anchor osimertinib forms with C797 is lost; this type is invisible to structure/topology alone."),
}


def analyze(choice):
    struct_id, site = CHOICES[choice]
    haus = geo.per_residue_hausdorff(struct_id, WT, _resnums)
    fp = haus["hausdorff"].get(site, float("nan"))
    background = float(np.median([v for n, v in haus["hausdorff"].items() if n not in (790, 797, 858)]))
    kind, why = _INTERP[site]

    wt_cloud = data.pocket_cloud(WT, _resnums)
    mut_cloud = data.pocket_cloud(struct_id, _resnums)
    dgm_wt, dgm_mut = topo.diagram(wt_cloud), topo.diagram(mut_cloud)

    fig, ax = plt.subplots(figsize=(4.2, 4.2))
    for dg, c, lab in [(dgm_wt, "#7f8c8d", "WT"), (dgm_mut, "#c0392b", choice.split(" ")[0])]:
        h1 = dg[1]
        if len(h1):
            ax.scatter(h1[:, 0], h1[:, 1], s=14, alpha=0.6, color=c, label=lab)
    ax.plot([0, 12], [0, 12], color="grey", lw=0.7)
    ax.set_xlabel("birth"); ax.set_ylabel("death"); ax.set_title("Pocket H1 persistence"); ax.legend()
    fig.tight_layout()

    verdict = "detected structurally ✅" if fp > background else "structurally silent — chemical resistance ⚠️"
    md = (
        f"## {choice}\n"
        f"- Structural footprint at the mutation site (res {site}): **{fp:.2f} Å**\n"
        f"- Pocket structural background (median): {background:.2f} Å → **{verdict}**\n"
        f"- Resistance type: **{kind}**\n\n"
        f"> {why}\n\n"
        f"> Structure PDB: {struct_id} (osimertinib-bound) · WT: {WT}. "
        f"Preliminary, non-clinical research demo."
    )
    return md, fig


with gr.Blocks(title="EGFR resistance structure & topology analysis") as demo:
    gr.Markdown(
        "# 🧬 Explainable EGFR targeted-therapy resistance analysis\n"
        "Quantifies the **structural footprint** and **pocket topology** that resistance mutations to "
        "EGFR targeted therapy in lung cancer (osimertinib) leave in the binding pocket. Only structures bound "
        "to the same drug are compared, controlling for ligand confounding. (Public PDB · persistent homology · preliminary study)"
    )
    choice = gr.Radio(list(CHOICES), value=list(CHOICES)[2], label="Select resistance mutation")
    btn = gr.Button("Analyze structure & topology", variant="primary")
    out_md = gr.Markdown()
    out_fig = gr.Plot()
    btn.click(analyze, choice, [out_md, out_fig])

if __name__ == "__main__":
    demo.launch()
