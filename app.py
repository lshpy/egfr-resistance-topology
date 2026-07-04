"""설명 가능한 EGFR 내성 구조 분석 데모 (Gradio).

돌연변이를 고르면, 그 돌연변이가 오시머티닙 결합 포켓에 남기는
구조적 발자국(Hausdorff)과 포켓 위상 다이어그램, 그리고 '입체 vs 화학' 해석을 보여준다.
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
    "L858R (활성화 루프, 민감화)": ("6JWL", 858),
    "T790M (게이트키퍼, 1·2세대 내성)": ("6JX0", 790),
    "L858R/T790M/C797S (오시머티닙 내성)": ("6LUD", 797),
}
_resnums = data.pocket_residue_numbers()
_INTERP = {
    790: ("입체(steric)", "게이트키퍼 Thr→Met로 곁사슬이 커져 약물 결합 자리가 좁아진다. 구조 변화가 뚜렷해 구조/위상 기술로 감지된다."),
    858: ("입체(steric)", "활성화 루프 Leu→Arg로 크고 하전된 곁사슬이 들어와 큰 구조 재배열을 만든다. 발자국이 가장 크다."),
    797: ("화학(chemical)", "Cys→Ser는 거의 같은 부피(isosteric)라 모양이 거의 안 바뀐다. 오시머티닙이 C797과 맺던 공유결합 앵커가 사라져 내성이 생기며, 이 유형은 구조/위상만으로는 보이지 않는다."),
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

    verdict = "구조로 감지됨 ✅" if fp > background else "구조적으로 침묵 — 화학적 내성 ⚠️"
    md = (
        f"## {choice}\n"
        f"- 돌연변이 자리(res {site}) 구조 발자국: **{fp:.2f} Å**\n"
        f"- 포켓 구조적 배경(중앙값): {background:.2f} Å → **{verdict}**\n"
        f"- 내성 유형: **{kind}**\n\n"
        f"> {why}\n\n"
        f"> 구조 PDB: {struct_id} (오시머티닙 결합) · WT: {WT}. "
        f"예비·비임상 연구 데모입니다."
    )
    return md, fig


with gr.Blocks(title="EGFR 내성 구조·위상 분석") as demo:
    gr.Markdown(
        "# 🧬 설명 가능한 EGFR 표적치료 내성 분석\n"
        "폐암 EGFR 표적치료(오시머티닙)의 내성 돌연변이가 결합 포켓에 남기는 "
        "**구조적 발자국**과 **포켓 위상**을 정량화합니다. 같은 약물이 결합한 구조들만 비교해 "
        "리간드 교란을 통제했습니다. (공개 PDB·지속 호몰로지·예비연구)"
    )
    choice = gr.Radio(list(CHOICES), value=list(CHOICES)[2], label="내성 돌연변이 선택")
    btn = gr.Button("구조·위상 분석", variant="primary")
    out_md = gr.Markdown()
    out_fig = gr.Plot()
    btn.click(analyze, choice, [out_md, out_fig])

if __name__ == "__main__":
    demo.launch()
