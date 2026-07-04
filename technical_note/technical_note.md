# Structural Footprints of Targeted-Therapy Resistance Mutations in the EGFR Kinase Pocket: A Topology-Aware Preliminary Study

**Author:** Seunghyun Lee — ORCID [0009-0006-1926-653X](https://orcid.org/0009-0006-1926-653X)
**Type:** Preliminary technical report (not peer-reviewed) · **Date:** 2026-07

> ⚠️ 예비 연구/technical report. 임상 유효성·인과 주장 아님. 공개 PDB 구조 사용. 블라인드 자소서 본문엔 실명·ORCID 미기재(이 공개본에만).

## Abstract
폐암 EGFR 표적치료(오시머티닙)에서 내성은 표적 단백질의 특정 잔기 돌연변이로 생긴다. 본 연구는 **같은 약물(osimertinib, PDB 리간드 YY3)이 결합한** EGFR 키나아제 구조 4개(WT 6JXT, L858R 6JWL, T790M 6JX0, L858R/T790M/C797S 6LUD)를 사용해 리간드 유도적합 교란을 통제하고, 결합 포켓(리간드 8 Å 이내 42잔기)의 구조·위상 변화를 정량화했다. 지속 호몰로지(persistent homology)로 계산한 **포켓 전체 위상 거리는 좌표 잡음 기준선(Wasserstein p95 ≈ 23.7)보다 작아**(모든 변이 12–14) 점돌연변이를 구별하지 못했다 — 전역 위상이 열적·결정학적 잡음에 지배됨을 뜻한다. 그러나 잔기별 국소 변형(중첩 후 Hausdorff)은 임상 돌연변이 자리에 뚜렷이 국소화되었고, 그 크기가 내성 유형을 **등급화**했다: **L858R 9.4 Å ≫ T790M 2.4 Å ≫ C797S 1.08 Å**(구조적 배경 중앙값 1.32 Å). 특히 **C797S는 배경보다도 작아 구조적으로 '침묵'** 했는데, 이는 Cys→Ser가 거의 등입체(isosteric)여서 오시머티닙 내성이 **모양이 아니라 화학(공유결합 앵커 소실)** 에서 옴을 반영한다. 요약하면, 구조 기반 기술은 입체적 내성(T790M·L858R)은 짚어내지만 화학적 내성(C797S)에는 원리적으로 눈이 멀며, 이는 구조·위상 특징과 화학·공유결합 특징을 함께 써야 함을 보여준다.

## 1. Background
EGFR 활성화 돌연변이(L858R, exon19del)를 가진 비소세포폐암은 표적치료에 잘 반응하지만, 대부분 내성으로 재발한다. 내성의 대표 경로는 위상적으로 다르다: (i) **T790M** 게이트키퍼는 곁사슬을 키워 1·2세대 저해제의 결합을 방해하고(3세대 오시머티닙으로 극복), (ii) **C797S** 는 오시머티닙이 공유결합을 맺는 시스테인을 세린으로 바꿔 3세대 약까지 무력화한다. 지원자의 선행 연구(매듭이론·지속 호몰로지)를 표적치료 내성이라는 임상 문제로 옮겨, **"돌연변이가 결합 포켓의 구조·위상에 남기는 발자국"** 을 정량화할 수 있는지 예비 검증한다.

## 2. Data & Method
- **구조**: RCSB PDB, 오시머티닙(YY3) 결합 EGFR 키나아제 도메인 — 6JXT(WT), 6JWL(L858R), 6JX0(T790M), 6LUD(L858R/T790M/C797S). 해상도 2.05–2.55 Å. 같은 리간드로 유도적합 통제. (보조: gefitinib 결합 2ITY/4I22.)
- **포켓 정의**: WT 기준 리간드 8 Å 이내 42잔기(P-loop 716–727, 촉매 K745, 게이트키퍼 T790, 힌지 Q791–M793, 공유앵커 C797, DFG D855, 활성화 L858 포함). 동일 잔기집합을 전 구조에 적용.
- **위상**: 포켓 중원자 점구름의 Vietoris–Rips 지속 호몰로지(H0/H1/H2, `ripser`). 변이 간 Wasserstein 거리(`persim`). 좌표에 σ=0.4 Å 가우시안 잡음을 준 자기-자기 분포로 잡음 기준선 추정.
- **국소 변형**: 포켓 Cα로 중첩(Biopython Superimposer) 후, 잔기별 WT↔변이 중원자 **대칭 Hausdorff 거리**로 국소 변형 정량화. 비돌연변이 잔기의 중앙값을 구조적 배경으로.

## 3. Results
**(A) 전역 위상은 잡음에 지배됨.** 포켓 전체 Wasserstein(H1+H2) 거리는 WT 대비 L858R 13.8, T790M 12.3, 삼중변이 13.2로 모두 잡음 기준선 p95(23.7)보다 작았다. 같은 단백질의 점돌연변이는 backbone을 거의 바꾸지 않으므로(포켓 Cα RMSD ≈ 1.25 Å) 전역 위상만으로는 유전형을 가르지 못한다. (그림: `fig_global_noise.png`, `fig_persistence_H1.png`)

**(B) 변형은 임상 돌연변이 자리에 국소화.** 중첩 후 잔기별 Hausdorff는 돌연변이 잔기에서 최대였다(그림 `fig_localization.png`).

**(C) 발자국 크기가 내성 유형을 등급화(머니 피규어 `fig_site_footprint.png`):**

| 돌연변이 | 자리 | 구조 발자국(Hausdorff) | 유형 | 감지 |
|---|---|---|---|---|
| L858R | 858 (활성화 루프) | **9.42 Å** | 입체 | ✅ |
| T790M | 790 (게이트키퍼) | **2.39 Å** | 입체 | ✅ |
| C797S | 797 (공유앵커) | **1.08 Å** | 화학 | ⚠️ 배경(1.32 Å) 이하 — 침묵 |

## 4. Interpretation
구조/위상 기반 기술은 **입체적 내성(T790M·L858R)** 을 잔기 수준에서 국소화·등급화한다. 반면 **C797S 는 등입체 치환**이라 구조적 발자국이 배경 이하로 사라진다 — 오시머티닙 내성이 포켓 **모양이 아니라 공유결합 화학**의 소실에서 오기 때문이다. 즉 "구조가 침묵한다"는 관찰 자체가 내성 유형을 진단하는 신호가 된다. 이는 서로 다른 물리적 기전을 가진 내성(결합자리 화학의 변화 vs 입체·구조의 변화)을 구별하려면 구조·위상 특징에 **화학·공유결합 특징**을 결합해야 함을 시사한다.

## 5. Limitations
유전형당 결정구조 n=1(예비), 단일 스냅샷(동역학 아님), in-crystallo 상태, 인과 아님. Hausdorff·지속 호몰로지는 좌표 잡음에 민감하므로 절대값보다 상대·국소 비교로 해석. 임상 적용 주장 없음.

## 6. Conclusion
같은 약물이 결합한 통제 구조에서, 잔기 국소 변형이 EGFR 표적치료 내성 돌연변이를 유형별로 짚어내되 화학적 내성(C797S)에는 원리적으로 눈이 먼다는 것을 예비적으로 보였다. 후속: 앙상블/MD 기반 위상, 공유결합·전자적 특징 결합, 다중 표적(BCR-ABL T315I, KRAS G12C)으로 확장.

## Artifacts
Code + figures: [GitHub](https://github.com/lshpy/egfr-resistance-topology) · This note: [Zenodo DOI 발급 예정] · ORCID: 0009-0006-1926-653X

## References
[1] Xia & Wei (2014), *Persistent homology analysis of protein structure, flexibility and folding*, Int. J. Numer. Methods Biomed. Eng.
[2] Cang & Wei (2017), *TopologyNet*, PLoS Comput. Biol.
[3] RCSB PDB: 6JXT, 6JWL, 6JX0, 6LUD, 2ITY, 4I22.
[4] Osimertinib covalent binding to EGFR C797; C797S resistance — see NSCLC targeted-therapy literature.
