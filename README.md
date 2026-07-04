# 🧬 EGFR Resistance Topology — 표적치료 내성의 구조·위상 발자국

폐암 EGFR 표적치료(오시머티닙)의 **내성 돌연변이**가 약물 결합 포켓에 남기는 구조·위상 발자국을 정량화하는 예비연구. 지속 호몰로지(persistent homology)와 잔기별 구조 변형으로 내성 **유형**을 구별한다.

> ⚠️ 예비 연구/technical report — 임상 도구 아님. 공개 PDB 구조 기반.

## 핵심 발견
같은 약물(osimertinib)이 결합한 EGFR 구조들만 비교해 리간드 교란을 통제한 뒤, 내성 돌연변이의 **구조적 발자국**을 측정:

| 돌연변이 | 자리 | 발자국(Hausdorff) | 내성 유형 |
|---|---|---|---|
| **L858R** | 활성화 루프 | 9.4 Å | 입체 — 감지됨 |
| **T790M** | 게이트키퍼 | 2.4 Å | 입체 — 감지됨 |
| **C797S** | 공유앵커 | 1.08 Å | 화학 — **구조적으로 침묵** |

**C797S는 배경(1.3 Å)보다도 작다.** Cys→Ser가 거의 등입체라 오시머티닙 내성이 *모양이 아니라 공유결합 화학의 소실*에서 오기 때문. → 구조/위상 기술은 입체적 내성은 짚지만 화학적 내성엔 원리적으로 눈이 먼다.

![money figure](results/fig_site_footprint.png)

## 방법
- **구조**: PDB 6JXT(WT)·6JWL(L858R)·6JX0(T790M)·6LUD(삼중변이), 모두 osimertinib 결합, 2.05–2.55 Å.
- **포켓**: 리간드 8 Å 이내 42잔기(게이트키퍼 790·공유앵커 797·활성화 858 포함).
- **위상**: 포켓 중원자 점구름의 Vietoris–Rips 지속 호몰로지(ripser) + Wasserstein 거리(persim), 좌표잡음 기준선.
- **국소 변형**: Cα 중첩 후 잔기별 대칭 Hausdorff 거리.

## 구조
```
src/data.py       PDB 다운로드·포켓 정의·점구름 추출
src/topology.py   지속 호몰로지·Wasserstein·잡음 기준선·잔기 기여도
src/geometry.py   구조 중첩·잔기별 Hausdorff·국소 위상
src/analyze.py    전체 파이프라인 → results/ (숫자 + 그림)
app.py            Gradio 데모(돌연변이 선택 → 발자국·위상·해석)
technical_note/   기술노트(Zenodo 업로드용)
scripts/probe_pdb.py  후보 PDB 구조 판별 유틸
```

## 실행
```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
python -m src.analyze     # 분석 재현 (results/ 생성)
python app.py             # 데모 실행
```

## 데이터 출처
RCSB PDB (https://www.rcsb.org) — 공개 구조. 각 PDB 원저작권은 해당 기탁자/문헌에 있음.

## 인용 / DOI
릴리즈는 Zenodo에 아카이빙되어 DOI가 발급됩니다. (발급 후 배지 추가) · ORCID 0009-0006-1926-653X

## 라이선스
코드 MIT. 구조 데이터는 RCSB PDB 정책을 따름.
