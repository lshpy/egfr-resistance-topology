# Structural Footprints of Targeted-Therapy Resistance Mutations in the EGFR Kinase Pocket: A Topology-Aware Preliminary Study

**Author:** Seunghyun Lee — ORCID [0009-0006-1926-653X](https://orcid.org/0009-0006-1926-653X)
**Type:** Preliminary technical report (not peer-reviewed) · **Date:** 2026-07

> ⚠️ Preliminary study / technical report. No claims of clinical validity or causality. Uses public PDB structures. The author's name and ORCID appear only in this public version.

## Abstract
In EGFR targeted therapy for lung cancer (osimertinib), resistance arises from point mutations at specific residues of the target protein. This study uses four EGFR kinase structures **bound to the same drug (osimertinib, PDB ligand YY3)** (WT 6JXT, L858R 6JWL, T790M 6JX0, L858R/T790M/C797S 6LUD) to control for ligand induced-fit confounding, and quantifies structural and topological changes of the binding pocket (42 residues within 8 Å of the ligand). The **whole-pocket topological distance** computed with persistent homology was **smaller than the coordinate-noise baseline (Wasserstein p95 ≈ 23.7)** (12–14 for every variant) and could not distinguish the point mutations, which means global topology is dominated by thermal and crystallographic noise. Per-residue local deformation (Hausdorff after superposition), however, was clearly localized at the clinical mutation sites, and its magnitude **graded** the resistance types: **L858R 9.4 Å ≫ T790M 2.4 Å ≫ C797S 1.08 Å** (structural background median 1.32 Å). Notably, **C797S was smaller than the background, i.e. structurally "silent"**, reflecting that Cys→Ser is nearly isosteric and osimertinib resistance therefore comes **from chemistry (loss of the covalent anchor), not shape**. In short, structure-based descriptors pick up steric resistance (T790M, L858R) but are in principle blind to chemical resistance (C797S), which shows that structural/topological features must be combined with chemical/covalent features.

## 1. Background
Non-small-cell lung cancers carrying EGFR activating mutations (L858R, exon19del) respond well to targeted therapy, but most relapse through resistance. The main resistance routes differ in kind: (i) the **T790M** gatekeeper enlarges the side chain and blocks binding of 1st/2nd-generation inhibitors (overcome by the 3rd-generation osimertinib), and (ii) **C797S** replaces the cysteine with which osimertinib forms a covalent bond by a serine, disabling even the 3rd-generation drug. This work carries the author's prior work (knot theory, persistent homology) over to the clinical problem of targeted-therapy resistance and tests, as a preliminary study, whether **"the footprint a mutation leaves in the structure and topology of the binding pocket"** can be quantified.

## 2. Data & Method
- **Structures**: RCSB PDB, osimertinib (YY3)-bound EGFR kinase domain: 6JXT (WT), 6JWL (L858R), 6JX0 (T790M), 6LUD (L858R/T790M/C797S). Resolution 2.05–2.55 Å. Induced fit is controlled by using the same ligand. (Auxiliary: gefitinib-bound 2ITY/4I22.)
- **Pocket definition**: 42 residues within 8 Å of the ligand in the WT reference (including P-loop 716–727, catalytic K745, gatekeeper T790, hinge Q791–M793, covalent anchor C797, DFG D855, activation-loop L858). The same residue set is applied to all structures.
- **Topology**: Vietoris–Rips persistent homology (H0/H1/H2, `ripser`) of the pocket heavy-atom point cloud. Wasserstein distances between variants (`persim`). The noise baseline is estimated from the self-vs-self distribution after adding σ = 0.4 Å Gaussian noise to the coordinates.
- **Local deformation**: after superposition on pocket Cα atoms (Biopython Superimposer), local deformation is quantified by the per-residue WT↔mutant heavy-atom **symmetric Hausdorff distance**. The median over non-mutated residues serves as the structural background.

## 3. Results
**(A) Global topology is noise-dominated.** Whole-pocket Wasserstein (H1+H2) distances from WT were 13.8 for L858R, 12.3 for T790M and 13.2 for the triple mutant, all below the noise-baseline p95 (23.7). Point mutations in the same protein barely change the backbone (pocket Cα RMSD ≈ 1.25 Å), so global topology alone cannot separate genotypes. (Figures: `fig_global_noise.png`, `fig_persistence_H1.png`)

**(B) Deformation localizes at the clinical mutation sites.** Per-residue Hausdorff after superposition was largest at the mutated residues (figure `fig_localization.png`).

**(C) Footprint size grades the resistance type (money figure `fig_site_footprint.png`):**

| Mutation | Site | Structural footprint (Hausdorff) | Type | Detected |
|---|---|---|---|---|
| L858R | 858 (activation loop) | **9.42 Å** | steric | ✅ |
| T790M | 790 (gatekeeper) | **2.39 Å** | steric | ✅ |
| C797S | 797 (covalent anchor) | **1.08 Å** | chemical | ⚠️ below background (1.32 Å): silent |

## 4. Interpretation
Structure/topology-based descriptors localize and grade **steric resistance (T790M, L858R)** at the residue level. **C797S, by contrast, is an isosteric substitution**, so its structural footprint vanishes below the background, because osimertinib resistance comes from the loss of **covalent chemistry, not pocket shape**. The observation that "the structure is silent" is thus itself a diagnostic signal of the resistance type. This suggests that distinguishing resistance with different physical mechanisms (a change in binding-site chemistry vs a steric/structural change) requires combining structural/topological features with **chemical/covalent features**.

## 5. Limitations
n = 1 crystal structure per genotype (preliminary), a single snapshot (not dynamics), in-crystallo state, not causal. Hausdorff distances and persistent homology are sensitive to coordinate noise, so they should be read as relative, local comparisons rather than absolute values. No claim of clinical applicability.

## 6. Conclusion
In controlled structures bound to the same drug, we show preliminarily that local residue deformation identifies EGFR targeted-therapy resistance mutations by type, while being in principle blind to chemical resistance (C797S). Next steps: ensemble/MD-based topology, combining covalent and electronic features, and extension to other targets (BCR-ABL T315I, KRAS G12C).

## Artifacts
Code + figures: [GitHub](https://github.com/lshpy/egfr-resistance-topology) · This note: [Zenodo DOI 10.5281/zenodo.21193665](https://doi.org/10.5281/zenodo.21193665) · ORCID: [0009-0006-1926-653X](https://orcid.org/0009-0006-1926-653X)

## References
[1] Xia & Wei (2014), *Persistent homology analysis of protein structure, flexibility and folding*, Int. J. Numer. Methods Biomed. Eng.
[2] Cang & Wei (2017), *TopologyNet*, PLoS Comput. Biol.
[3] RCSB PDB: 6JXT, 6JWL, 6JX0, 6LUD, 2ITY, 4I22.
[4] Osimertinib covalent binding to EGFR C797; C797S resistance — see NSCLC targeted-therapy literature.
