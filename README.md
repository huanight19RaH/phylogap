# PhyloGap: Evidence-Aware Genome Gap Reconstruction

[![Status: Scoping & Novelty Validation](https://img.shields.io/badge/Status-Scoping%20%26%20Novelty%20Validation-blue.svg)](file:///d:/THStudy/UniversityStudy/Research/phylogap/RESEARCH_PLAN.md)
[![Scope: Computational Research](https://img.shields.io/badge/Scope-In%20Silico%20Dry%20Lab-success.svg)](file:///d:/THStudy/UniversityStudy/Research/phylogap/RESEARCH_PLAN.md)
[![Domain: Avian Genomes](https://img.shields.io/badge/Domain-Avian%20Genomes%20(B10K)-orange.svg)](file:///d:/THStudy/UniversityStudy/Research/phylogap/RESEARCH_PLAN.md)

**PhyloGap** is a computational research project aimed at benchmark-driven and evidence-aware reconstruction of genome assembly gaps under constrained reference and sequencing read scenarios, with a primary application domain in avian genomes (*Aves*).

---

## 📌 Project Overview

* **Working Description:** Genome Gap Reconstruction through Cross-Species Evidence and Reliable Prediction
* **Provisional Research Title:** *Cross-Species Benchmarking of Genome Assembly Gap-Filling Methods under Reference and Read-Data Constraints*
* **Core Philosophy:**
  $$\text{Reproduce Baselines} \longrightarrow \text{Identify Empirical Failure Modes} \longrightarrow \text{Targeted Algorithmic Innovation}$$
* **Research Plan Master File:** [RESEARCH_PLAN.md](file:///d:/THStudy/UniversityStudy/Research/phylogap/RESEARCH_PLAN.md)
* **Detailed Execution Plan (SOP):** [docs/DETAILED_EXECUTION_PLAN.md](file:///d:/THStudy/UniversityStudy/Research/phylogap/docs/DETAILED_EXECUTION_PLAN.md)

---

## 🔬 Research Questions (RQ1–RQ5)

* **RQ1 (Baseline Fidelity):** How do existing read-based (DENTIST, Sealer) and learning-based (GapPredict, DLGapCloser, DL-GapFilling) gap-filling methods perform on standardized benchmarks with verifiable ground truth? *(Mandatory Phase 1)*
* **RQ2 (Failure Mode Taxonomy):** What structural gap categories are most prone to failure (gap length, repetitive elements, low sequence conservation, or lack of read spanning)? *(Mandatory Phase 1)*
* **RQ3 (Phylogenetic Distance Impact):** How does evolutionary divergence between target and reference species correlate with nucleotide-level reconstruction accuracy? *(Mandatory Phase 1)*
* **RQ4 (Multi-Source Evidence Fusion):** Can combining orthogonal evidence streams (reads, homology, DNA language priors) outperform individual modalities? *(Conditional Phase 2)*
* **RQ5 (Selective Prediction & Safety):** Can confidence scoring and prediction rejection mechanisms suppress catastrophic false closures? *(Conditional Phase 2)*

---

## 📂 Repository Layout

```text
phylogap/
├── README.md                      # Project overview & quickstart (this file)
├── RESEARCH_PLAN.md               # Master research plan v0.1 (detailed blueprint)
├── docs/                          # Detailed research documentation & notes
│   ├── DETAILED_EXECUTION_PLAN.md # Phase-by-phase execution protocol & Done Criteria
│   ├── literature_matrix.md       # SOTA comparison & literature audit (7 dimensions)
│   └── data_manifest.md           # Accessions, species list, quality metrics & gap stats
├── data/                          # Datasets (gitignored large files)
│   └── benchmarks/dentist_pilot/  # Pilot benchmark assemblies & ground truth BED/JSON
├── scripts/                       # Automated data download, evaluation & audit scripts
│   ├── download_dentist_metadata.py
│   ├── download_pilot_assemblies.py
│   ├── audit_gap_distribution.py
│   └── analyze_gap_sequences.py   # Biological EDA & GC/repeat analysis
└── src/                           # Custom algorithm codebase (conditional on Gates 4-5)
```

---

## 🧭 Roadmap & Decision Gates

The project follows a rigorous 14-week schedule with hard **Go/No-Go Decision Gates**:
1. **Gate 1 (Week 2):** Literature differentiation vs. recent SOTA (DL-GapFilling 2026, GapSuite 2026).
2. **Gate 2 (Week 4):** Verifiable benchmark with reliable ground truth and reproducibility scripts.
3. **Gate 3 (Week 6):** Successful reproduction of DENTIST and at least one AI baseline.
4. **Gate 4 (Week 8):** Concrete statistical identification of a critical failure mode. *If absent, do NOT build a new DL model.*
5. **Gate 5 (Week 12):** Validated Pareto improvement (Accuracy vs. False-Closure Rate) on an unseen held-out species test.

For full technical specifications, metric formulations, and risk mitigations, consult the [RESEARCH_PLAN.md](file:///d:/THStudy/UniversityStudy/Research/phylogap/RESEARCH_PLAN.md).
