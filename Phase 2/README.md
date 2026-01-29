

# RAIDers Pipeline: Phase 1 → Phase 2 Transition

**Internal Reference** **Last Updated:** 16 January 2026

**Status:** 🟢 Phase 1 Complete | 🟡 Phase 2 In Development

---

## 🏗️ Phase 1 Recap: Infrastructure & Severity Modeling

Phase 1 focused on building a synthetic biobank infrastructure from the ground up, moving from raw ClinVar data to a federated clustering validation.

### Key Achievements

* **Data Standardization:** Normalized coordinates and extracted gene symbols from raw ClinVar variant data.
* **Population Simulation:** Generated **15,000 synthetic patients** (3,000 per population: AFR, AMR, EAS, EUR, SAS) using Hardy-Weinberg equilibrium and gnomAD-style frequency constraints.
* **Contextual Penetrance Model:** Developed a novel scoring system where ancestral modifiers modulate disease severity.
* **Federated Infrastructure:** Created privacy-partitioned client datasets and successfully ran federated k-means clustering.

### The Interaction Score Formula

To simulate biological complexity, we implemented the following formula:

Where:

* I: Base Impact
* M: Ancestral Modifier (Protective/Neutral/Aggravating)
* C: Consequence Multiplier (e.g., LoF vs. Missense)
* S: Stochastic Noise

---

## 🔍 The Problem: Severity vs. Biology

While Phase 1 produced statistically significant clusters, they were **severity-stratified** rather than **pathway-based**.

### Phase 1 Clustering Results

| Cluster | Severity | Progression | Top Genes |
| --- | --- | --- | --- |
| **C0** | Mild (5.0) | 74% Slow | SETX, ALS2 |
| **C4** | Moderate (6.3) | 91% Moderate | ALS2, OPTN |
| **C1** | Mod-Severe (7.3) | 93% Moderate | SOD1 |
| **C2** | Severe (9.3) | 96.5% Fast | TBK1, FUS |

**The Bottleneck:** In Phase 1, **FUS** and **TBK1** cluster together because they are both "severe." However, they require radically different treatments:

* **FUS** (RNA Metabolism)  Requires ASO therapy.
* **TBK1** (Autophagy)  Requires autophagy inducers.

---

## 🎯 Phase 2 Goal: Pathway-Driven Subtyping

The objective of Phase 2 is to make **biological pathways** the primary clustering signal. We are moving from "How sick is the patient?" to "Why is the patient sick?"

### The Strategy

1. **Explicit Pathway Mapping:** Map every gene to its dysfunctional pathway (RNA Metabolism, Autophagy, Proteostasis, etc.).
2. **One-Hot Encoding:** Convert pathway labels into binary features for the clustering algorithm.
3. **Severity Decoupling:** Remove severity as a clustering feature and treat it as an **outcome variable** for post-hoc analysis.

### Expected Outcome (Target Clustering)

| Cluster | Dominant Pathway | Severity Range | Therapeutic Implication |
| --- | --- | --- | --- |
| **C0** | RNA Metabolism | 4.0 – 9.0 | ASOs, splicing modulators |
| **C1** | Autophagy | 5.0 – 8.0 | Rapamycin, trehalose |
| **C2** | Proteostasis | 6.0 – 9.0 | HSP inducers, CuATSM |
| **C3** | Vesicle Trafficking | 4.0 – 7.0 | Emerging targets |

---

## 📊 Summary of Shift

| Feature | Phase 1 (Legacy) | Phase 2 (Current) |
| --- | --- | --- |
| **Primary Signal** | Severity | Pathway |
| **Cluster Logic** | "How severe is the patient?" | "Which pathway is broken?" |
| **Utility** | Prognosis / Triage | Treatment Selection |
| **Severity Role** | Input Feature | Outcome Variable |

---

## 🛠️ Implementation Roadmap

* [ ] **Literature Validation:** Verify gene-to-pathway mappings against primary sources.
* [ ] **Feature Engineering:** Implement `GENE_TO_PATHWAY` dictionary and one-hot encoding functions.
* [ ] **Pipeline Update:** Modify the federated k-means script to utilize pathway columns.
* [ ] **Validation:** Confirm clusters show within-pathway severity variation.
* [ ] **Bonus:** Integrate a privacy hashing layer for client-side data.
