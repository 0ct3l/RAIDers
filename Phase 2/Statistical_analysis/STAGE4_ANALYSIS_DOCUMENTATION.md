# Stage 4: Comprehensive Pathway Analysis - Documentation

## Overview

This script performs pathway co-occurrence and correlation analysis using **both** count-based and burden-based scoring methods, enabling direct comparison of their predictive power.

---

## Input Requirements

**Required file**: `patients_with_pathways_weighted.csv` from Phase 2 weighted script

**Required columns**:
- `pathway_{name}` — Binary (0/1)
- `pathway_{name}_count` — Gene count
- `pathway_{name}_weighted` — Sum of gene evidence weights
- `pathway_{name}_interaction` — Gene interaction multiplier
- `pathway_{name}_burden` — Normalized burden score
- `composite_score` — Patient's total weighted score
- `severity_category` — Clinical severity classification
- `severity_score` — Numeric severity (if available)
- `all_genes` — List of patient's genes

---

## Output Files

### 4.1: Pathway Prevalence
**File**: `4.1_pathway_prevalence_by_severity.csv`

| Column | Description |
|--------|-------------|
| severity | Mild/Moderate/Severe |
| pathway | Pathway name |
| count | Number of patients with pathway |
| pct | Percentage of severity category |
| mean_burden | Mean burden score for affected patients |

**Use case**: Identify which pathways are enriched in severe patients.

---

### 4.2A: Co-occurrence
**Files**: 
- `4.2A_cooccurrence_matrices.xlsx` — Count matrices by severity
- `4.2A_cooccurrence_statistics.csv` — Statistical tests

| Column | Description |
|--------|-------------|
| pathway_1, pathway_2 | Pathway pair |
| cooccurrence_count | Patients with both |
| expected_count | Expected under independence |
| obs_exp_ratio | Observed/Expected (>1 = enriched) |
| jaccard | Jaccard similarity coefficient |
| odds_ratio | Effect size |
| fisher_neglog10p | -log10(p-value) from Fisher's exact |

**Use case**: Identify pathway pairs that co-occur more than expected.

---

### 4.2B: Correlations
**Files**: 
- `4.2B_correlation_matrices_count.xlsx` — Count-based correlations
- `4.2B_correlation_matrices_burden.xlsx` — Burden-based correlations

**Method**: Spearman rank correlation between pathway scores.

**Use case**: Compare whether burden scores show different correlation patterns than simple counts.

---

### 4.3: Frequency vs Intensity
**File**: `4.3_frequency_vs_intensity.csv`

Compares multiple scoring methods side-by-side:

| Column | Description |
|--------|-------------|
| frequency_pct | % of patients with pathway |
| mean_count | Mean gene count (affected only) |
| mean_weighted | Mean weighted sum |
| mean_interaction | Mean interaction multiplier |
| mean_burden | Mean normalized burden |

**Use case**: Identify pathways where intensity matters more than frequency.

---

### 4.4: Network Visualization
**Files**:
- `4.4_network_data.json` — Full network for visualization tools
- `4.4_network_nodes.csv` — Node attributes
- `4.4_network_edges.csv` — Edge attributes

**Edge types (4 total)**:

| Edge Weight | Description | Interpretation |
|-------------|-------------|----------------|
| `cooccurrence_count` | Patients with both pathways | Raw frequency |
| `jaccard` | Jaccard similarity (0-1) | Normalized overlap |
| `burden_correlation` | Spearman r of burden scores | Score relationship |
| `mean_interaction_multiplier` | Mean gene interaction when co-occurring | Synergy potential |

**Use case**: Load into Cytoscape, D3.js, or similar for visualization.

<img width="729" height="866" alt="Screenshot 2026-01-29 at 10 43 28 am" src="https://github.com/user-attachments/assets/ac0485dc-6a5d-40dd-bda7-0ff763c3b02f" />


---

### 4.5: Cross-Severity Comparison
**File**: `4.5_cross_severity_comparison.csv`

| Column | Description |
|--------|-------------|
| chi2_*, cramers_v | Chi-square test on binary presence |
| kw_count_* | Kruskal-Wallis on gene counts |
| kw_burden_* | Kruskal-Wallis on burden scores |

**Effect size interpretation** (Cramér's V, epsilon²):
- < 0.1: Negligible
- 0.1 - 0.3: Small
- 0.3 - 0.5: Medium
- > 0.5: Large

**Use case**: Identify which pathways differentiate severity categories.

---

### 4.6: Validation (NEW)
**Files**:
- `4.6_validation_count_vs_burden.csv` — Detailed test results
- `4.6_validation_summary.txt` — Summary and conclusion

**Tests performed**:

| Test | What it measures | Better = Higher |
|------|------------------|-----------------|
| Spearman r with severity | Correlation with clinical score | ✓ |
| AUC (Mild vs Severe) | Discrimination ability | ✓ |
| Cohen's d | Effect size between groups | ✓ |

**Use case**: Determine if weighted scoring improves predictions over simple counts.

---

### 4.7: Interaction Effects (NEW)
**Files**:
- `4.7_interaction_effects_population.csv` — Per-interaction statistics
- `4.7_interaction_effects_by_severity.csv` — Interaction rates by severity

**Population-level columns**:

| Column | Description |
|--------|-------------|
| interaction_pair | Gene pair (e.g., OPTN-TBK1) |
| n_patients | How many patients have this pair |
| multiplier | Synergy (>1) or redundancy (<1) |
| type | synergy/redundancy/epistatic |
| mean_severity | Mean severity for patients with interaction |
| severity_difference | Compared to patients without interaction |

**Severity-level columns**:

| Column | Description |
|--------|-------------|
| pct_with_synergy | % of patients with synergistic interactions |
| mean_interaction_multiplier | Mean multiplier by severity |

**Use case**: Identify gene pairs that contribute to worse outcomes.

---

## Patient-Level Interaction Annotations

The script adds these columns to the carrier data:

| Column | Description |
|--------|-------------|
| n_interactions | Count of known interactions |
| interaction_pairs | Which pairs (e.g., "OPTN-TBK1; SOD1-SIGMAR1") |
| combined_interaction_multiplier | Product of all multipliers |
| n_synergies | Count of synergistic interactions |
| n_redundancies | Count of redundant interactions |
| net_interaction_effect | synergy/redundancy/neutral/none |

---

## Interpretation Guide

### When Burden > Count (in validation):
- The evidence weights and interaction multipliers capture biological reality
- High-penetrance genes appropriately dominate the score
- Synergistic interactions amplify risk

### When Count ≈ Burden:
- For this dataset, simple counts may be sufficient
- Or: the weight estimates need refinement
- Or: the patient population lacks variation in gene severity

### When Count > Burden:
- Weighted scoring may be overcorrecting
- Weight estimates may not match this population
- Simple approach may be more robust

---

## Network Visualization Tips

**For Cytoscape**:
1. Import `4.4_network_nodes.csv` as node table
2. Import `4.4_network_edges.csv` as edge table
3. Map `size` → Node Size
4. Map `cooccurrence_count` OR `burden_correlation` → Edge Width
5. Map `mean_interaction_multiplier` → Edge Color (red = synergy)

**Edge weight choice depends on question**:
- "Which pathways co-occur?" → Use `cooccurrence_count`
- "Which pathways are biologically related?" → Use `burden_correlation`
- "Which combinations have synergistic risk?" → Use `mean_interaction_multiplier`

---

## Statistical Notes

1. **P-values reported as -log10(p)**: 
   - Value of 2 = p=0.01
   - Value of 10 = p=1e-10
   - Value of 300 = underflow (p < 1e-300)

2. **Focus on effect sizes**: With large N, all p-values will be significant. Effect sizes (Cramér's V, Cohen's d, AUC) indicate clinical meaningfulness.

3. **Multiple testing**: No correction applied. Use effect sizes for prioritization.

---

## Citation

If using this analysis, cite the methodology documentation and note that:
- Gene weights are estimates from ALS literature (ClinGen, ALSoD)
- Interaction multipliers are based on published mechanistic studies
- Validation against clinical severity is essential for any application
