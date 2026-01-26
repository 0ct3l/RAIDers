#!/usr/bin/env python3
"""
Stage 4: Pathway Co-occurrence and Correlation Analysis (COMPREHENSIVE VERSION)
================================================================================
Analyzes pathway patterns using BOTH count-based and burden-based scoring.

REQUIREMENTS:
  - Input from Phase 2 WEIGHTED script (patients_with_pathways_weighted.csv)
  - Contains both _count and _burden columns

ANALYSES:
  4.1: Pathway prevalence by severity (binary)
  4.2A: Co-occurrence matrices (binary)
  4.2B: Correlation matrices (count-based AND burden-based)
  4.3: Frequency vs intensity comparison (all scoring methods)
  4.4: Network visualization (4 edge types: co-occurrence, Jaccard, correlation, interaction)
  4.5: Cross-severity comparison (with effect sizes)
  4.6: VALIDATION - Count vs Burden scoring comparison
  4.7: INTERACTION EFFECT ANALYSIS - Gene pair impacts

OUTPUT STRUCTURE:
  /4.1_pathway_prevalence_by_severity.csv
  /4.2A_cooccurrence_matrices.xlsx
  /4.2A_cooccurrence_statistics.csv
  /4.2B_correlation_matrices_count.xlsx
  /4.2B_correlation_matrices_burden.xlsx
  /4.3_frequency_vs_intensity.csv
  /4.4_network_data.json
  /4.4_network_nodes.csv
  /4.4_network_edges.csv
  /4.5_cross_severity_comparison.csv
  /4.6_validation_count_vs_burden.csv
  /4.6_validation_summary.txt
  /4.7_interaction_effects_population.csv
  /4.7_interaction_effects_by_severity.csv
"""

import pandas as pd
import numpy as np
from scipy import stats
from sklearn.metrics import roc_auc_score, classification_report
from sklearn.preprocessing import LabelEncoder
from itertools import combinations
import json
import os
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning)

# =============================================================================
# CONFIGURATION
# =============================================================================

INPUT_FILE = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data/phase2_output_weighted/patients_with_pathways_weighted.csv")
OUTPUT_DIR = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/stage4_output")

PATHWAYS = [
    'Proteostasis',
    'RNA_Metabolism', 
    'Cytoskeletal_Axonal_Transport',
    'Mitochondrial',
    'Excitotoxicity',
    'Vesicle_Trafficking',
    'DNA_Damage',
]

SEVERITY_ORDER = ['Unaffected', 'Mild', 'Moderate', 'Severe']
SEVERITY_NUMERIC = {'Unaffected': 0, 'Mild': 1, 'Moderate': 2, 'Severe': 3}

# =============================================================================
# GENE-GENE INTERACTIONS (from Phase 2 - duplicated here for analysis)
# =============================================================================

GENE_INTERACTIONS = {
    ('TARDBP', 'FUS'): {'multiplier': 1.3, 'type': 'synergy', 'mechanism': 'Co-aggregation, overlapping RNA targets'},
    ('OPTN', 'TBK1'): {'multiplier': 1.5, 'type': 'synergy', 'mechanism': 'TBK1 activates OPTN; dual loss = autophagy collapse'},
    ('SOD1', 'SIGMAR1'): {'multiplier': 1.4, 'type': 'synergy', 'mechanism': 'Both impair mitochondrial Ca2+ handling'},
    ('C9ORF72', 'TARDBP'): {'multiplier': 1.3, 'type': 'synergy', 'mechanism': 'DPR proteins enhance TDP-43 mislocalization'},
    ('C9ORF72', 'FUS'): {'multiplier': 1.2, 'type': 'synergy', 'mechanism': 'DPR proteins affect FUS nuclear transport'},
    ('SOD1', 'CHCHD10'): {'multiplier': 1.3, 'type': 'synergy', 'mechanism': 'Convergent mitochondrial toxicity'},
    ('FUS', 'CHCHD10'): {'multiplier': 1.2, 'type': 'synergy', 'mechanism': 'Both target mitochondrial energy production'},
    ('NEK1', 'TBK1'): {'multiplier': 1.2, 'type': 'synergy', 'mechanism': 'Both kinases regulate stress responses'},
    ('HNRNPA1', 'HNRNPA2B1'): {'multiplier': 0.7, 'type': 'redundancy', 'mechanism': 'Same stress granule mechanism'},
    ('OPTN', 'SQSTM1'): {'multiplier': 0.8, 'type': 'redundancy', 'mechanism': 'Both autophagy receptors'},
    ('TBK1', 'SQSTM1'): {'multiplier': 0.85, 'type': 'epistatic', 'mechanism': 'TBK1 activates p62'},
    ('DCTN1', 'KIF5A'): {'multiplier': 0.9, 'type': 'partial_redundancy', 'mechanism': 'Opposite transport directions'},
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def safe_log10_pvalue(p):
    """Convert p-value to -log10 scale, handling underflow."""
    if p is None or np.isnan(p):
        return np.nan
    if p <= 0:
        return 300.0
    if p >= 1:
        return 0.0
    return min(-np.log10(p), 300.0)


def cramers_v(contingency_table):
    """Calculate Cramér's V effect size for chi-square test."""
    chi2 = stats.chi2_contingency(contingency_table)[0]
    n = contingency_table.sum()
    min_dim = min(contingency_table.shape) - 1
    if min_dim == 0 or n == 0:
        return 0.0
    return np.sqrt(chi2 / (n * min_dim))


def epsilon_squared(h_stat, n):
    """Calculate epsilon squared effect size for Kruskal-Wallis test."""
    if n <= 1:
        return 0.0
    return h_stat / (n - 1)


def interpret_effect_size(value, metric='cramers_v'):
    """Return interpretation of effect size."""
    if pd.isna(value):
        return 'N/A'
    thresholds = {
        'cramers_v': [(0.1, 'negligible'), (0.3, 'small'), (0.5, 'medium'), (1.0, 'large')],
        'epsilon_sq': [(0.01, 'negligible'), (0.06, 'small'), (0.14, 'medium'), (1.0, 'large')],
        'correlation': [(0.1, 'negligible'), (0.3, 'small'), (0.5, 'medium'), (1.0, 'large')],
    }
    for threshold, label in thresholds.get(metric, [(1.0, 'unknown')]):
        if abs(value) < threshold:
            return label
    return 'large'


def parse_genes(all_genes_str):
    """Parse the 'all_genes' column into a list of gene names."""
    if pd.isna(all_genes_str) or all_genes_str == '':
        return []
    if ';' in str(all_genes_str):
        genes = [g.strip() for g in str(all_genes_str).split(';')]
    else:
        genes = [g.strip() for g in str(all_genes_str).split(',')]
    return [g.upper() for g in genes if g]


def get_patient_interactions(genes):
    """
    Identify which known gene-gene interactions a patient has.
    
    Returns:
        list of dicts with interaction details
    """
    if len(genes) < 2:
        return []
    
    interactions_found = []
    gene_set = set(g.upper() for g in genes)
    
    for (g1, g2), info in GENE_INTERACTIONS.items():
        if g1 in gene_set and g2 in gene_set:
            interactions_found.append({
                'gene_pair': f"{g1}-{g2}",
                'multiplier': info['multiplier'],
                'type': info['type'],
                'mechanism': info['mechanism']
            })
    
    return interactions_found


# =============================================================================
# DATA LOADING
# =============================================================================

def load_data():
    """Load and prepare the annotated patient data."""
    print("Loading data...")
    df = pd.read_csv(INPUT_FILE)
    print(f"  Loaded {len(df)} patients")
    print(f"  Columns: {len(df.columns)}")
    
    # Verify we have the weighted columns
    required_cols = ['pathway_Proteostasis_burden', 'composite_score']
    missing = [c for c in required_cols if c not in df.columns]
    if missing:
        print(f"\n⚠️  WARNING: Missing weighted columns: {missing}")
        print("    Make sure to use output from Phase 2 WEIGHTED script!")
        print("    Falling back to count-only analysis...")
    
    # Get carriers only for most analyses
    carriers = df[df['n_pathways_affected'] > 0].copy()
    print(f"  Carriers: {len(carriers)}")
    
    # Add patient-level interaction annotations
    print("  Adding patient-level interaction annotations...")
    carriers = annotate_patient_interactions(carriers)
    
    return df, carriers


def annotate_patient_interactions(df):
    """Add patient-level interaction columns."""
    
    interaction_data = []
    
    for idx, row in df.iterrows():
        genes = parse_genes(row.get('all_genes', ''))
        interactions = get_patient_interactions(genes)
        
        if interactions:
            pairs = [i['gene_pair'] for i in interactions]
            multipliers = [i['multiplier'] for i in interactions]
            types = [i['type'] for i in interactions]
            
            # Calculate combined multiplier
            combined_mult = np.prod(multipliers)
            
            # Count synergies vs redundancies
            n_synergies = sum(1 for t in types if t == 'synergy')
            n_redundancies = sum(1 for t in types if t in ['redundancy', 'partial_redundancy'])
            
            interaction_data.append({
                'n_interactions': len(interactions),
                'interaction_pairs': '; '.join(pairs),
                'interaction_multipliers': '; '.join([f"{m:.2f}" for m in multipliers]),
                'combined_interaction_multiplier': round(combined_mult, 3),
                'n_synergies': n_synergies,
                'n_redundancies': n_redundancies,
                'net_interaction_effect': 'synergy' if combined_mult > 1.05 else ('redundancy' if combined_mult < 0.95 else 'neutral')
            })
        else:
            interaction_data.append({
                'n_interactions': 0,
                'interaction_pairs': '',
                'interaction_multipliers': '',
                'combined_interaction_multiplier': 1.0,
                'n_synergies': 0,
                'n_redundancies': 0,
                'net_interaction_effect': 'none'
            })
    
    interaction_df = pd.DataFrame(interaction_data)
    
    # Merge with original dataframe
    result = pd.concat([df.reset_index(drop=True), interaction_df], axis=1)
    
    # Summary
    n_with_interactions = (result['n_interactions'] > 0).sum()
    print(f"    Patients with known interactions: {n_with_interactions} ({100*n_with_interactions/len(result):.1f}%)")
    
    return result


# =============================================================================
# 4.1: PATHWAY PREVALENCE BY SEVERITY
# =============================================================================

def analyze_pathway_prevalence(df, carriers):
    """Calculate pathway prevalence by severity category."""
    print("\n" + "="*70)
    print("4.1: PATHWAY PREVALENCE BY SEVERITY")
    print("="*70)
    
    results = []
    
    for severity in SEVERITY_ORDER:
        sev_carriers = carriers[carriers['severity_category'] == severity]
        n_carriers = len(sev_carriers)
        
        if n_carriers == 0:
            continue
            
        row = {'severity': severity, 'n_carriers': n_carriers}
        
        for pathway in PATHWAYS:
            col = f'pathway_{pathway}'
            count = sev_carriers[col].sum()
            pct = 100 * count / n_carriers if n_carriers > 0 else 0
            row[f'{pathway}_count'] = count
            row[f'{pathway}_pct'] = round(pct, 1)
            
            # Add mean burden for affected patients
            burden_col = f'pathway_{pathway}_burden'
            if burden_col in sev_carriers.columns:
                affected = sev_carriers[sev_carriers[col] == 1]
                row[f'{pathway}_mean_burden'] = round(affected[burden_col].mean(), 3) if len(affected) > 0 else 0
        
        results.append(row)
        
        print(f"\n{severity} (n={n_carriers}):")
        for pathway in PATHWAYS:
            burden_str = f", mean burden={row.get(f'{pathway}_mean_burden', 'N/A')}" if f'{pathway}_mean_burden' in row else ""
            print(f"  {pathway}: {row[f'{pathway}_count']} ({row[f'{pathway}_pct']:.1f}%){burden_str}")
    
    prevalence_df = pd.DataFrame(results)
    
    # Create pivot version
    pivot_data = []
    for row in results:
        for pathway in PATHWAYS:
            pivot_data.append({
                'severity': row['severity'],
                'pathway': pathway,
                'count': row[f'{pathway}_count'],
                'pct': row[f'{pathway}_pct'],
                'mean_burden': row.get(f'{pathway}_mean_burden', np.nan),
                'n_carriers': row['n_carriers']
            })
    
    pivot_df = pd.DataFrame(pivot_data)
    pivot_df.to_csv(f"{OUTPUT_DIR}/4.1_pathway_prevalence_by_severity.csv", index=False)
    print(f"\nSaved: 4.1_pathway_prevalence_by_severity.csv")
    
    return prevalence_df, pivot_df


# =============================================================================
# 4.2A: CO-OCCURRENCE MATRICES (BINARY)
# =============================================================================

def analyze_cooccurrence(carriers):
    """Calculate pathway co-occurrence matrices with statistics."""
    print("\n" + "="*70)
    print("4.2A: PATHWAY CO-OCCURRENCE MATRICES")
    print("="*70)
    
    all_results = {}
    
    # Overall co-occurrence
    print("\n--- Overall Co-occurrence ---")
    cooc_matrix = calculate_cooccurrence_matrix(carriers)
    all_results['Overall'] = cooc_matrix
    print(cooc_matrix.to_string())
    
    # Calculate statistics
    print("\n--- Co-occurrence Statistics ---")
    cooc_stats = calculate_cooccurrence_statistics(carriers)
    all_results['Overall_Stats'] = cooc_stats
    
    # By severity
    for severity in ['Mild', 'Moderate', 'Severe']:
        sev_carriers = carriers[carriers['severity_category'] == severity]
        if len(sev_carriers) < 10:
            continue
        print(f"\n--- {severity} Severity (n={len(sev_carriers)}) ---")
        cooc_matrix = calculate_cooccurrence_matrix(sev_carriers)
        all_results[severity] = cooc_matrix
    
    # Save
    with pd.ExcelWriter(f"{OUTPUT_DIR}/4.2A_cooccurrence_matrices.xlsx") as writer:
        for name, matrix in all_results.items():
            if isinstance(matrix, pd.DataFrame):
                matrix.to_excel(writer, sheet_name=name[:31])
    
    cooc_stats.to_csv(f"{OUTPUT_DIR}/4.2A_cooccurrence_statistics.csv", index=False)
    print(f"\nSaved: 4.2A_cooccurrence_matrices.xlsx, 4.2A_cooccurrence_statistics.csv")
    
    return all_results


def calculate_cooccurrence_matrix(df):
    """Calculate pairwise co-occurrence counts."""
    n = len(PATHWAYS)
    matrix = np.zeros((n, n), dtype=int)
    
    for i, p1 in enumerate(PATHWAYS):
        col1 = f'pathway_{p1}'
        for j, p2 in enumerate(PATHWAYS):
            col2 = f'pathway_{p2}'
            count = ((df[col1] == 1) & (df[col2] == 1)).sum()
            matrix[i, j] = count
    
    return pd.DataFrame(matrix, index=PATHWAYS, columns=PATHWAYS)


def calculate_cooccurrence_statistics(df):
    """Calculate co-occurrence statistics with multiple metrics."""
    n_total = len(df)
    results = []
    
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i >= j:
                continue
            
            col1 = f'pathway_{p1}'
            col2 = f'pathway_{p2}'
            
            # Build 2x2 contingency table
            both = ((df[col1] == 1) & (df[col2] == 1)).sum()
            p1_only = ((df[col1] == 1) & (df[col2] == 0)).sum()
            p2_only = ((df[col1] == 0) & (df[col2] == 1)).sum()
            neither = ((df[col1] == 0) & (df[col2] == 0)).sum()
            
            n1 = both + p1_only
            n2 = both + p2_only
            
            # Jaccard similarity
            union = both + p1_only + p2_only
            jaccard = both / union if union > 0 else 0
            
            # Expected under independence
            expected = (n1 * n2) / n_total if n_total > 0 else 0
            obs_exp_ratio = both / expected if expected > 0 else np.nan
            
            # Fisher's exact test
            contingency = np.array([[both, p1_only], [p2_only, neither]])
            try:
                odds_ratio, fisher_p = stats.fisher_exact(contingency)
            except:
                odds_ratio, fisher_p = np.nan, np.nan
            
            results.append({
                'pathway_1': p1,
                'pathway_2': p2,
                'cooccurrence_count': both,
                'expected_count': round(expected, 1),
                'obs_exp_ratio': round(obs_exp_ratio, 2) if not np.isnan(obs_exp_ratio) else None,
                'jaccard': round(jaccard, 3),
                'odds_ratio': round(odds_ratio, 2) if not np.isnan(odds_ratio) else None,
                'fisher_neglog10p': round(safe_log10_pvalue(fisher_p), 1),
            })
    
    return pd.DataFrame(results)


# =============================================================================
# 4.2B: CORRELATION MATRICES (COUNT AND BURDEN)
# =============================================================================

def analyze_correlations(carriers):
    """Calculate pathway correlations for BOTH count and burden scores."""
    print("\n" + "="*70)
    print("4.2B: PATHWAY CORRELATION MATRICES (Count AND Burden)")
    print("="*70)
    
    count_results = {}
    burden_results = {}
    
    # Check if burden columns exist
    has_burden = f'pathway_{PATHWAYS[0]}_burden' in carriers.columns
    
    # COUNT-BASED CORRELATIONS
    print("\n--- COUNT-BASED Correlations ---")
    count_cols = [f'pathway_{p}_count' for p in PATHWAYS]
    
    if all(c in carriers.columns for c in count_cols):
        print("\nOverall (Count):")
        corr_count = carriers[count_cols].corr(method='spearman')
        corr_count.index = PATHWAYS
        corr_count.columns = PATHWAYS
        count_results['Overall'] = corr_count
        print(corr_count.round(3).to_string())
    
    # BURDEN-BASED CORRELATIONS
    if has_burden:
        print("\n--- BURDEN-BASED Correlations ---")
        burden_cols = [f'pathway_{p}_burden' for p in PATHWAYS]
        
        print("\nOverall (Burden):")
        corr_burden = carriers[burden_cols].corr(method='spearman')
        corr_burden.index = PATHWAYS
        corr_burden.columns = PATHWAYS
        burden_results['Overall'] = corr_burden
        print(corr_burden.round(3).to_string())
        
        # By severity
        for severity in ['Mild', 'Moderate', 'Severe']:
            sev_carriers = carriers[carriers['severity_category'] == severity]
            if len(sev_carriers) < 30:
                continue
            
            corr_burden_sev = sev_carriers[burden_cols].corr(method='spearman')
            corr_burden_sev.index = PATHWAYS
            corr_burden_sev.columns = PATHWAYS
            burden_results[severity] = corr_burden_sev
    
    # Save
    if count_results:
        with pd.ExcelWriter(f"{OUTPUT_DIR}/4.2B_correlation_matrices_count.xlsx") as writer:
            for name, matrix in count_results.items():
                matrix.to_excel(writer, sheet_name=name)
    
    if burden_results:
        with pd.ExcelWriter(f"{OUTPUT_DIR}/4.2B_correlation_matrices_burden.xlsx") as writer:
            for name, matrix in burden_results.items():
                matrix.to_excel(writer, sheet_name=name)
    
    print(f"\nSaved: 4.2B_correlation_matrices_count.xlsx, 4.2B_correlation_matrices_burden.xlsx")
    
    return count_results, burden_results


# =============================================================================
# 4.3: FREQUENCY VS INTENSITY COMPARISON
# =============================================================================

def analyze_frequency_vs_intensity(carriers):
    """Compare pathway metrics across all scoring methods."""
    print("\n" + "="*70)
    print("4.3: FREQUENCY VS INTENSITY (All Scoring Methods)")
    print("="*70)
    
    has_burden = f'pathway_{PATHWAYS[0]}_burden' in carriers.columns
    
    results = []
    
    for severity in ['Mild', 'Moderate', 'Severe']:
        sev_carriers = carriers[carriers['severity_category'] == severity]
        n = len(sev_carriers)
        if n == 0:
            continue
            
        print(f"\n--- {severity} (n={n}) ---")
        
        for pathway in PATHWAYS:
            binary_col = f'pathway_{pathway}'
            count_col = f'pathway_{pathway}_count'
            
            # Frequency
            freq = 100 * sev_carriers[binary_col].sum() / n
            
            # Count-based intensity
            affected = sev_carriers[sev_carriers[binary_col] == 1]
            mean_count = affected[count_col].mean() if len(affected) > 0 else 0
            
            row = {
                'severity': severity,
                'pathway': pathway,
                'n_carriers': n,
                'frequency_pct': round(freq, 1),
                'mean_count': round(mean_count, 2),
            }
            
            # Burden-based (if available)
            if has_burden:
                burden_col = f'pathway_{pathway}_burden'
                weighted_col = f'pathway_{pathway}_weighted'
                interaction_col = f'pathway_{pathway}_interaction'
                
                row['mean_weighted'] = round(affected[weighted_col].mean(), 3) if len(affected) > 0 else 0
                row['mean_interaction'] = round(affected[interaction_col].mean(), 3) if len(affected) > 0 else 0
                row['mean_burden'] = round(affected[burden_col].mean(), 3) if len(affected) > 0 else 0
            
            results.append(row)
            
            burden_str = f", burden={row.get('mean_burden', 'N/A')}" if has_burden else ""
            print(f"  {pathway}: {freq:.1f}% freq, count={mean_count:.2f}{burden_str}")
    
    freq_int_df = pd.DataFrame(results)
    freq_int_df.to_csv(f"{OUTPUT_DIR}/4.3_frequency_vs_intensity.csv", index=False)
    print(f"\nSaved: 4.3_frequency_vs_intensity.csv")
    
    return freq_int_df


# =============================================================================
# 4.4: NETWORK VISUALIZATION DATA (4 EDGE TYPES)
# =============================================================================

def generate_network_data(carriers, cooc_matrices, burden_correlations):
    """Generate network with 4 edge types: co-occurrence, Jaccard, correlation, interaction."""
    print("\n" + "="*70)
    print("4.4: NETWORK VISUALIZATION DATA (4 Edge Types)")
    print("="*70)
    
    n_total = len(carriers)
    has_burden = f'pathway_{PATHWAYS[0]}_burden' in carriers.columns
    
    # Nodes
    nodes = []
    for pathway in PATHWAYS:
        col = f'pathway_{pathway}'
        count = carriers[col].sum()
        
        node = {
            'id': pathway,
            'label': pathway.replace('_', ' '),
            'size': int(count),
            'pct': round(100 * count / n_total, 1),
        }
        
        # Add mean burden if available
        if has_burden:
            burden_col = f'pathway_{pathway}_burden'
            affected = carriers[carriers[col] == 1]
            node['mean_burden'] = round(affected[burden_col].mean(), 3) if len(affected) > 0 else 0
        
        nodes.append(node)
    
    # Edges with 4 weight types
    edges = []
    cooc_overall = cooc_matrices['Overall']
    
    # Get burden correlation matrix if available
    burden_corr = burden_correlations.get('Overall', None) if burden_correlations else None
    
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i >= j:
                continue
            
            count = cooc_overall.loc[p1, p2]
            if count == 0:
                continue
            
            n1 = cooc_overall.loc[p1, p1]
            n2 = cooc_overall.loc[p2, p2]
            
            # (a) Co-occurrence count
            cooccurrence = int(count)
            
            # (b) Jaccard similarity
            jaccard = count / (n1 + n2 - count) if (n1 + n2 - count) > 0 else 0
            
            # (c) Burden correlation (if available)
            if burden_corr is not None:
                burden_correlation = burden_corr.loc[p1, p2]
            else:
                burden_correlation = np.nan
            
            # (d) Mean interaction multiplier when pathways co-occur
            # Get patients with both pathways and calculate mean interaction
            col1 = f'pathway_{p1}'
            col2 = f'pathway_{p2}'
            both_affected = carriers[(carriers[col1] == 1) & (carriers[col2] == 1)]
            
            if len(both_affected) > 0 and 'combined_interaction_multiplier' in both_affected.columns:
                mean_interaction = both_affected['combined_interaction_multiplier'].mean()
            else:
                mean_interaction = 1.0
            
            edge = {
                'source': p1,
                'target': p2,
                # 4 edge weight types
                'cooccurrence_count': cooccurrence,
                'jaccard': round(jaccard, 3),
                'burden_correlation': round(burden_correlation, 3) if not np.isnan(burden_correlation) else None,
                'mean_interaction_multiplier': round(mean_interaction, 3),
                # Interpretation
                'cooccurrence_strength': 'strong' if cooccurrence > 500 else ('moderate' if cooccurrence > 100 else 'weak'),
                'correlation_strength': interpret_effect_size(burden_correlation, 'correlation') if not np.isnan(burden_correlation) else 'N/A',
            }
            
            edges.append(edge)
    
    network_data = {'nodes': nodes, 'edges': edges}
    
    # Type conversion for JSON
    def convert_types(obj):
        if isinstance(obj, dict):
            return {k: convert_types(v) for k, v in obj.items()}
        elif isinstance(obj, list):
            return [convert_types(i) for i in obj]
        elif isinstance(obj, (np.integer, np.int64)):
            return int(obj)
        elif isinstance(obj, (np.floating, np.float64)):
            return float(obj) if not np.isnan(obj) else None
        return obj
    
    # Save
    with open(f"{OUTPUT_DIR}/4.4_network_data.json", 'w') as f:
        json.dump(convert_types(network_data), f, indent=2)
    
    pd.DataFrame(nodes).to_csv(f"{OUTPUT_DIR}/4.4_network_nodes.csv", index=False)
    pd.DataFrame(edges).to_csv(f"{OUTPUT_DIR}/4.4_network_edges.csv", index=False)
    
    print(f"\nNodes: {len(nodes)}")
    print(f"Edges: {len(edges)}")
    
    print("\n--- Edge Weight Types ---")
    print("  (a) cooccurrence_count: Number of patients with both pathways")
    print("  (b) jaccard: Jaccard similarity coefficient")
    print("  (c) burden_correlation: Spearman correlation of burden scores")
    print("  (d) mean_interaction_multiplier: Mean gene-gene interaction multiplier when co-occurring")
    
    print("\n--- Top 5 by each metric ---")
    edges_df = pd.DataFrame(edges)
    
    print("\nBy Co-occurrence:")
    for _, e in edges_df.nlargest(5, 'cooccurrence_count').iterrows():
        print(f"  {e['source']} ↔ {e['target']}: {e['cooccurrence_count']}")
    
    print("\nBy Jaccard:")
    for _, e in edges_df.nlargest(5, 'jaccard').iterrows():
        print(f"  {e['source']} ↔ {e['target']}: {e['jaccard']:.3f}")
    
    if burden_corr is not None:
        print("\nBy Burden Correlation:")
        valid_corr = edges_df[edges_df['burden_correlation'].notna()]
        for _, e in valid_corr.nlargest(5, 'burden_correlation').iterrows():
            print(f"  {e['source']} ↔ {e['target']}: {e['burden_correlation']:.3f}")
    
    print("\nBy Mean Interaction Multiplier:")
    for _, e in edges_df.nlargest(5, 'mean_interaction_multiplier').iterrows():
        print(f"  {e['source']} ↔ {e['target']}: {e['mean_interaction_multiplier']:.3f}")
    
    print(f"\nSaved: 4.4_network_data.json, 4.4_network_nodes.csv, 4.4_network_edges.csv")
    
    return network_data


# =============================================================================
# 4.5: CROSS-SEVERITY COMPARISON (with effect sizes)
# =============================================================================

def analyze_cross_severity(carriers):
    """Statistical comparison of pathways across severity levels."""
    print("\n" + "="*70)
    print("4.5: CROSS-SEVERITY COMPARISON (Count AND Burden)")
    print("="*70)
    
    has_burden = f'pathway_{PATHWAYS[0]}_burden' in carriers.columns
    
    results = []
    
    for pathway in PATHWAYS:
        binary_col = f'pathway_{pathway}'
        count_col = f'pathway_{pathway}_count'
        
        # Prevalence by severity
        prev = {}
        for severity in ['Mild', 'Moderate', 'Severe']:
            sev_df = carriers[carriers['severity_category'] == severity]
            prev[severity] = sev_df[binary_col].mean() if len(sev_df) > 0 else 0
        
        # Chi-square test (binary)
        contingency_data = []
        for severity in ['Mild', 'Moderate', 'Severe']:
            sev_df = carriers[carriers['severity_category'] == severity]
            has_pathway = sev_df[binary_col].sum()
            no_pathway = len(sev_df) - has_pathway
            contingency_data.append([has_pathway, no_pathway])
        
        contingency = np.array(contingency_data)
        if contingency.min() >= 5:
            chi2, chi_p, dof, expected = stats.chi2_contingency(contingency)
            cramers = cramers_v(contingency)
        else:
            chi2, chi_p, cramers = np.nan, np.nan, np.nan
        
        # Kruskal-Wallis test (COUNT-based)
        groups_count = []
        for severity in ['Mild', 'Moderate', 'Severe']:
            sev_df = carriers[carriers['severity_category'] == severity]
            if len(sev_df) > 0:
                groups_count.append(sev_df[count_col].values)
        
        if len(groups_count) >= 2:
            h_count, kw_p_count = stats.kruskal(*groups_count)
            eps_count = epsilon_squared(h_count, len(carriers))
        else:
            h_count, kw_p_count, eps_count = np.nan, np.nan, np.nan
        
        row = {
            'pathway': pathway,
            'prevalence_mild_pct': round(100 * prev['Mild'], 1),
            'prevalence_moderate_pct': round(100 * prev['Moderate'], 1),
            'prevalence_severe_pct': round(100 * prev['Severe'], 1),
            # Chi-square
            'chi2_statistic': round(chi2, 2) if not np.isnan(chi2) else None,
            'chi2_neglog10p': round(safe_log10_pvalue(chi_p), 1) if not np.isnan(chi_p) else None,
            'cramers_v': round(cramers, 3) if not np.isnan(cramers) else None,
            # Kruskal-Wallis (COUNT)
            'kw_count_H': round(h_count, 2) if not np.isnan(h_count) else None,
            'kw_count_neglog10p': round(safe_log10_pvalue(kw_p_count), 1) if not np.isnan(kw_p_count) else None,
            'kw_count_epsilon_sq': round(eps_count, 4) if not np.isnan(eps_count) else None,
        }
        
        # Kruskal-Wallis test (BURDEN-based)
        if has_burden:
            burden_col = f'pathway_{pathway}_burden'
            groups_burden = []
            for severity in ['Mild', 'Moderate', 'Severe']:
                sev_df = carriers[carriers['severity_category'] == severity]
                if len(sev_df) > 0:
                    groups_burden.append(sev_df[burden_col].values)
            
            if len(groups_burden) >= 2:
                h_burden, kw_p_burden = stats.kruskal(*groups_burden)
                eps_burden = epsilon_squared(h_burden, len(carriers))
            else:
                h_burden, kw_p_burden, eps_burden = np.nan, np.nan, np.nan
            
            row['kw_burden_H'] = round(h_burden, 2) if not np.isnan(h_burden) else None
            row['kw_burden_neglog10p'] = round(safe_log10_pvalue(kw_p_burden), 1) if not np.isnan(kw_p_burden) else None
            row['kw_burden_epsilon_sq'] = round(eps_burden, 4) if not np.isnan(eps_burden) else None
        
        results.append(row)
        
        print(f"\n{pathway}:")
        print(f"  Prevalence: Mild={row['prevalence_mild_pct']}%, Mod={row['prevalence_moderate_pct']}%, Sev={row['prevalence_severe_pct']}%")
        print(f"  Chi-square: Cramér's V = {row['cramers_v']} ({interpret_effect_size(cramers, 'cramers_v')})")
        print(f"  K-W (count): ε² = {row['kw_count_epsilon_sq']} ({interpret_effect_size(eps_count, 'epsilon_sq')})")
        if has_burden:
            print(f"  K-W (burden): ε² = {row.get('kw_burden_epsilon_sq')} ({interpret_effect_size(row.get('kw_burden_epsilon_sq'), 'epsilon_sq')})")
    
    comparison_df = pd.DataFrame(results)
    comparison_df.to_csv(f"{OUTPUT_DIR}/4.5_cross_severity_comparison.csv", index=False)
    print(f"\nSaved: 4.5_cross_severity_comparison.csv")
    
    return comparison_df


# =============================================================================
# 4.6: VALIDATION - COUNT VS BURDEN SCORING
# =============================================================================

def validation_count_vs_burden(carriers):
    """Compare predictive power of count-based vs burden-based scoring."""
    print("\n" + "="*70)
    print("4.6: VALIDATION - Count vs Burden Scoring Comparison")
    print("="*70)
    
    has_burden = 'composite_score' in carriers.columns
    
    if not has_burden:
        print("\n⚠️  Burden columns not found. Skipping validation.")
        return None
    
    # Create severity numeric
    carriers_valid = carriers[carriers['severity_category'].isin(['Mild', 'Moderate', 'Severe'])].copy()
    carriers_valid['severity_numeric'] = carriers_valid['severity_category'].map(SEVERITY_NUMERIC)
    
    results = []
    
    # -------------------------------------------------------------------------
    # Test 1: Correlation with clinical severity
    # -------------------------------------------------------------------------
    print("\n--- Test 1: Correlation with Clinical Severity Score ---")
    
    if 'severity_score' in carriers_valid.columns:
        # Count-based: sum of pathway counts
        count_cols = [f'pathway_{p}_count' for p in PATHWAYS]
        carriers_valid['total_count'] = carriers_valid[count_cols].sum(axis=1)
        
        # Correlation: count vs severity
        corr_count, p_count = stats.spearmanr(
            carriers_valid['total_count'], 
            carriers_valid['severity_score'],
            nan_policy='omit'
        )
        
        # Correlation: burden vs severity
        corr_burden, p_burden = stats.spearmanr(
            carriers_valid['composite_score'], 
            carriers_valid['severity_score'],
            nan_policy='omit'
        )
        
        print(f"\n  Total Count vs Severity Score:")
        print(f"    Spearman r = {corr_count:.3f}, p = {p_count:.2e}")
        print(f"    Interpretation: {interpret_effect_size(corr_count, 'correlation')}")
        
        print(f"\n  Composite Burden vs Severity Score:")
        print(f"    Spearman r = {corr_burden:.3f}, p = {p_burden:.2e}")
        print(f"    Interpretation: {interpret_effect_size(corr_burden, 'correlation')}")
        
        improvement = corr_burden - corr_count
        print(f"\n  Improvement (burden over count): {improvement:+.3f}")
        
        results.append({
            'test': 'Correlation with severity_score',
            'metric': 'Spearman r',
            'count_based': round(corr_count, 3),
            'burden_based': round(corr_burden, 3),
            'improvement': round(improvement, 3),
            'better_method': 'burden' if corr_burden > corr_count else 'count'
        })
    
    # -------------------------------------------------------------------------
    # Test 2: Discrimination between severity categories (AUC)
    # -------------------------------------------------------------------------
    print("\n--- Test 2: Discrimination - Mild vs Severe (AUC) ---")
    
    mild_severe = carriers_valid[carriers_valid['severity_category'].isin(['Mild', 'Severe'])].copy()
    if len(mild_severe) > 50:
        y_true = (mild_severe['severity_category'] == 'Severe').astype(int)
        
        # AUC for count
        auc_count = roc_auc_score(y_true, mild_severe['total_count'])
        
        # AUC for burden
        auc_burden = roc_auc_score(y_true, mild_severe['composite_score'])
        
        print(f"\n  Total Count AUC: {auc_count:.3f}")
        print(f"  Composite Burden AUC: {auc_burden:.3f}")
        print(f"  Improvement: {auc_burden - auc_count:+.3f}")
        
        results.append({
            'test': 'AUC Mild vs Severe',
            'metric': 'ROC AUC',
            'count_based': round(auc_count, 3),
            'burden_based': round(auc_burden, 3),
            'improvement': round(auc_burden - auc_count, 3),
            'better_method': 'burden' if auc_burden > auc_count else 'count'
        })
    
    # -------------------------------------------------------------------------
    # Test 3: Mean scores by severity category
    # -------------------------------------------------------------------------
    print("\n--- Test 3: Mean Scores by Severity Category ---")
    
    print("\n  Count-based (total_count):")
    for sev in ['Mild', 'Moderate', 'Severe']:
        mean_val = carriers_valid[carriers_valid['severity_category'] == sev]['total_count'].mean()
        print(f"    {sev}: {mean_val:.2f}")
    
    print("\n  Burden-based (composite_score):")
    for sev in ['Mild', 'Moderate', 'Severe']:
        mean_val = carriers_valid[carriers_valid['severity_category'] == sev]['composite_score'].mean()
        print(f"    {sev}: {mean_val:.3f}")
    
    # Effect size: difference between Mild and Severe
    mild_count = carriers_valid[carriers_valid['severity_category'] == 'Mild']['total_count']
    severe_count = carriers_valid[carriers_valid['severity_category'] == 'Severe']['total_count']
    mild_burden = carriers_valid[carriers_valid['severity_category'] == 'Mild']['composite_score']
    severe_burden = carriers_valid[carriers_valid['severity_category'] == 'Severe']['composite_score']
    
    # Cohen's d for count
    pooled_std_count = np.sqrt((mild_count.std()**2 + severe_count.std()**2) / 2)
    cohens_d_count = (severe_count.mean() - mild_count.mean()) / pooled_std_count if pooled_std_count > 0 else 0
    
    # Cohen's d for burden
    pooled_std_burden = np.sqrt((mild_burden.std()**2 + severe_burden.std()**2) / 2)
    cohens_d_burden = (severe_burden.mean() - mild_burden.mean()) / pooled_std_burden if pooled_std_burden > 0 else 0
    
    print(f"\n  Cohen's d (Mild vs Severe):")
    print(f"    Count-based: {cohens_d_count:.3f}")
    print(f"    Burden-based: {cohens_d_burden:.3f}")
    
    results.append({
        'test': 'Cohen d Mild vs Severe',
        'metric': 'Cohen d',
        'count_based': round(cohens_d_count, 3),
        'burden_based': round(cohens_d_burden, 3),
        'improvement': round(cohens_d_burden - cohens_d_count, 3),
        'better_method': 'burden' if cohens_d_burden > cohens_d_count else 'count'
    })
    
    # -------------------------------------------------------------------------
    # Test 4: Interaction effects
    # -------------------------------------------------------------------------
    print("\n--- Test 4: Interaction Effects on Severity ---")
    
    if 'combined_interaction_multiplier' in carriers_valid.columns:
        # Does interaction multiplier correlate with severity?
        corr_int, p_int = stats.spearmanr(
            carriers_valid['combined_interaction_multiplier'],
            carriers_valid['severity_score'],
            nan_policy='omit'
        )
        
        print(f"\n  Interaction Multiplier vs Severity Score:")
        print(f"    Spearman r = {corr_int:.3f}, p = {p_int:.2e}")
        
        # Mean interaction by severity
        print("\n  Mean Interaction Multiplier by Severity:")
        for sev in ['Mild', 'Moderate', 'Severe']:
            mean_int = carriers_valid[carriers_valid['severity_category'] == sev]['combined_interaction_multiplier'].mean()
            print(f"    {sev}: {mean_int:.3f}")
        
        results.append({
            'test': 'Interaction multiplier vs severity',
            'metric': 'Spearman r',
            'count_based': None,
            'burden_based': round(corr_int, 3),
            'improvement': None,
            'better_method': 'N/A'
        })
    
    # Save results
    results_df = pd.DataFrame(results)
    results_df.to_csv(f"{OUTPUT_DIR}/4.6_validation_count_vs_burden.csv", index=False)
    
    # Summary
    summary_lines = [
        "="*70,
        "VALIDATION SUMMARY: Count vs Burden Scoring",
        "="*70,
        "",
    ]
    
    n_burden_wins = sum(1 for r in results if r.get('better_method') == 'burden')
    n_count_wins = sum(1 for r in results if r.get('better_method') == 'count')
    
    summary_lines.append(f"Tests where BURDEN performed better: {n_burden_wins}")
    summary_lines.append(f"Tests where COUNT performed better: {n_count_wins}")
    summary_lines.append("")
    
    if n_burden_wins > n_count_wins:
        summary_lines.append("CONCLUSION: Weighted burden scoring shows improvement over simple counts.")
        summary_lines.append("            The evidence weights and interaction multipliers add predictive value.")
    elif n_count_wins > n_burden_wins:
        summary_lines.append("CONCLUSION: Simple gene counts perform as well or better than weighted scoring.")
        summary_lines.append("            Consider whether the added complexity is justified.")
    else:
        summary_lines.append("CONCLUSION: Results are mixed. Both methods have merit for different analyses.")
    
    summary_lines.append("")
    summary_lines.append("Detailed results in: 4.6_validation_count_vs_burden.csv")
    
    summary_text = "\n".join(summary_lines)
    print(f"\n{summary_text}")
    
    with open(f"{OUTPUT_DIR}/4.6_validation_summary.txt", 'w') as f:
        f.write(summary_text)
    
    print(f"\nSaved: 4.6_validation_count_vs_burden.csv, 4.6_validation_summary.txt")
    
    return results_df


# =============================================================================
# 4.7: INTERACTION EFFECT ANALYSIS
# =============================================================================

def analyze_interaction_effects(carriers):
    """Analyze which gene-gene interactions impact clinical outcomes."""
    print("\n" + "="*70)
    print("4.7: INTERACTION EFFECT ANALYSIS")
    print("="*70)
    
    # -------------------------------------------------------------------------
    # Population-level: Which interactions are most common and impactful?
    # -------------------------------------------------------------------------
    print("\n--- Population-Level Interaction Analysis ---")
    
    # Count patients with each interaction
    interaction_counts = {}
    interaction_severity = {}
    
    for (g1, g2), info in GENE_INTERACTIONS.items():
        pair_name = f"{g1}-{g2}"
        
        # Find patients with this interaction
        patients_with_pair = carriers[carriers['interaction_pairs'].str.contains(pair_name, na=False)]
        n_patients = len(patients_with_pair)
        
        if n_patients > 0:
            interaction_counts[pair_name] = n_patients
            
            # Mean severity for patients with this interaction
            if 'severity_score' in patients_with_pair.columns:
                mean_severity = patients_with_pair['severity_score'].mean()
            else:
                mean_severity = np.nan
            
            # Mean severity for patients WITHOUT this interaction (for comparison)
            patients_without = carriers[~carriers['interaction_pairs'].str.contains(pair_name, na=False)]
            baseline_severity = patients_without['severity_score'].mean() if 'severity_score' in patients_without.columns else np.nan
            
            interaction_severity[pair_name] = {
                'n_patients': n_patients,
                'pct_of_carriers': round(100 * n_patients / len(carriers), 1),
                'multiplier': info['multiplier'],
                'type': info['type'],
                'mechanism': info['mechanism'],
                'mean_severity': round(mean_severity, 2) if not np.isnan(mean_severity) else None,
                'baseline_severity': round(baseline_severity, 2) if not np.isnan(baseline_severity) else None,
                'severity_difference': round(mean_severity - baseline_severity, 2) if not (np.isnan(mean_severity) or np.isnan(baseline_severity)) else None
            }
    
    # Create population-level dataframe
    pop_results = []
    for pair, data in interaction_severity.items():
        pop_results.append({'interaction_pair': pair, **data})
    
    pop_df = pd.DataFrame(pop_results)
    if len(pop_df) > 0:
        pop_df = pop_df.sort_values('n_patients', ascending=False)
    
    print(f"\nInteractions found in population: {len(pop_df)}")
    
    if len(pop_df) > 0:
        print("\nTop interactions by frequency:")
        for _, row in pop_df.head(10).iterrows():
            sev_diff = row['severity_difference']
            sev_str = f", severity diff={sev_diff:+.2f}" if sev_diff is not None else ""
            print(f"  {row['interaction_pair']}: {row['n_patients']} patients ({row['pct_of_carriers']}%), "
                  f"mult={row['multiplier']}{sev_str}")
    
    pop_df.to_csv(f"{OUTPUT_DIR}/4.7_interaction_effects_population.csv", index=False)
    
    # -------------------------------------------------------------------------
    # By Severity: Do synergistic interactions correlate with worse outcomes?
    # -------------------------------------------------------------------------
    print("\n--- Interaction Effects by Severity Category ---")
    
    severity_results = []
    
    for severity in ['Mild', 'Moderate', 'Severe']:
        sev_carriers = carriers[carriers['severity_category'] == severity]
        n_total = len(sev_carriers)
        
        if n_total == 0:
            continue
        
        # Count interactions
        n_with_interactions = (sev_carriers['n_interactions'] > 0).sum()
        n_with_synergy = (sev_carriers['n_synergies'] > 0).sum()
        n_with_redundancy = (sev_carriers['n_redundancies'] > 0).sum()
        
        mean_interaction_mult = sev_carriers['combined_interaction_multiplier'].mean()
        
        severity_results.append({
            'severity': severity,
            'n_carriers': n_total,
            'n_with_interactions': n_with_interactions,
            'pct_with_interactions': round(100 * n_with_interactions / n_total, 1),
            'n_with_synergy': n_with_synergy,
            'pct_with_synergy': round(100 * n_with_synergy / n_total, 1),
            'n_with_redundancy': n_with_redundancy,
            'pct_with_redundancy': round(100 * n_with_redundancy / n_total, 1),
            'mean_interaction_multiplier': round(mean_interaction_mult, 3),
        })
        
        print(f"\n  {severity} (n={n_total}):")
        print(f"    With any interaction: {n_with_interactions} ({100*n_with_interactions/n_total:.1f}%)")
        print(f"    With synergy: {n_with_synergy} ({100*n_with_synergy/n_total:.1f}%)")
        print(f"    With redundancy: {n_with_redundancy} ({100*n_with_redundancy/n_total:.1f}%)")
        print(f"    Mean interaction multiplier: {mean_interaction_mult:.3f}")
    
    severity_df = pd.DataFrame(severity_results)
    severity_df.to_csv(f"{OUTPUT_DIR}/4.7_interaction_effects_by_severity.csv", index=False)
    
    # Statistical test: is synergy rate different across severity?
    if len(severity_df) >= 3:
        print("\n--- Statistical Test: Synergy Rate vs Severity ---")
        
        contingency_data = []
        for _, row in severity_df.iterrows():
            has_synergy = row['n_with_synergy']
            no_synergy = row['n_carriers'] - row['n_with_synergy']
            contingency_data.append([has_synergy, no_synergy])
        
        contingency = np.array(contingency_data)
        if contingency.min() >= 5:
            chi2, p, dof, expected = stats.chi2_contingency(contingency)
            cramers = cramers_v(contingency)
            print(f"  Chi-square: χ² = {chi2:.2f}, p = {p:.4f}")
            print(f"  Cramér's V = {cramers:.3f} ({interpret_effect_size(cramers, 'cramers_v')})")
        else:
            print("  Insufficient cell counts for chi-square test")
    
    print(f"\nSaved: 4.7_interaction_effects_population.csv, 4.7_interaction_effects_by_severity.csv")
    
    return pop_df, severity_df


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("STAGE 4: COMPREHENSIVE PATHWAY ANALYSIS")
    print("(Count-based AND Burden-based, with Interactions)")
    print("="*70)
    
    # Create output directory
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    
    # Load data
    df, carriers = load_data()
    
    # Run analyses
    print("\n" + "="*70)
    print("RUNNING ANALYSES")
    print("="*70)
    
    # 4.1: Prevalence
    prevalence_df, pivot_df = analyze_pathway_prevalence(df, carriers)
    
    # 4.2A: Co-occurrence
    cooc_matrices = analyze_cooccurrence(carriers)
    
    # 4.2B: Correlations (count AND burden)
    count_correlations, burden_correlations = analyze_correlations(carriers)
    
    # 4.3: Frequency vs Intensity
    freq_int_df = analyze_frequency_vs_intensity(carriers)
    
    # 4.4: Network (4 edge types)
    network_data = generate_network_data(carriers, cooc_matrices, burden_correlations)
    
    # 4.5: Cross-severity
    comparison_df = analyze_cross_severity(carriers)
    
    # 4.6: Validation
    validation_df = validation_count_vs_burden(carriers)
    
    # 4.7: Interaction effects
    interaction_pop_df, interaction_sev_df = analyze_interaction_effects(carriers)
    
    # Final summary
    print("\n" + "="*70)
    print("STAGE 4 COMPLETE - OUTPUT FILES")
    print("="*70)
    
    files = [
        "4.1_pathway_prevalence_by_severity.csv",
        "4.2A_cooccurrence_matrices.xlsx",
        "4.2A_cooccurrence_statistics.csv",
        "4.2B_correlation_matrices_count.xlsx",
        "4.2B_correlation_matrices_burden.xlsx",
        "4.3_frequency_vs_intensity.csv",
        "4.4_network_data.json",
        "4.4_network_nodes.csv",
        "4.4_network_edges.csv",
        "4.5_cross_severity_comparison.csv",
        "4.6_validation_count_vs_burden.csv",
        "4.6_validation_summary.txt",
        "4.7_interaction_effects_population.csv",
        "4.7_interaction_effects_by_severity.csv",
    ]
    
    print("\nGenerated files:")
    for f in files:
        filepath = os.path.join(OUTPUT_DIR, f)
        exists = "✓" if os.path.exists(filepath) else "✗"
        print(f"  {exists} {f}")
    
    print("\n" + "="*70)
    print("KEY FINDINGS SUMMARY")
    print("="*70)
    
    # Summarize validation results
    if validation_df is not None and len(validation_df) > 0:
        n_burden_better = sum(1 for _, r in validation_df.iterrows() if r.get('better_method') == 'burden')
        n_count_better = sum(1 for _, r in validation_df.iterrows() if r.get('better_method') == 'count')
        print(f"\nValidation: Burden better in {n_burden_better}/{len(validation_df)} tests, "
              f"Count better in {n_count_better}/{len(validation_df)} tests")
    
    # Summarize interaction findings
    if interaction_pop_df is not None and len(interaction_pop_df) > 0:
        n_synergies = len(interaction_pop_df[interaction_pop_df['type'] == 'synergy'])
        print(f"\nInteractions: {len(interaction_pop_df)} interaction types observed, "
              f"{n_synergies} synergistic")
    
    print("\n✓ Stage 4 comprehensive analysis complete!")
    
    return {
        'prevalence': prevalence_df,
        'cooccurrence': cooc_matrices,
        'correlations_count': count_correlations,
        'correlations_burden': burden_correlations,
        'freq_intensity': freq_int_df,
        'network': network_data,
        'cross_severity': comparison_df,
        'validation': validation_df,
        'interactions_population': interaction_pop_df,
        'interactions_severity': interaction_sev_df,
    }


if __name__ == "__main__":
    main()
