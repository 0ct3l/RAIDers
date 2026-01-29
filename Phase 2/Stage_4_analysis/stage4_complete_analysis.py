#!/usr/bin/env python3
"""
Stage 4: Complete Pathway Analysis
===================================
Generates ALL data needed for Figures 4.1-4.6 and Tables 4.1-4.3

Input: patients_with_pathways_weighted.csv
Output: JSON file with all analysis results for visualization

Analysis Groups:
- Overall (all carriers with n_pathways_affected > 0)
- Mild, Moderate, Severe (by severity_category)
"""

import pandas as pd
import numpy as np
from scipy import stats
from itertools import combinations
import json
import sys
import os

# =============================================================================
# CONFIGURATION
# =============================================================================

PATHWAYS = [
    'Proteostasis',
    'RNA_Metabolism',
    'Cytoskeletal_Axonal_Transport',
    'Mitochondrial',
    'Excitotoxicity',
    'Vesicle_Trafficking',
    'DNA_Damage',
]

PATHWAY_LABELS = {
    'Proteostasis': 'Proteostasis',
    'RNA_Metabolism': 'RNA Metabolism',
    'Cytoskeletal_Axonal_Transport': 'Cytoskeletal/Axonal',
    'Mitochondrial': 'Mitochondrial',
    'Excitotoxicity': 'Excitotoxicity',
    'Vesicle_Trafficking': 'Vesicle Trafficking',
    'DNA_Damage': 'DNA Damage',
}

# Gene-pathway mapping for Table 4.2
GENE_PATHWAY_MAP = {
    'Proteostasis': ['SOD1', 'C9ORF72', 'VCP', 'UBQLN2', 'OPTN', 'SQSTM1', 'TBK1', 'CCNF'],
    'RNA_Metabolism': ['TARDBP', 'FUS', 'MATR3', 'HNRNPA1', 'HNRNPA2B1', 'ANG', 'ELP3', 'C9ORF72'],
    'Cytoskeletal_Axonal_Transport': ['TUBA4A', 'PFN1', 'NEFH', 'DCTN1', 'KIF5A'],
    'Mitochondrial': ['SOD1', 'FUS', 'CHCHD10', 'SIGMAR1', 'ATXN2', 'C19orf12'],
    'Excitotoxicity': ['SOD1', 'C9ORF72', 'TARDBP', 'UNC13A', 'DAO'],
    'Vesicle_Trafficking': ['ALS2', 'CHMP2B', 'VAPB', 'FIG4', 'SPG11'],
    'DNA_Damage': ['NEK1', 'C21orf2', 'SETX', 'SPG11'],
}

PATHWAY_MECHANISMS = {
    'Proteostasis': 'Protein misfolding, aggregation, autophagy/proteasome dysfunction',
    'RNA_Metabolism': 'RNA processing defects, stress granule dysregulation, nuclear transport failure',
    'Cytoskeletal_Axonal_Transport': 'Microtubule instability, motor protein dysfunction, axonal transport failure',
    'Mitochondrial': 'Mitochondrial dysfunction, ROS accumulation, energy deficit',
    'Excitotoxicity': 'Glutamate toxicity, calcium dysregulation, AMPA/NMDA receptor dysfunction',
    'Vesicle_Trafficking': 'Endosomal/lysosomal dysfunction, vesicle transport defects',
    'DNA_Damage': 'DNA repair defects, R-loop accumulation, genomic instability',
}

SEVERITY_ORDER = ['Mild', 'Moderate', 'Severe']
GROUPS = ['Overall', 'Mild', 'Moderate', 'Severe']

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def safe_float(x):
    """Convert numpy types to Python native for JSON serialization."""
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    elif isinstance(x, (np.integer, np.int64)):
        return int(x)
    elif isinstance(x, (np.floating, np.float64)):
        return float(x) if not np.isnan(x) else None
    elif isinstance(x, np.ndarray):
        return x.tolist()
    return x

def get_group_data(df, group_name):
    """Get data for a specific analysis group."""
    carriers = df[df['n_pathways_affected'] > 0].copy()
    
    if group_name == 'Overall':
        return carriers
    else:
        return carriers[carriers['severity_category'] == group_name]

# =============================================================================
# ANALYSIS 4.1: PATHWAY PREVALENCE BY SEVERITY
# =============================================================================

def analyze_prevalence(df):
    """Calculate pathway prevalence and scores by severity."""
    print("Analyzing 4.1: Pathway Prevalence...")
    
    carriers = df[df['n_pathways_affected'] > 0]
    
    results = {
        'by_severity': [],
        'heatmap_data': {
            'rows': SEVERITY_ORDER,
            'cols': PATHWAYS,
            'prevalence': [],
            'mean_count': [],
            'mean_burden': [],
        }
    }
    
    for severity in SEVERITY_ORDER:
        sev_data = carriers[carriers['severity_category'] == severity]
        n = len(sev_data)
        
        if n == 0:
            continue
        
        row = {
            'severity': severity,
            'n_patients': n,
            'pathways': {}
        }
        
        prev_row = []
        count_row = []
        burden_row = []
        
        for pathway in PATHWAYS:
            binary_col = f'pathway_{pathway}'
            count_col = f'pathway_{pathway}_count'
            burden_col = f'pathway_{pathway}_burden'
            
            count = sev_data[binary_col].sum()
            prevalence = 100 * count / n
            mean_count = sev_data[count_col].mean()
            mean_burden = sev_data[burden_col].mean()
            
            # Among those with pathway
            affected = sev_data[sev_data[binary_col] == 1]
            mean_count_affected = affected[count_col].mean() if len(affected) > 0 else 0
            mean_burden_affected = affected[burden_col].mean() if len(affected) > 0 else 0
            
            row['pathways'][pathway] = {
                'count': safe_float(count),
                'prevalence': safe_float(round(prevalence, 1)),
                'mean_count': safe_float(round(mean_count, 3)),
                'mean_burden': safe_float(round(mean_burden, 3)),
                'mean_count_affected': safe_float(round(mean_count_affected, 3)),
                'mean_burden_affected': safe_float(round(mean_burden_affected, 3)),
            }
            
            prev_row.append(safe_float(round(prevalence, 1)))
            count_row.append(safe_float(round(mean_count, 3)))
            burden_row.append(safe_float(round(mean_burden, 3)))
        
        results['by_severity'].append(row)
        results['heatmap_data']['prevalence'].append(prev_row)
        results['heatmap_data']['mean_count'].append(count_row)
        results['heatmap_data']['mean_burden'].append(burden_row)
    
    return results

# =============================================================================
# ANALYSIS 4.2A: CO-OCCURRENCE MATRICES
# =============================================================================

def calculate_cooccurrence(df, group_name):
    """Calculate co-occurrence matrix for a group."""
    data = get_group_data(df, group_name)
    n_total = len(data)
    
    if n_total < 10:
        return None
    
    n_pathways = len(PATHWAYS)
    matrix = np.zeros((n_pathways, n_pathways))
    
    # Calculate pairwise co-occurrence
    for i, p1 in enumerate(PATHWAYS):
        col1 = f'pathway_{p1}'
        n1 = data[col1].sum()
        
        for j, p2 in enumerate(PATHWAYS):
            col2 = f'pathway_{p2}'
            
            if i == j:
                # Diagonal: prevalence %
                matrix[i, j] = 100 * n1 / n_total if n_total > 0 else 0
            else:
                # Off-diagonal: P(p2 | p1) = "Of those with p1, what % also have p2?"
                with_p1 = data[data[col1] == 1]
                if len(with_p1) > 0:
                    both = (with_p1[col2] == 1).sum()
                    matrix[i, j] = 100 * both / len(with_p1)
                else:
                    matrix[i, j] = 0
    
    return {
        'group': group_name,
        'n_patients': safe_float(n_total),
        'matrix': [[safe_float(round(x, 1)) for x in row] for row in matrix],
        'pathways': PATHWAYS,
    }

def analyze_cooccurrence(df):
    """Calculate co-occurrence matrices for all groups."""
    print("Analyzing 4.2A: Co-occurrence Matrices...")
    
    results = {}
    for group in GROUPS:
        result = calculate_cooccurrence(df, group)
        if result:
            results[group] = result
    
    return results

# =============================================================================
# ANALYSIS 4.2B: CORRELATION MATRICES
# =============================================================================

def calculate_correlation(df, group_name, score_type='count'):
    """Calculate correlation matrix for a group."""
    data = get_group_data(df, group_name)
    n_total = len(data)
    
    if n_total < 30:
        return None
    
    # Get score columns
    if score_type == 'count':
        cols = [f'pathway_{p}_count' for p in PATHWAYS]
    else:
        cols = [f'pathway_{p}_burden' for p in PATHWAYS]
    
    # Calculate Spearman correlations
    score_data = data[cols]
    corr_matrix = score_data.corr(method='spearman')
    
    # Also calculate p-values
    n_pathways = len(PATHWAYS)
    pval_matrix = np.zeros((n_pathways, n_pathways))
    
    for i in range(n_pathways):
        for j in range(n_pathways):
            if i == j:
                pval_matrix[i, j] = 0
            else:
                r, p = stats.spearmanr(score_data.iloc[:, i], score_data.iloc[:, j])
                pval_matrix[i, j] = p
    
    return {
        'group': group_name,
        'n_patients': safe_float(n_total),
        'score_type': score_type,
        'matrix': [[safe_float(round(x, 3)) if not np.isnan(x) else 0 for x in row] for row in corr_matrix.values],
        'pvalues': [[safe_float(round(x, 4)) if not np.isnan(x) else 1 for x in row] for row in pval_matrix],
        'pathways': PATHWAYS,
    }

def analyze_correlations(df):
    """Calculate correlation matrices for all groups."""
    print("Analyzing 4.2B: Correlation Matrices...")
    
    results = {
        'count': {},
        'burden': {}
    }
    
    for group in GROUPS:
        count_result = calculate_correlation(df, group, 'count')
        burden_result = calculate_correlation(df, group, 'burden')
        
        if count_result:
            results['count'][group] = count_result
        if burden_result:
            results['burden'][group] = burden_result
    
    return results

# =============================================================================
# ANALYSIS 4.3: CO-OCCURRENCE VS CORRELATION COMPARISON
# =============================================================================

def analyze_cooccurrence_vs_correlation(cooc_results, corr_results):
    """Compare co-occurrence (frequency) vs correlation (intensity)."""
    print("Analyzing 4.3: Co-occurrence vs Correlation...")
    
    results = {}
    
    for group in GROUPS:
        if group not in cooc_results or group not in corr_results['count']:
            continue
        
        cooc_matrix = np.array(cooc_results[group]['matrix'])
        corr_matrix = np.array(corr_results['count'][group]['matrix'])
        
        # Get all unique pairs
        pairs = []
        n = len(PATHWAYS)
        
        for i in range(n):
            for j in range(i + 1, n):
                # Average the asymmetric co-occurrence
                cooc_ij = (cooc_matrix[i, j] + cooc_matrix[j, i]) / 2
                corr_ij = corr_matrix[i, j]
                
                pairs.append({
                    'pathway1': PATHWAYS[i],
                    'pathway2': PATHWAYS[j],
                    'cooccurrence': safe_float(round(cooc_ij, 1)),
                    'correlation': safe_float(round(corr_ij, 3)),
                    'label': f"{PATHWAY_LABELS[PATHWAYS[i]]} + {PATHWAY_LABELS[PATHWAYS[j]]}"
                })
        
        # Classify patterns
        for pair in pairs:
            cooc = pair['cooccurrence']
            corr = pair['correlation']
            
            if cooc > 30 and corr > 0.5:
                pair['pattern'] = 'High Freq + High Corr'
                pair['interpretation'] = 'Strong co-dependent pathways'
            elif cooc > 30 and corr < 0.3:
                pair['pattern'] = 'High Freq + Low Corr'
                pair['interpretation'] = 'Frequently co-occur but independent intensity'
            elif cooc < 20 and corr > 0.5:
                pair['pattern'] = 'Low Freq + High Corr'
                pair['interpretation'] = 'Rare but when present, co-escalate'
            else:
                pair['pattern'] = 'Low Freq + Low Corr'
                pair['interpretation'] = 'Independent pathways'
        
        results[group] = {
            'pairs': pairs,
            'n_patients': cooc_results[group]['n_patients']
        }
    
    return results

# =============================================================================
# ANALYSIS 4.4: NETWORK DATA
# =============================================================================

def calculate_network(df, group_name, cooc_results, corr_results):
    """Generate network visualization data for a group."""
    data = get_group_data(df, group_name)
    n_total = len(data)
    
    if n_total < 10:
        return None
    
    # Nodes
    nodes = []
    for pathway in PATHWAYS:
        col = f'pathway_{pathway}'
        count = data[col].sum()
        nodes.append({
            'id': pathway,
            'label': PATHWAY_LABELS[pathway],
            'size': safe_float(count),
            'pct': safe_float(round(100 * count / n_total, 1)) if n_total > 0 else 0
        })
    
    # Edges
    edges = []
    
    cooc_matrix = np.array(cooc_results[group_name]['matrix']) if group_name in cooc_results else None
    corr_matrix = np.array(corr_results['count'][group_name]['matrix']) if group_name in corr_results['count'] else None
    
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i >= j:
                continue
            
            col1 = f'pathway_{p1}'
            col2 = f'pathway_{p2}'
            
            # Co-occurrence count
            cooc_count = ((data[col1] == 1) & (data[col2] == 1)).sum()
            
            if cooc_count == 0:
                continue
            
            # Expected under independence
            n1 = data[col1].sum()
            n2 = data[col2].sum()
            expected = (n1 * n2) / n_total if n_total > 0 else 0
            
            # Jaccard
            union = ((data[col1] == 1) | (data[col2] == 1)).sum()
            jaccard = cooc_count / union if union > 0 else 0
            
            # Odds ratio
            a = cooc_count  # both
            b = n1 - cooc_count  # p1 only
            c = n2 - cooc_count  # p2 only
            d = n_total - n1 - n2 + cooc_count  # neither
            
            if b > 0 and c > 0 and d > 0:
                odds_ratio = (a * d) / (b * c)
            else:
                odds_ratio = None
            
            # Effect size classification
            if odds_ratio:
                if odds_ratio > 3:
                    effect_size = 'large'
                elif odds_ratio > 1.5:
                    effect_size = 'medium'
                else:
                    effect_size = 'negligible'
            else:
                effect_size = 'unknown'
            
            # Correlation
            correlation = corr_matrix[i, j] if corr_matrix is not None else None
            
            edges.append({
                'source': p1,
                'target': p2,
                'weight': safe_float(cooc_count),
                'expected': safe_float(round(expected, 1)),
                'jaccard': safe_float(round(jaccard, 3)),
                'odds_ratio': safe_float(round(odds_ratio, 2)) if odds_ratio else None,
                'correlation': safe_float(round(correlation, 3)) if correlation else None,
                'effect_size': effect_size,
            })
    
    return {
        'group': group_name,
        'n_patients': safe_float(n_total),
        'nodes': nodes,
        'edges': edges,
    }

def analyze_networks(df, cooc_results, corr_results):
    """Generate network data for all groups."""
    print("Analyzing 4.4: Network Data...")
    
    results = {}
    for group in GROUPS:
        result = calculate_network(df, group, cooc_results, corr_results)
        if result:
            results[group] = result
    
    return results

# =============================================================================
# ANALYSIS 4.5: CROSS-SEVERITY COMPARISON
# =============================================================================

def analyze_cross_severity(df, cooc_results, corr_results):
    """Compare pathway relationships across severity levels."""
    print("Analyzing 4.5: Cross-Severity Comparison...")
    
    carriers = df[df['n_pathways_affected'] > 0]
    
    results = {
        'pathway_trends': [],
        'relationship_stability': [],
        'statistical_tests': []
    }
    
    # 1. Pathway prevalence trends across severity
    for pathway in PATHWAYS:
        binary_col = f'pathway_{pathway}'
        count_col = f'pathway_{pathway}_count'
        burden_col = f'pathway_{pathway}_burden'
        
        trend = {
            'pathway': pathway,
            'label': PATHWAY_LABELS[pathway],
            'by_severity': {}
        }
        
        prevalences = []
        burdens = []
        
        for severity in SEVERITY_ORDER:
            sev_data = carriers[carriers['severity_category'] == severity]
            n = len(sev_data)
            
            if n > 0:
                prev = 100 * sev_data[binary_col].sum() / n
                mean_burden = sev_data[burden_col].mean()
                
                trend['by_severity'][severity] = {
                    'prevalence': safe_float(round(prev, 1)),
                    'mean_burden': safe_float(round(mean_burden, 3)),
                    'n': safe_float(n)
                }
                
                prevalences.append(prev)
                burdens.append(mean_burden)
        
        # Trend direction
        if len(prevalences) >= 2:
            if prevalences[-1] > prevalences[0] + 5:
                trend['trend'] = 'increasing'
            elif prevalences[-1] < prevalences[0] - 5:
                trend['trend'] = 'decreasing'
            else:
                trend['trend'] = 'stable'
        
        results['pathway_trends'].append(trend)
    
    # 2. Relationship stability (which pathway pairs have consistent vs changing relationships)
    for i, p1 in enumerate(PATHWAYS):
        for j, p2 in enumerate(PATHWAYS):
            if i >= j:
                continue
            
            relationship = {
                'pathway1': p1,
                'pathway2': p2,
                'cooccurrence_by_severity': {},
                'correlation_by_severity': {},
            }
            
            coocs = []
            corrs = []
            
            for severity in SEVERITY_ORDER:
                if severity in cooc_results:
                    cooc_matrix = np.array(cooc_results[severity]['matrix'])
                    avg_cooc = (cooc_matrix[i, j] + cooc_matrix[j, i]) / 2
                    relationship['cooccurrence_by_severity'][severity] = safe_float(round(avg_cooc, 1))
                    coocs.append(avg_cooc)
                
                if severity in corr_results['count']:
                    corr_matrix = np.array(corr_results['count'][severity]['matrix'])
                    relationship['correlation_by_severity'][severity] = safe_float(round(corr_matrix[i, j], 3))
                    corrs.append(corr_matrix[i, j])
            
            # Classify stability
            if len(coocs) >= 2:
                cooc_range = max(coocs) - min(coocs)
                corr_range = max(corrs) - min(corrs) if corrs else 0
                
                if cooc_range < 10 and corr_range < 0.2:
                    relationship['stability'] = 'stable'
                elif cooc_range >= 20 or corr_range >= 0.4:
                    relationship['stability'] = 'severity-dependent'
                else:
                    relationship['stability'] = 'moderate-change'
            
            results['relationship_stability'].append(relationship)
    
    # 3. Statistical tests (chi-square for prevalence, Kruskal-Wallis for burden)
    for pathway in PATHWAYS:
        binary_col = f'pathway_{pathway}'
        burden_col = f'pathway_{pathway}_burden'
        
        test_result = {
            'pathway': pathway,
            'label': PATHWAY_LABELS[pathway],
        }
        
        # Chi-square for prevalence
        contingency = []
        for severity in SEVERITY_ORDER:
            sev_data = carriers[carriers['severity_category'] == severity]
            has = sev_data[binary_col].sum()
            no = len(sev_data) - has
            contingency.append([has, no])
        
        contingency = np.array(contingency)
        if contingency.min() >= 5:
            chi2, p_chi, dof, expected = stats.chi2_contingency(contingency)
            test_result['chi2'] = safe_float(round(chi2, 2))
            test_result['chi2_pvalue'] = safe_float(round(p_chi, 4))
            test_result['chi2_significant'] = bool(p_chi < 0.05)
        
        # Kruskal-Wallis for burden
        groups = []
        for severity in SEVERITY_ORDER:
            sev_data = carriers[carriers['severity_category'] == severity]
            burdens = sev_data[burden_col].values
            if len(burdens) > 0:
                groups.append(burdens)
        
        if len(groups) >= 2:
            h_stat, p_kw = stats.kruskal(*groups)
            test_result['kruskal_h'] = safe_float(round(h_stat, 2))
            test_result['kruskal_pvalue'] = safe_float(round(p_kw, 4))
            test_result['kruskal_significant'] = bool(p_kw < 0.05)
        
        results['statistical_tests'].append(test_result)
    
    return results

# =============================================================================
# TABLES
# =============================================================================

def generate_tables():
    """Generate Table 4.2 (pathway definitions) and base for Table 4.3."""
    print("Generating Tables...")
    
    # Table 4.2: Pathway Module Definitions
    table_4_2 = []
    for pathway in PATHWAYS:
        table_4_2.append({
            'pathway': pathway,
            'label': PATHWAY_LABELS[pathway],
            'genes': GENE_PATHWAY_MAP.get(pathway, []),
            'mechanism': PATHWAY_MECHANISMS.get(pathway, ''),
            'n_genes': len(GENE_PATHWAY_MAP.get(pathway, [])),
        })
    
    return {
        'table_4_2_pathway_definitions': table_4_2,
    }

# =============================================================================
# MAIN
# =============================================================================

def main(input_file, output_file):
    print("=" * 70)
    print("Stage 4: Complete Pathway Analysis")
    print("=" * 70)
    
    # Load data
    print(f"\nLoading: {input_file}")
    df = pd.read_csv(input_file)
    print(f"Total patients: {len(df)}")
    
    carriers = df[df['n_pathways_affected'] > 0]
    print(f"Carriers (n_pathways > 0): {len(carriers)}")
    
    # Distribution
    print("\nSeverity distribution (carriers only):")
    for sev in SEVERITY_ORDER:
        n = len(carriers[carriers['severity_category'] == sev])
        print(f"  {sev}: {n}")
    
    # Run all analyses
    prevalence_results = analyze_prevalence(df)
    cooc_results = analyze_cooccurrence(df)
    corr_results = analyze_correlations(df)
    cooc_vs_corr = analyze_cooccurrence_vs_correlation(cooc_results, corr_results)
    network_results = analyze_networks(df, cooc_results, corr_results)
    cross_severity = analyze_cross_severity(df, cooc_results, corr_results)
    tables = generate_tables()
    
    # Compile all results
    all_results = {
        'metadata': {
            'total_patients': safe_float(len(df)),
            'total_carriers': safe_float(len(carriers)),
            'pathways': PATHWAYS,
            'pathway_labels': PATHWAY_LABELS,
            'severity_levels': SEVERITY_ORDER,
            'analysis_groups': GROUPS,
        },
        'fig_4_1_prevalence': prevalence_results,
        'fig_4_2_cooccurrence': cooc_results,
        'fig_4_3_correlation': corr_results,
        'fig_4_4_cooc_vs_corr': cooc_vs_corr,
        'fig_4_5_networks': network_results,
        'fig_4_6_cross_severity': cross_severity,
        'tables': tables,
    }
    
    # Save
    print(f"\nSaving: {output_file}")
    with open(output_file, 'w') as f:
        json.dump(all_results, f, indent=2)
    
    print("\n" + "=" * 70)
    print("ANALYSIS COMPLETE")
    print("=" * 70)
    print(f"\nOutput: {output_file}")
    print("\nContents:")
    print("  - fig_4_1_prevalence: Pathway prevalence heatmap data")
    print("  - fig_4_2_cooccurrence: Co-occurrence matrices (4 groups)")
    print("  - fig_4_3_correlation: Correlation matrices (4 groups × 2 score types)")
    print("  - fig_4_4_cooc_vs_corr: Scatter plot data")
    print("  - fig_4_5_networks: Network diagram data (4 groups)")
    print("  - fig_4_6_cross_severity: Cross-severity comparison")
    print("  - tables: Pathway definitions (Table 4.2)")
    
    return all_results
home = os.path.expanduser("~")


if __name__ == "__main__":
    if len(sys.argv) >= 3:
        input_file = sys.argv[1]
        output_file = sys.argv[2]
    else:
        input_file = os.path.join(home, "Desktop", "future", "2601CMUxNVIDIA_hackathon", "ALS_Synthetic_Data", "phase2_output_weighted", "patients_with_pathways_weighted.csv")
        output_file = "stage4_analysis_results.json"
    
    main(input_file, output_file)
