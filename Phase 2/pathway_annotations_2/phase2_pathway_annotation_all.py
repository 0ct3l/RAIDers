#!/usr/bin/env python3
"""
Phase 2: Pathway Annotation Script (All Populations)
=====================================================
Annotates Phase 1 patients with pathway binary + score columns.

Input:  All client_*.csv files from Phase 1 (excluding *_carriers.csv)
Output: 
  - patients_with_pathways.csv (ALL patients combined with superpopulation column)
  - client_AFR_pathways.csv, client_AMR_pathways.csv, etc. (per population)

Usage:
  python phase2_pathway_annotation_all.py
"""

import pandas as pd
import numpy as np
from pathlib import Path
import glob
import os

# =============================================================================
# CONFIGURATION - UPDATE THIS PATH IF NEEDED
# =============================================================================

INPUT_FOLDER = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data")
OUTPUT_FOLDER = os.path.expanduser("~/Desktop/future/2601CMUxNVIDIA_hackathon/ALS_Synthetic_Data/phase2_output")

# =============================================================================
# PATHWAY MAPPING (from gene_pathway_mechanism_mapping.py)
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

GENE_PATHWAY_MECHANISM = {
    # Proteostasis
    'SOD1': {'pathways': ['Proteostasis', 'Mitochondrial', 'Excitotoxicity']},
    'C9ORF72': {'pathways': ['Proteostasis', 'RNA_Metabolism', 'Excitotoxicity']},
    'VCP': {'pathways': ['Proteostasis']},
    'UBQLN2': {'pathways': ['Proteostasis']},
    'OPTN': {'pathways': ['Proteostasis']},
    'SQSTM1': {'pathways': ['Proteostasis']},
    'TBK1': {'pathways': ['Proteostasis']},
    'CCNF': {'pathways': ['Proteostasis']},
    # RNA Metabolism
    'TARDBP': {'pathways': ['RNA_Metabolism', 'Excitotoxicity']},
    'FUS': {'pathways': ['RNA_Metabolism', 'Mitochondrial']},
    'MATR3': {'pathways': ['RNA_Metabolism']},
    'HNRNPA1': {'pathways': ['RNA_Metabolism']},
    'HNRNPA2B1': {'pathways': ['RNA_Metabolism']},
    'ANG': {'pathways': ['RNA_Metabolism']},
    'ELP3': {'pathways': ['RNA_Metabolism']},
    # Cytoskeletal/Axonal Transport
    'TUBA4A': {'pathways': ['Cytoskeletal_Axonal_Transport']},
    'PFN1': {'pathways': ['Cytoskeletal_Axonal_Transport']},
    'NEFH': {'pathways': ['Cytoskeletal_Axonal_Transport']},
    'DCTN1': {'pathways': ['Cytoskeletal_Axonal_Transport']},
    'KIF5A': {'pathways': ['Cytoskeletal_Axonal_Transport']},
    # Mitochondrial
    'CHCHD10': {'pathways': ['Mitochondrial']},
    'SIGMAR1': {'pathways': ['Mitochondrial']},
    'ATXN2': {'pathways': ['Mitochondrial']},
    'C19orf12': {'pathways': ['Mitochondrial']},
    # Vesicle Trafficking
    'ALS2': {'pathways': ['Vesicle_Trafficking']},
    'CHMP2B': {'pathways': ['Vesicle_Trafficking']},
    'VAPB': {'pathways': ['Vesicle_Trafficking']},
    'FIG4': {'pathways': ['Vesicle_Trafficking']},
    'SPG11': {'pathways': ['Vesicle_Trafficking', 'DNA_Damage']},
    # DNA Damage
    'NEK1': {'pathways': ['DNA_Damage']},
    'C21orf2': {'pathways': ['DNA_Damage']},
    'SETX': {'pathways': ['DNA_Damage']},
    # Excitotoxicity
    'UNC13A': {'pathways': ['Excitotoxicity']},
    'DAO': {'pathways': ['Excitotoxicity']},
}

# Population code to name mapping
POPULATION_NAMES = {
    'AFR': 'African',
    'AMR': 'American',
    'EAS': 'East Asian',
    'EUR': 'European',
    'SAS': 'South Asian'
}

# =============================================================================
# SCORING FUNCTIONS
# =============================================================================

def parse_genes(all_genes_str):
    """Parse the 'all_genes' column into a list of gene names."""
    if pd.isna(all_genes_str) or all_genes_str == '':
        return []
    # Handle both semicolon and comma separators
    if ';' in str(all_genes_str):
        genes = [g.strip() for g in str(all_genes_str).split(';')]
    else:
        genes = [g.strip() for g in str(all_genes_str).split(',')]
    return [g for g in genes if g]  # Remove empty strings


def get_patient_pathway_scores(patient_genes):
    """
    Calculate pathway scores for a single patient.
    
    Scoring: Count of unique genes affecting each pathway.
    A patient with SOD1 + SIGMAR1 gets Mitochondrial_score = 2 (two genes).
    
    Returns dict with:
      - pathway_{name}: binary (0/1) 
      - pathway_{name}_score: gene count
      - n_pathways_affected: count of pathways with score > 0
      - primary_pathway: pathway with highest score
    """
    unique_genes = set(patient_genes)
    
    # Track which genes contribute to each pathway
    pathway_gene_sets = {pathway: set() for pathway in PATHWAYS}
    
    for gene in unique_genes:
        gene_info = GENE_PATHWAY_MECHANISM.get(gene)
        if gene_info:
            for pathway in gene_info['pathways']:
                pathway_gene_sets[pathway].add(gene)
    
    result = {}
    
    # Binary and score columns for each pathway
    for pathway in PATHWAYS:
        gene_count = len(pathway_gene_sets[pathway])
        result[f'pathway_{pathway}'] = 1 if gene_count > 0 else 0
        result[f'pathway_{pathway}_score'] = gene_count
    
    # Summary columns
    scores = {p: len(pathway_gene_sets[p]) for p in PATHWAYS}
    result['n_pathways_affected'] = sum(1 for s in scores.values() if s > 0)
    
    if result['n_pathways_affected'] > 0:
        result['primary_pathway'] = max(scores, key=scores.get)
    else:
        result['primary_pathway'] = 'None'
    
    return result


def process_population_file(filepath):
    """Process a single population CSV file and return annotated dataframe."""
    
    print(f"\nProcessing: {filepath}")
    df = pd.read_csv(filepath)
    print(f"  Loaded {len(df)} patients")
    
    # Extract population code from filename (e.g., 'client_AFR.csv' -> 'AFR')
    filename = os.path.basename(filepath)
    pop_code = filename.replace('client_', '').replace('.csv', '')
    
    # Add superpopulation column
    df['superpopulation'] = pop_code
    
    # Identify columns to keep (non-variant columns)
    base_cols = [
        'patient_id', 'superpopulation', 'n_variants_carried', 'primary_gene',
        'severity_score', 'severity_category', 'predicted_progression', 'all_genes'
    ]
    
    # Keep only columns that exist
    keep_cols = [c for c in base_cols if c in df.columns]
    
    # Calculate pathway scores for each patient
    print("  Calculating pathway scores...")
    pathway_data = []
    
    for idx, row in df.iterrows():
        genes = parse_genes(row.get('all_genes', ''))
        scores = get_patient_pathway_scores(genes)
        pathway_data.append(scores)
    
    # Create pathway dataframe and merge
    pathway_df = pd.DataFrame(pathway_data)
    
    # Combine base columns with pathway columns
    result_df = pd.concat([df[keep_cols].reset_index(drop=True), pathway_df], axis=1)
    
    # Count carriers
    carriers = result_df[result_df['n_pathways_affected'] > 0]
    print(f"  Carriers with pathway annotations: {len(carriers)}")
    
    return result_df


def main():
    print("=" * 70)
    print("Phase 2: Pathway Annotation - All Populations")
    print("=" * 70)
    
    # Create output folder
    os.makedirs(OUTPUT_FOLDER, exist_ok=True)
    print(f"\nInput folder:  {INPUT_FOLDER}")
    print(f"Output folder: {OUTPUT_FOLDER}")
    
    # Find all client_*.csv files (exclude *_carriers.csv)
    pattern = os.path.join(INPUT_FOLDER, "client_*.csv")
    all_files = glob.glob(pattern)
    
    # Filter out carriers files
    input_files = [f for f in all_files if '_carriers' not in f]
    
    print(f"\nFound {len(input_files)} population files to process:")
    for f in sorted(input_files):
        print(f"  - {os.path.basename(f)}")
    
    if len(input_files) == 0:
        print("\nERROR: No client_*.csv files found!")
        print(f"Please check the INPUT_FOLDER path: {INPUT_FOLDER}")
        return
    
    # Process each population
    all_dfs = []
    
    for filepath in sorted(input_files):
        df = process_population_file(filepath)
        all_dfs.append(df)
        
        # Save individual population file
        pop_code = os.path.basename(filepath).replace('client_', '').replace('.csv', '')
        output_path = os.path.join(OUTPUT_FOLDER, f"client_{pop_code}_pathways.csv")
        df.to_csv(output_path, index=False)
        print(f"  Saved: {output_path}")
    
    # Combine all populations
    print("\n" + "=" * 70)
    print("Combining all populations...")
    combined_df = pd.concat(all_dfs, ignore_index=True)
    
    # Save combined file
    combined_path = os.path.join(OUTPUT_FOLDER, "patients_with_pathways.csv")
    combined_df.to_csv(combined_path, index=False)
    print(f"Saved combined file: {combined_path}")
    
    # Print summary statistics
    print("\n" + "=" * 70)
    print("SUMMARY STATISTICS")
    print("=" * 70)
    
    print(f"\nTotal patients: {len(combined_df)}")
    print(f"Total carriers (n_pathways > 0): {len(combined_df[combined_df['n_pathways_affected'] > 0])}")
    
    print("\n--- By Superpopulation ---")
    for pop in sorted(combined_df['superpopulation'].unique()):
        pop_df = combined_df[combined_df['superpopulation'] == pop]
        carriers = pop_df[pop_df['n_pathways_affected'] > 0]
        print(f"  {pop}: {len(pop_df)} total, {len(carriers)} carriers")
    
    print("\n--- Pathway Prevalence (carriers only) ---")
    carriers_df = combined_df[combined_df['n_pathways_affected'] > 0]
    for pathway in PATHWAYS:
        col = f'pathway_{pathway}'
        count = carriers_df[col].sum()
        pct = 100 * count / len(carriers_df)
        print(f"  {pathway}: {count} ({pct:.1f}%)")
    
    print("\n--- Multi-Pathway Patients ---")
    for n in range(1, 8):
        count = len(carriers_df[carriers_df['n_pathways_affected'] == n])
        if count > 0:
            pct = 100 * count / len(carriers_df)
            print(f"  {n} pathway(s): {count} ({pct:.1f}%)")
    
    print("\n--- Primary Pathway Distribution ---")
    primary_counts = carriers_df['primary_pathway'].value_counts()
    for pathway, count in primary_counts.items():
        pct = 100 * count / len(carriers_df)
        print(f"  {pathway}: {count} ({pct:.1f}%)")
    
    print("\n--- Pathway Prevalence by Severity ---")
    for sev in ['Mild', 'Moderate', 'Severe']:
        sev_df = carriers_df[carriers_df['severity_category'] == sev]
        if len(sev_df) > 0:
            print(f"\n  {sev} (n={len(sev_df)}):")
            for pathway in PATHWAYS:
                col = f'pathway_{pathway}'
                count = sev_df[col].sum()
                pct = 100 * count / len(sev_df)
                if pct > 0:
                    print(f"    {pathway}: {count} ({pct:.1f}%)")
    
    print("\n" + "=" * 70)
    print("OUTPUT FILES:")
    print("=" * 70)
    print(f"\n1. Combined (all populations):")
    print(f"   {combined_path}")
    print(f"\n2. Per-population files:")
    for pop in sorted(combined_df['superpopulation'].unique()):
        print(f"   {os.path.join(OUTPUT_FOLDER, f'client_{pop}_pathways.csv')}")
    
    print(f"\nColumns in output ({len(combined_df.columns)} total):")
    print(f"  {list(combined_df.columns)}")
    
    print("\n✓ Phase 2 annotation complete!")


if __name__ == "__main__":
    main()
