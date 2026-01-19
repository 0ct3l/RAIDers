#!/usr/bin/env python3
"""
Gene-Pathway-Mechanism Mapping for ALS
Extracted from RAIDers Phase 2 Literature Review

This file contains the curated gene-pathway-mechanism relationships
based on ALS-specific literature, NOT generic database annotations.

Each gene is mapped to:
- pathways: List of dysfunctional pathways (can be multiple)
- mechanism: ALS-specific disease mechanism
- pathway_detail: More specific pathway category
"""

import pandas as pd
import numpy as np

# =============================================================================
# PATHWAY DEFINITIONS
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

# =============================================================================
# GENE-PATHWAY-MECHANISM MAPPING
# Extracted from RAIDers Phase 2 Literature Review
# =============================================================================

GENE_PATHWAY_MECHANISM = {
    # =========================================================================
    # PROTEOSTASIS PATHWAY GENES
    # =========================================================================
    'SOD1': {
        'pathways': ['Proteostasis', 'Mitochondrial', 'Excitotoxicity'],
        'mechanism': 'Mutations impair conformational stability causing misfolding and toxic oligomers/inclusion bodies; misfolded aggregates in mitochondria cause vacuolation and impair ROS elimination; enhanced caspase-3 cleavage of EAAT2 impairs glutamate uptake',
        'mechanism_short': 'Protein misfolding; mitochondrial toxicity; glutamate transporter cleavage',
    },
    'C9ORF72': {
        'pathways': ['Proteostasis', 'RNA_Metabolism', 'Excitotoxicity'],
        'mechanism': 'GGGGCC repeat expansion produces toxic dipeptide repeat (DPR) proteins that inhibit proteasome and interfere with nuclear-cytoplasmic transport; repeat RNA forms intranuclear foci sequestering RNA-binding proteins and altering splicing of 1000+ transcripts; upregulates GluA1 AMPA receptors increasing calcium permeability',
        'mechanism_short': 'DPR toxicity; RNA foci; AMPA receptor dysregulation',
    },
    'VCP': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Mutations impair autophagosome maturation and proteasomal degradation of ubiquitinated proteins, leading to TDP-43 and other toxic species accumulation',
        'mechanism_short': 'Autophagosome maturation failure; TDP-43 accumulation',
    },
    'UBQLN2': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Normally directs ubiquitinated substrates to proteasome; mutations impair clearance pathway resulting in inclusions that trap p62 and TDP-43',
        'mechanism_short': 'Ubiquitin-proteasome clearance failure',
    },
    'OPTN': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Acts as autophagy receptor; mutations prevent sequestration and transport of damaged organelles and aggregated proteins toward autophagosome for degradation',
        'mechanism_short': 'Autophagy receptor dysfunction',
    },
    'SQSTM1': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Acts as autophagy receptor (p62); mutations prevent sequestration and transport of damaged organelles and aggregated proteins toward autophagosome for degradation',
        'mechanism_short': 'Autophagy receptor dysfunction (p62)',
    },
    'TBK1': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Haploinsufficiency reduces ability to activate optineurin and p62, stalling autophagic clearance of protein aggregates and damaged mitochondria',
        'mechanism_short': 'Autophagy initiation kinase deficiency',
    },
    'CCNF': {
        'pathways': ['Proteostasis'],
        'mechanism': 'Mutations affect E3 ubiquitin-ligase complex leading to aberrant protein tagging with ubiquitin and failure to clear them, particularly increasing ubiquitin-tagged TDP-43',
        'mechanism_short': 'E3 ubiquitin-ligase dysfunction',
    },
    
    # =========================================================================
    # RNA METABOLISM PATHWAY GENES
    # =========================================================================
    'TARDBP': {
        'pathways': ['RNA_Metabolism', 'Excitotoxicity'],
        'mechanism': 'Mutations cause mislocalization from nucleus to cytoplasm resulting in loss of nuclear functions (pre-mRNA splicing) and gain of cytoplasmic toxicity through persistent aggregates; pathological loss downregulates ADAR2 leading to unedited GluA2 and increased AMPA receptor calcium permeability',
        'mechanism_short': 'Nuclear-cytoplasmic mislocalization; ADAR2 downregulation',
    },
    'FUS': {
        'pathways': ['RNA_Metabolism', 'Mitochondrial'],
        'mechanism': 'Mutations cause mislocalization from nucleus to cytoplasm with loss of nuclear functions and cytoplasmic aggregate toxicity; mutant forms interact with mitochondrial ATP synthase disrupting assembly and reducing ATP synthesis',
        'mechanism_short': 'Nuclear-cytoplasmic mislocalization; ATP synthase disruption',
    },
    'MATR3': {
        'pathways': ['RNA_Metabolism'],
        'mechanism': 'Mutations impede nuclear export of mRNA, sequestering transcripts within nucleus and disrupting global gene expression regulation',
        'mechanism_short': 'mRNA nuclear export block',
    },
    'HNRNPA1': {
        'pathways': ['RNA_Metabolism'],
        'mechanism': 'Mutations in prion-like domains promote excessive stress granule assembly preventing proper RNA processing and translation',
        'mechanism_short': 'Stress granule dysregulation',
    },
    'HNRNPA2B1': {
        'pathways': ['RNA_Metabolism'],
        'mechanism': 'Mutations in prion-like domains promote excessive stress granule assembly preventing proper RNA processing and translation',
        'mechanism_short': 'Stress granule dysregulation',
    },
    'ANG': {
        'pathways': ['RNA_Metabolism'],
        'mechanism': 'Loss-of-function mutations diminish ribonucleolytic activity and nuclear localization critical for ribosomal biogenesis and motor neuron survival',
        'mechanism_short': 'Ribosomal biogenesis impairment',
    },
    'ELP3': {
        'pathways': ['RNA_Metabolism'],
        'mechanism': 'Mutations trigger proteome impairment by altering tRNA modifications resulting in motor neuron axon shortening and abnormal branching',
        'mechanism_short': 'tRNA modification defects',
    },
    
    # =========================================================================
    # CYTOSKELETAL / AXONAL TRANSPORT PATHWAY GENES
    # =========================================================================
    'TUBA4A': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism': 'Missense mutations interfere with tubulin dimerization leading to destabilized microtubule network and reduced repolymerization capability',
        'mechanism_short': 'Microtubule destabilization',
    },
    'PFN1': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism': 'Mutations disrupt ATP-mediated actin polymerization arresting axonal growth and causing fragmented mitochondria accumulation',
        'mechanism_short': 'Actin polymerization failure',
    },
    'NEFH': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism': 'Variants cause abnormal neurofilament accumulation blocking movement of molecular motors and organelles',
        'mechanism_short': 'Neurofilament accumulation',
    },
    'DCTN1': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism': 'Mutations disrupt dynein-dynactin complex preventing retrograde movement of vesicles and organelles along microtubules',
        'mechanism_short': 'Retrograde transport failure',
    },
    'KIF5A': {
        'pathways': ['Cytoskeletal_Axonal_Transport'],
        'mechanism': 'Mutations cause loss of kinesin motor function disrupting anterograde transport of mitochondria and signaling molecules to distal synaptic sites',
        'mechanism_short': 'Anterograde transport failure',
    },
    
    # =========================================================================
    # MITOCHONDRIAL / OXIDATIVE STRESS PATHWAY GENES
    # =========================================================================
    'CHCHD10': {
        'pathways': ['Mitochondrial'],
        'mechanism': 'Mutations disrupt cristae integrity and mitochondrial genome stability leading to mitochondrial network fragmentation and respiratory chain deficiency',
        'mechanism_short': 'Cristae disruption; mtDNA instability',
    },
    'SIGMAR1': {
        'pathways': ['Mitochondrial'],
        'mechanism': 'Mutations prevent binding to IP3R at mitochondria-associated ER membranes causing calcium dysregulation and inducing apoptosis',
        'mechanism_short': 'MAM calcium dysregulation',
    },
    'ATXN2': {
        'pathways': ['Mitochondrial'],
        'mechanism': 'Expanded polyglutamine stretches interact with NADPH oxidase leading to massive increase in ROS production and DNA damage',
        'mechanism_short': 'NADPH oxidase activation; ROS surge',
    },
    
    # =========================================================================
    # VESICLE TRAFFICKING PATHWAY GENES
    # =========================================================================
    'ALS2': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism': 'Encodes Rab5 factor; mutations impair endosome motility and fusion with lysosomes reducing cargo degradation and increasing motor neuron vulnerability',
        'mechanism_short': 'Endosome-lysosome fusion failure',
    },
    'CHMP2B': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism': 'Disrupts ESCRT-III complex leading to dysmorphic endosomes and failure of autophagosomes to fuse with lysosomes for clearance',
        'mechanism_short': 'ESCRT-III dysfunction',
    },
    'VAPB': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism': 'Mutations lead to ER-Golgi transport defects and formation of reticular aggregates that disrupt nuclear envelope structure',
        'mechanism_short': 'ER-Golgi transport defects',
    },
    'FIG4': {
        'pathways': ['Vesicle_Trafficking'],
        'mechanism': 'Disrupts phosphoinositide metabolism affecting endolysosomal trafficking and autophagy',
        'mechanism_short': 'Phosphoinositide metabolism defects',
    },
    'SPG11': {
        'pathways': ['Vesicle_Trafficking', 'DNA_Damage'],
        'mechanism': 'Mutations disrupt lysosomal-autophagy pathway leading to lipid accumulation and DNA damage in juvenile ALS forms',
        'mechanism_short': 'Lysosomal-autophagy disruption',
    },
    
    # =========================================================================
    # DNA DAMAGE RESPONSE PATHWAY GENES
    # =========================================================================
    'NEK1': {
        'pathways': ['DNA_Damage'],
        'mechanism': 'Interacts with C21orf2 to repair DNA double-strand breaks; mutations result in failure to repair DNA damage leading to motor neuron death',
        'mechanism_short': 'DNA double-strand break repair failure',
    },
    'C21orf2': {
        'pathways': ['DNA_Damage'],
        'mechanism': 'Interacts with NEK1 to repair DNA double-strand breaks; mutations result in failure to repair DNA damage leading to motor neuron death',
        'mechanism_short': 'DNA double-strand break repair failure',
    },
    'SETX': {
        'pathways': ['DNA_Damage'],
        'mechanism': 'Encodes DNA/RNA helicase; mutations impair ability to cope with oxidative stress-induced DNA damage',
        'mechanism_short': 'Oxidative DNA damage sensitivity',
    },
    
    # =========================================================================
    # ADDITIONAL GENES FROM REFERENCE SECTION
    # =========================================================================
    'UNC13A': {
        'pathways': ['Excitotoxicity'],
        'mechanism': 'Risk variant affecting synaptic vesicle release and glutamatergic transmission',
        'mechanism_short': 'Synaptic transmission defects',
    },
    'DAO': {
        'pathways': ['Excitotoxicity'],
        'mechanism': 'D-amino acid oxidase variants affect D-serine metabolism and NMDA receptor modulation',
        'mechanism_short': 'D-serine metabolism; NMDA modulation',
    },
    'C19orf12': {
        'pathways': ['Mitochondrial'],
        'mechanism': 'Mutations affect mitochondrial integrity and iron metabolism',
        'mechanism_short': 'Mitochondrial iron dysregulation',
    },
}


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_pathway_one_hot(gene_name):
    """
    Convert gene name to one-hot encoded pathway dictionary.
    
    Parameters
    ----------
    gene_name : str
        Gene symbol (e.g., 'SOD1', 'C9ORF72')
    
    Returns
    -------
    dict
        Dictionary with pathway columns as keys and 0/1 as values
    """
    # Initialize all pathways to 0
    one_hot = {f'pathway_{pathway}': 0 for pathway in PATHWAYS}
    
    # Get pathways for this gene
    gene_clean = str(gene_name).upper().strip()
    
    # Handle C9orf72 variations
    if gene_clean in ['C9ORF72', 'C9ORF72']:
        gene_clean = 'C9ORF72'
    
    gene_info = GENE_PATHWAY_MECHANISM.get(gene_clean, None)
    
    if gene_info:
        for pathway in gene_info['pathways']:
            col_name = f'pathway_{pathway}'
            if col_name in one_hot:
                one_hot[col_name] = 1
    else:
        one_hot['pathway_Unknown'] = 1
    
    return one_hot


def get_gene_info(gene_name):
    """
    Get full pathway and mechanism info for a gene.
    
    Returns
    -------
    dict with keys: pathways, mechanism, mechanism_short
    """
    gene_clean = str(gene_name).upper().strip()
    
    if gene_clean in ['C9ORF72', 'C9ORF72']:
        gene_clean = 'C9ORF72'
    
    return GENE_PATHWAY_MECHANISM.get(gene_clean, {
        'pathways': ['Unknown'],
        'mechanism': 'Gene not in curated ALS pathway mapping',
        'mechanism_short': 'Unknown mechanism',
    })


def get_primary_pathway(gene_name):
    """
    Get the primary (first-listed) pathway for a gene.
    """
    info = get_gene_info(gene_name)
    return info['pathways'][0]


def get_mechanism(gene_name, short=True):
    """
    Get the mechanism description for a gene.
    """
    info = get_gene_info(gene_name)
    if short:
        return info.get('mechanism_short', info.get('mechanism', 'Unknown'))
    return info.get('mechanism', 'Unknown')


def create_one_hot_matrix():
    """
    Create a complete one-hot encoded DataFrame for all genes.
    """
    rows = []
    for gene, info in GENE_PATHWAY_MECHANISM.items():
        row = {
            'gene': gene,
            'primary_pathway': info['pathways'][0],
            'all_pathways': '; '.join(info['pathways']),
            'mechanism_short': info.get('mechanism_short', ''),
            'n_pathways': len(info['pathways']),
        }
        # Add one-hot columns
        one_hot = get_pathway_one_hot(gene)
        row.update(one_hot)
        rows.append(row)
    
    df = pd.DataFrame(rows)
    return df


def create_pathway_summary():
    """
    Create a summary of genes per pathway.
    """
    pathway_genes = {p: [] for p in PATHWAYS}
    
    for gene, info in GENE_PATHWAY_MECHANISM.items():
        for pathway in info['pathways']:
            if pathway in pathway_genes:
                pathway_genes[pathway].append(gene)
    
    summary = []
    for pathway, genes in pathway_genes.items():
        summary.append({
            'pathway': pathway,
            'n_genes': len(genes),
            'genes': ', '.join(sorted(genes)),
        })
    
    return pd.DataFrame(summary)


# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print("=" * 80)
    print("GENE-PATHWAY-MECHANISM MAPPING FOR ALS")
    print("Extracted from RAIDers Phase 2 Literature Review")
    print("=" * 80)
    
    # Create one-hot matrix
    one_hot_df = create_one_hot_matrix()
    
    print(f"\n[1] Total genes mapped: {len(GENE_PATHWAY_MECHANISM)}")
    print(f"[2] Total pathways: {len(PATHWAYS)}")
    
    # Multi-pathway genes
    multi_pathway = one_hot_df[one_hot_df['n_pathways'] > 1]
    print(f"\n[3] Multi-pathway genes ({len(multi_pathway)}):")
    for _, row in multi_pathway.iterrows():
        print(f"    • {row['gene']}: {row['all_pathways']}")
    
    # Pathway summary
    print("\n[4] Genes per pathway:")
    summary = create_pathway_summary()
    for _, row in summary.iterrows():
        print(f"    • {row['pathway']}: {row['n_genes']} genes")
    
    # Save files
    one_hot_df.to_csv("gene_pathway_one_hot_matrix.csv", index=False)
    summary.to_csv("pathway_summary.csv", index=False)
    
    print("\n[5] Files saved:")
    print("    • gene_pathway_one_hot_matrix.csv")
    print("    • pathway_summary.csv")
    
    # Preview
    print("\n[6] One-Hot Matrix Preview:")
    print(one_hot_df[['gene', 'primary_pathway', 'n_pathways', 
                       'pathway_Proteostasis', 'pathway_RNA_Metabolism', 
                       'pathway_Mitochondrial']].head(10).to_string(index=False))
