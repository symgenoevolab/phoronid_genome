##################### Quantify ALG mixing within a chromosome #####################

####### Thomas D. Lewin
####### Sym Geno Evo Lab
####### Last edited 24/09/2025

### FUNCTION
### This script counts the average size of blocks of genes from each ALG on a set of chromososomes

### SIMPLIFIED EXAMPLE
# Where A, B and C are genes from ALGs A, B and C:
# Chr AAABBBCCC: total genes = 9; transitions = 2; blocks = 3; average genes per block = 3
# Chr AAAABBBB: total genes = 8; transitions = 1; blocks = 2; average genes per block = 4
# Chr ABCABCABC: total genes = 9; transitions = 8; blocks = 9; average genes per block = 1

### INPUT
### This script takes as input a directory of files named *_coordinates.tsv
### These are produced as an output of the SyntenyFinder pipeline
### SyntenyFinder can be found at: https://github.com/symgenoevolab/SyntenyFinder

### OUTPUT
### tsv file containing: Chromosome Num_Genes Num_Transitions Avg_Genes_per_Block Dominant_ALGs (>5%)

# Import packages
import os
import pandas as pd
from glob import glob
from collections import Counter

# Collapse certain ALG variants
def normalize_alg(alg):
    if alg in {"Ea", "Eb"}:
        return "E"
    elif alg in {"A1a", "A1b"}:
        return "A1"
    elif alg in {"Qa", "Qb", "Qc", "Qd"}:
        return "Q"
    return alg

# Count transitions between ALGs
def count_transitions(algs):
    transitions = 0
    for i in range(1, len(algs)):
        if algs[i] != algs[i - 1]:
            transitions += 1
    return transitions

# Get dominant ALGs with >5% abundance
def get_dominant_algs(algs):
    count = Counter(algs)
    total = len(algs)
    dominant = []
    for alg, freq in count.items():
        proportion = freq / total
        if proportion > 0.05:
            dominant.append(f"{alg} ({proportion*100:.1f}%)")
    return ", ".join(sorted(dominant, key=lambda x: -float(x.split('(')[-1].rstrip('%)'))))

# Process a single file
def process_file(input_file):
    prefix = os.path.basename(input_file).replace("_coordinates.tsv", "")
    output_file = f"{prefix}_shuffling_stats.tsv"

    df = pd.read_csv(input_file, sep="\t", header=None,
                     names=["idx", "status", "chrom", "start", "end", "ALG"])

    # Collapse ALGs
    df["ALG"] = df["ALG"].apply(normalize_alg)

    # Sort by chromosome then start position
    df_sorted = df.sort_values(by=["chrom", "start"])

    results = []

    for chrom, group in df_sorted.groupby("chrom"):
        algs = list(group["ALG"])
        num_genes = len(algs)
        num_transitions = count_transitions(algs)
        avg_genes_per_block = num_genes / num_transitions if num_transitions != 0 else num_genes
        dominant_algs = get_dominant_algs(algs)

        results.append([
            chrom,
            num_genes,
            num_transitions,
            round(avg_genes_per_block, 2),
            dominant_algs
        ])

    out_df = pd.DataFrame(results, columns=[
        "Chromosome", "Num_Genes", "Num_Transitions",
        "Avg_Genes_per_Block", "Dominant_ALGs (>5%)"
    ])
    out_df.to_csv(output_file, sep="\t", index=False)
    print(f"Wrote {output_file}")

# Run on all *_coordinates.tsv files in current directory
if __name__ == "__main__":
    input_files = glob("*_coordinates.tsv")
    if not input_files:
        print("No *_coordinates.tsv files found.")
    for file in input_files:
        process_file(file)

### END OF SCRIPT       
