##################### Get the longest isoforms from a Trinity output #####################

####### Thomas D. Lewin
####### Sym Geno Evo Lab
####### Last edited 25/09/2025

from collections import defaultdict
from Bio import SeqIO

# Set input and output files here
input_file = "trinity_out_dir.Trinity.fasta"
output_file = "trinity_out_dir.Trinity_longest_isos.fasta"

# For stats
input_total_seqs = 0
output_total_seqs = 0

# Store the longest isoform per gene
longest_isoforms = {}  # gene_id -> (length, record)

# Parse the input FASTA
for record in SeqIO.parse(input_file, "fasta"):
    input_total_seqs += 1

    header_parts = record.description.split()
    trinity_id = header_parts[0] 
    gene_id = "_".join(trinity_id.split("_")[:4])  

    # Extract length from 'len=XXX'
    try:
        length_str = [part for part in header_parts if part.startswith("len=")][0]
        length = int(length_str.split("=")[1])
    except IndexError:
        print(f"Could not parse length from: {record.description}")
        continue

    # Store if longest so far
    if gene_id not in longest_isoforms or length > longest_isoforms[gene_id][0]:
        longest_isoforms[gene_id] = (length, record)

# Write longest isoforms to output file
with open(output_file, "w") as out_fh:
    for length, record in longest_isoforms.values():
        SeqIO.write(record, out_fh, "fasta")
        output_total_seqs += 1

# Print stats
print(f"Input file: {input_total_seqs} total sequences")
print(f"Input file: {len(longest_isoforms)} unique genes")
print(f"Output file: {output_total_seqs} longest isoforms (1 per gene)")
