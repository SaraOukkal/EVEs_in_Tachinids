import os
import pandas as pd
import argparse
from collections import defaultdict

def parse_filtered_fasta(fasta_file):
    """
    Parse the filtered FASTA file to extract regions of interest.
    """
    regions_of_interest = []
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                header = line.strip()[1:]  # Remove ">"
                scaffold_info, region_strand, species = header.split(":")
                region, strand = region_strand.split("(")
                strand = strand.strip(")")  # Remove closing parenthesis
                start, end = map(int, region.split("-"))
                regions_of_interest.append({
                    "Species": species,
                    "Scaffold": scaffold_info.replace("|", "--"),
                    "Start": start,
                    "End": end,
                    "Strand": strand
                })
    return regions_of_interest

def extract_adjacent_genes(gff_file, region_start, region_end):
    """
    Extract the gene immediately upstream and downstream of the region of interest.
    """
    adjacent_genes = {"Upstream": None, "Downstream": None}
    with open(gff_file, 'r') as f:
        genes = []
        for line in f:
            if not line.startswith("#"):
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    gene_start = int(fields[3])
                    gene_end = int(fields[4])
                    attributes = fields[8]
                    uniref_id = None
                    for attr in attributes.split(";"):
                        if attr.startswith("Target_ID="):
                            uniref_id = attr.split("=")[1]
                            break
                    if uniref_id:
                        genes.append({"Start": gene_start, "End": gene_end, "UniRef_ID": uniref_id})
        
        # Sort genes by start position
        genes = sorted(genes, key=lambda x: x["Start"])
        
        # Find upstream and downstream genes
        for gene in genes:
            if gene["End"] < region_start:
                adjacent_genes["Upstream"] = gene
            elif gene["Start"] > region_end and not adjacent_genes["Downstream"]:
                adjacent_genes["Downstream"] = gene
                break
    return adjacent_genes

def compare_synteny(fasta_file, predicted_genes_dir, output_report, output_groups):
    """
    Compare the upstream and downstream genes for sequences in the same cluster.
    """
    regions_of_interest = parse_filtered_fasta(fasta_file)
    comparison_results = []
    group_map = defaultdict(list)

    for region in regions_of_interest:
        species = region["Species"]
        scaffold = region["Scaffold"]
        start = region["Start"]
        end = region["End"]

        # Format the GFF file path using the species and scaffold name
        gff_file = os.path.join(predicted_genes_dir, f"{species}:{scaffold}.gff")
        if not os.path.exists(gff_file):
            print(f"Warning: GFF file {gff_file} not found for {species}, {scaffold}.")
            continue

        adjacent_genes = extract_adjacent_genes(gff_file, start, end)
        comparison_results.append({
            "Cluster": os.path.basename(fasta_file).replace(".faa", ""),
            "Species": species,
            "Scaffold": scaffold,
            "Start": start,
            "End": end,
            "Upstream_Gene": adjacent_genes["Upstream"]["UniRef_ID"] if adjacent_genes["Upstream"] else None,
            "Downstream_Gene": adjacent_genes["Downstream"]["UniRef_ID"] if adjacent_genes["Downstream"] else None
        })

    # Group sequences with the same upstream and downstream genes
    for result in comparison_results:
        group_key = (result["Upstream_Gene"], result["Downstream_Gene"])
        group_map[group_key].append(f"{result['Species']}:{result['Scaffold']}")

    # Save comparison results
    comparison_df = pd.DataFrame(comparison_results)
    comparison_df.to_csv(output_report, index=False)

    # Save groups with at least two sequences
    group_results = []
    for group_key, sequences in group_map.items():
        if len(sequences) > 1:
            group_results.append({
                "Upstream_Gene": group_key[0],
                "Downstream_Gene": group_key[1],
                "Sequences": ", ".join(sequences)
            })

    group_df = pd.DataFrame(group_results)
    group_df.to_csv(output_groups, index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Compare synteny based on upstream and downstream genes.")
    parser.add_argument("--fasta", help="Input filtered FASTA file for a single cluster.", required=True)
    parser.add_argument("--predicted_genes_dir", help="Directory containing predicted gene files for the cluster.", required=True)
    parser.add_argument("--output", help="Output CSV file for synteny report.", required=True)
    parser.add_argument("--groups_output", help="Output CSV file for synteny groups.", required=True)
    args = parser.parse_args()

    compare_synteny(
        fasta_file=args.fasta,
        predicted_genes_dir=args.predicted_genes_dir,
        output_report=args.output,
        output_groups=args.groups_output
    )

