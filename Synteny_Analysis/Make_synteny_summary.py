import os
import pandas as pd
import glob
import argparse

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description="Generate a synteny summary from synteny group files and cluster FASTA files.")
    parser.add_argument("--synteny_dir", required=True, help="Directory containing synteny_groups_[Cluster].csv files.")
    parser.add_argument("--fasta_dir", required=True, help="Directory containing FASTA files for clusters.")
    parser.add_argument("--output", required=True, help="Output file path for the synteny summary.")
    return parser.parse_args()

def extract_species_from_fasta(fasta_file):
    """Extract unique species names from a cluster FASTA file."""
    species_set = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith(">"):
                parts = line.strip().split(":")
                if len(parts) >= 3:
                    species = parts[2]
                    species_set.add(species)
    return species_set

def generate_synteny_summary(synteny_dir, fasta_dir, output_file):
    """Generate summary report for synteny groups."""
    clusters = [os.path.basename(f).replace(".faa", "") for f in glob.glob(os.path.join(fasta_dir, "*.faa"))]
    synteny_results = []

    for cluster in clusters:
        fasta_file = os.path.join(fasta_dir, f"{cluster}.faa")
        synteny_group_file = os.path.join(synteny_dir, f"synteny_groups_{cluster}.csv")
        species_in_cluster = extract_species_from_fasta(fasta_file)
        total_species = len(species_in_cluster)

        # **Check if the synteny group file is empty or missing → No synteny**
        if not os.path.exists(synteny_group_file) or os.stat(synteny_group_file).st_size == 0:
            synteny_results.append([cluster, "No", "NA", "NA", total_species, "NA"])
            continue  # Skip reading

        # **Try reading the file safely, skipping errors**
        try:
            synteny_df = pd.read_csv(synteny_group_file)
        except pd.errors.EmptyDataError:
            synteny_results.append([cluster, "No", "NA", "NA", total_species, "NA"])
            continue

        # **Process synteny groups if the file contains data**
        synteny_group_number = 1  # Track group numbers for the cluster

        for _, row in synteny_df.iterrows():
            species_with_synteny = set(seq.split(":")[0] for seq in row["Sequences"].split(", "))
            species_count = len(species_with_synteny)
            percent_synteny = (species_count / total_species) * 100 if total_species > 0 else 0
            
            # Append the result for this synteny group
            synteny_results.append([
                cluster, "Yes", synteny_group_number, species_count, total_species, percent_synteny
            ])
            synteny_group_number += 1  # Increment group number

    # **Save the final summary**
    pd.DataFrame(
        synteny_results,
        columns=["Cluster", "Synteny", "Synteny_group", "Species_with_Synteny", "Total_Species", "Percent_Synteny"]
    ).to_csv(output_file, index=False)

if __name__ == "__main__":
    args = parse_arguments()
    generate_synteny_summary(args.synteny_dir, args.fasta_dir, args.output)

