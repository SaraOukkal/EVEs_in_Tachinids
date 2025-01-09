import os
import argparse
import subprocess

SEQKIT_PATH = "/beegfs/home/soukkal/miniconda3/bin/seqkit"

def extract_scaffolds(fasta_file, output_dir, genome_dir):
    """
    Extract scaffolds from genomes and replace '|' with '--' in filenames and headers.
    """
    os.makedirs(output_dir, exist_ok=True)

    with open(fasta_file, 'r') as infile:
        for line in infile:
            if line.startswith(">"):
                header = line.strip()
                scaffold = header.split(":")[0][1:]  # Extract scaffold
                species = header.split(":")[-1]  # Extract species

                # Replace '|' with '--' in the scaffold and species names
                scaffold_sanitized = scaffold.replace("|", "--")
                species_sanitized = species.replace("|", "--")
                scaffold_file = os.path.join(output_dir, f"{species_sanitized}:{scaffold_sanitized}.fa")

                # Path to the genome file
                genome_path = os.path.join(genome_dir, f"{species}.fa")
                if not os.path.exists(genome_path):
                    print(f"Warning: Genome file {genome_path} not found.")
                    continue

                # Extract scaffold using seqkit and sanitize headers in the output
                cmd = f"{SEQKIT_PATH} grep -p '{scaffold}' {genome_path} | sed 's/|/--/g' > \"{scaffold_file}\""
                try:
                    subprocess.run(cmd, shell=True, check=True)
                except subprocess.CalledProcessError as e:
                    print(f"Error during execution: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract scaffolds using seqkit and sanitize filenames/headers.")
    parser.add_argument("--fasta", help="Input filtered FASTA file.", required=True)
    parser.add_argument("--output_dir", help="Output directory for extracted scaffolds.", required=True)
    parser.add_argument("--genome_dir", help="Directory containing genome FASTA files.", required=True)
    args = parser.parse_args()

    extract_scaffolds(args.fasta, args.output_dir, args.genome_dir)

