import os
import subprocess
import argparse

def predict_genes(scaffolds_dir, output_dir, db, metaeuk_path, temp_folder, memory_limit):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(temp_folder, exist_ok=True)

    for scaffold_file in os.listdir(scaffolds_dir):
        scaffold_path = os.path.join(scaffolds_dir, scaffold_file)
        output_prefix = os.path.join(output_dir, os.path.splitext(scaffold_file)[0])

        cmd = [
            metaeuk_path,
            "easy-predict",
            scaffold_path,
            db,
            output_prefix,
            temp_folder,
            "--min-length", "33",
            "--compressed", "1",
            "--split-memory-limit", memory_limit  # Directly use the string with suffix
        ]

        subprocess.run(cmd, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Predict genes using Metaeuk.")
    parser.add_argument("--scaffolds_dir", help="Directory containing scaffold FASTA files.", required=True)
    parser.add_argument("--output_dir", help="Output directory for predicted genes.", required=True)
    parser.add_argument("--db", help="Path to Metaeuk database.", required=True)
    parser.add_argument("--metaeuk", help="Path to Metaeuk executable.", required=True)
    parser.add_argument("--temp_folder", help="Temporary folder for Metaeuk.", required=True)
    parser.add_argument("--memory_limit", help="Memory limit for Metaeuk (e.g., 64G).", default="64G")
    args = parser.parse_args()

    predict_genes(
        scaffolds_dir=args.scaffolds_dir,
        output_dir=args.output_dir,
        db=args.db,
        metaeuk_path=args.metaeuk,
        temp_folder=args.temp_folder,
        memory_limit=args.memory_limit  # Pass the string as-is
    )

