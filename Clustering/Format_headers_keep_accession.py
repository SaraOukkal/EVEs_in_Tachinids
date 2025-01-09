import re

def format_fasta_headers(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith('>'):
                header = re.search(r'\[protein_id=(.*?)\]', line)
                if header:
                    outfile.write(f">{header.group(1)}\n")
                else:
                    outfile.write(line.strip() + "\n")
            else:
                outfile.write(line)

# Exemple d'utilisation :
input_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Nuc_sequences/Schizophora.fna"
output_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Nuc_sequences/Schizophora_formated.fna"
format_fasta_headers(input_file, output_file)
