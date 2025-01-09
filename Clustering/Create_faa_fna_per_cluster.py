from Bio import SeqIO  # Importer le module SeqIO de Biopython pour travailler avec des fichiers FASTA
import os  # Importer le module os pour effectuer des opérations sur le système de fichiers

def Read_nt_file(nt_file_path):
    nt_sequences_dict = {}
    with open(nt_file_path, "r") as nt_file:
        for record in SeqIO.parse(nt_file, "fasta"):
            nt_sequences_dict[record.id] = str(record.seq)
    return nt_sequences_dict

def Read_aa_file(aa_file_path):
    aa_sequences_dict = {}
    with open(aa_file_path, "r") as aa_file:
        for record in SeqIO.parse(aa_file, "fasta"):
            # Séparation du header pour récupérer la première partie
            accession_number = record.id.split()[0]
            aa_sequences_dict[accession_number] = str(record.seq)
    return aa_sequences_dict

def filter_clusters(clusters_file_path, nt_sequences_dict):
    filtered_lines = {}
    removed_lines = {}
    with open(clusters_file_path, "r") as clusters_file:
        for line in clusters_file:
            # Lire le nom de cluster et le motif de nom de séquence à partir de chaque ligne
            cluster_name, sequence_name = line.strip().split("\t")
            # Vérifier si le nom de séquence est présent dans le dictionnaire des séquences nucléotidiques
            if sequence_name in nt_sequences_dict.keys():
                # Ajouter la ligne filtrée au dictionnaire des lignes filtrées
                if cluster_name in filtered_lines:
                    filtered_lines[cluster_name].append(sequence_name)
                else:
                    filtered_lines[cluster_name] = [sequence_name]
            else:
                # Ajouter la ligne retirée au dictionnaire des lignes retirées
                if cluster_name in removed_lines:
                    removed_lines[cluster_name].append(sequence_name)
                else:
                    removed_lines[cluster_name] = [sequence_name]
    return filtered_lines, removed_lines

def write_sequences_to_file(output_file_path, sequences_dict, sequence_names):
    with open(output_file_path, "w") as output_file:
        for sequence_name in sequence_names:
            if sequence_name in sequences_dict:
                header = sequence_name
                sequence = sequences_dict[sequence_name]
                output_file.write(f">{header}\n{sequence}\n")

def write_cluster_files(nt_sequences_dict, aa_sequences_dict, filtered_lines, output_dir):
    for cluster_name in filtered_lines.keys():
        nt_output_file_path = os.path.join(output_dir, f"{cluster_name}.fna")
        aa_output_file_path = os.path.join(output_dir, f"{cluster_name}.faa")
        write_sequences_to_file(nt_output_file_path, nt_sequences_dict, filtered_lines[cluster_name])
        write_sequences_to_file(aa_output_file_path, aa_sequences_dict, filtered_lines[cluster_name])

def write_table_files(filtered_lines, removed_lines, filtered_clusters_file_path, removed_candidates_file_path):
    with open(filtered_clusters_file_path, "w") as filtered_file:
        for cluster_name, sequences in filtered_lines.items():
            for sequence_name in sequences:
                filtered_file.write(f"{cluster_name}\t{sequence_name}\n")
    with open(removed_candidates_file_path, "w") as removed_file:
        for cluster_name, sequences in removed_lines.items():
            for sequence_name in sequences:
                removed_file.write(f"{cluster_name}\t{sequence_name}\n")

# Chemins des inputs et outputs
nt_file_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Nuc_sequences/Concatenated_Candidates_Schizophora_Viruses.fna"
aa_file_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Concatenated_Candidates_Schizophora_Viruses.faa"
clusters_file_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Candidate_clusters_kept.m8"
output_dir = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Clusters_fasta/"
filtered_clusters_file_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Candidate_clusters_kept_filtered.m8"
removed_candidates_file_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Removed_candidates_no_nuc.txt"

# Faire des dictionnaires à partir des fasta
nt_sequences_dict = Read_nt_file(nt_file_path)
print("Dictionnaire nucléotidique généré")
aa_sequences_dict = Read_aa_file(aa_file_path)
print("Dictionnaire protéique généré")

# Filtrer le tableau de clusters à partir des séquences dispo en nucléotidique
filtered_lines, removed_lines = filter_clusters(clusters_file_path, nt_sequences_dict)
print("Fichiers filtrés")

# Ecrire les tableaux de ce qu'on garde et ce qui est retiré  
write_table_files(filtered_lines, removed_lines, filtered_clusters_file_path, removed_candidates_file_path)
print("Fichier Candidate_clusters_kept_filtered.m8 généré")
print("Fichier Removed_candidates_no_nuc.txt généré")

# Ecrire fichiers faa et fna pour chaque cluster 
write_cluster_files(nt_sequences_dict, aa_sequences_dict, filtered_lines, output_dir)
print("Fichiers de clusters écrits")
