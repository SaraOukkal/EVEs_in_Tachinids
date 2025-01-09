import pandas as pd
from Bio import AlignIO

# Lire la liste des clusters
Clusters_Phy = pd.read_csv("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Domestication/Clusters_with_EVEs_list.txt", header=None)

# Filtrer les clusters avec au moins 3 séquences dans l'alignement
filtered_clusters = []

for cluster in Clusters_Phy[0]:
    alignment_file = f"/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Domestication/Codon_ali/{cluster}.fna.aln.trimmed"  # Modifier le chemin vers vos alignements
    try:
        alignment = AlignIO.read(alignment_file, "fasta")
        if len(alignment) >= 3:
            filtered_clusters.append(cluster)
    except Exception as e:
        print(f"Error reading {alignment_file}: {e}")

# Écrire la nouvelle liste des clusters dans un fichier
with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Domestication/Clusters_with_EVEs_list_min3seq.txt", "w") as f:
    for cluster in filtered_clusters:
        f.write(cluster + "\n")
