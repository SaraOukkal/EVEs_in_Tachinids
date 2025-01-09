import pandas as pd
import networkx as nx
from Bio import SeqIO

# Fichiers d'entrée
bv_clusters_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/BV_candidates_clusters.tsv"
iv_clusters_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/IV_candidates_clusters.tsv"
bv_candidates_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Homology_search/braconidae_headers.txt"
iv_candidates_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Homology_search/ichneumonidae_headers.txt"
bv_fasta_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/Nudivirus_IVSPER_BV_candidates.faa"
iv_fasta_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/Nudivirus_IVSPER_IV_candidates.faa"

# Fichiers de sortie
bv_filtered_clusters_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/BV_candidates_clusters_filtered.tsv"
iv_filtered_clusters_file = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/IV_candidates_clusters_filtered.tsv"
bv_output_dir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/BV_clusters/"
iv_output_dir = "/beegfs/data/soukkal/Thesis/Hymenoptera_Project/PDV/Results/Viral_genes/Clustering/IV_clusters/"

# Fonction pour lire les fichiers et créer un DataFrame des clusters
def read_cluster_file(cluster_file):
    df = pd.read_csv(cluster_file, delimiter='\t', header=None, names=["Col1", "Col2"])
    g = nx.from_pandas_edgelist(df, source='Col1', target='Col2', create_using=nx.Graph)
    data = [[f'Cluster{i}', element] for i, component in enumerate(nx.connected_components(g), 1) for element in component]
    return pd.DataFrame(data=data, columns=['Cluster', 'Names'])

# Fonction pour filtrer les clusters contenant au moins une séquence candidate
def filter_clusters_with_candidates(cluster_df, candidates_set):
    return cluster_df.groupby('Cluster').filter(lambda group: any(group['Names'].isin(candidates_set)))

# Fonction pour écrire les clusters filtrés dans un fichier
def write_filtered_clusters(filtered_clusters, output_file):
    filtered_clusters.to_csv(output_file, sep='\t', index=False, header=False)

# Fonction pour extraire les séquences FASTA par cluster et les écrire dans des fichiers, avec gestion des doublons
def extract_fasta_by_cluster(filtered_clusters, fasta_file, output_dir):
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id not in sequences:
            sequences[record.id] = record
        else:
            # Combine sequences with the same ID
            existing_record = sequences[record.id]
            existing_record.seq += record.seq

    # Créer des fichiers FASTA pour chaque cluster
    for cluster, group in filtered_clusters.groupby('Cluster'):
        cluster_file = f"{output_dir}/cluster_{cluster}.faa"
        with open(cluster_file, "w") as out_fasta:
            for name in group['Names']:
                if name in sequences:
                    SeqIO.write(sequences[name], out_fasta, "fasta")

# Lire les fichiers de clusters
bv_clusters = read_cluster_file(bv_clusters_file)
iv_clusters = read_cluster_file(iv_clusters_file)

# Lire les listes d'identifiants candidats
bv_candidates = set(line.strip() for line in open(bv_candidates_file))
iv_candidates = set(line.strip() for line in open(iv_candidates_file))

# Filtrer les clusters qui contiennent au moins un candidat pour BV
bv_filtered_clusters = filter_clusters_with_candidates(bv_clusters, bv_candidates)
write_filtered_clusters(bv_filtered_clusters, bv_filtered_clusters_file)

# Filtrer les clusters qui contiennent au moins un candidat pour IV
iv_filtered_clusters = filter_clusters_with_candidates(iv_clusters, iv_candidates)
write_filtered_clusters(iv_filtered_clusters, iv_filtered_clusters_file)

# Extraire les séquences FASTA pour chaque cluster BV
extract_fasta_by_cluster(bv_filtered_clusters, bv_fasta_file, bv_output_dir)

# Extraire les séquences FASTA pour chaque cluster IV
extract_fasta_by_cluster(iv_filtered_clusters, iv_fasta_file, iv_output_dir)

print("Traitement terminé.")
