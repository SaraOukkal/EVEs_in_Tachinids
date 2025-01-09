# Lire les IDs de schizophora à partir du fichier IDs_schizophora_in_clusters.txt
with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/IDs_schizophora_in_clusters.txt", "r") as file:
    schizophora_ids = set(file.read().splitlines())

# Lire le fichier Candidate_clusters.m8 et créer un dictionnaire pour stocker les informations sur les clusters
clusters = {}
clusters_kept = set()
clusters_removed = set()

with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidate_clusters.m8", "r") as file:
    for line in file:
        cluster, sequence = line.strip().split("\t")
        if cluster not in clusters:
            clusters[cluster] = {"total": 0, "schizophora": 0}
        clusters[cluster]["total"] += 1
        if sequence in schizophora_ids:
            clusters[cluster]["schizophora"] += 1

        total = clusters[cluster]["total"]
        schizophora = clusters[cluster]["schizophora"]
        percentage = schizophora / total

        if percentage <= 0.9:
            clusters_kept.add(cluster)
        else:
            clusters_removed.add(cluster)

# Écrire les clusters retenus dans le fichier Candidate_clusters_kept.m8
# Écrire les clusters retirés dans le fichier Candidate_clusters_removed.m8
with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Candidate_clusters_kept.m8", "w") as kept_file, open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Candidate_clusters_removed.m8", "w") as removed_file:
    with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidate_clusters.m8", "r") as file:
        for line in file:
            cluster, _ = line.strip().split("\t")
            if cluster in clusters_kept:
                kept_file.write(line)
            elif cluster in clusters_removed:
                removed_file.write(line)

# Lire le fichier Clusters_sizes.tsv et filtrer les lignes pour garder seulement les clusters gardés ou retirés
clusters_sizes_kept = []
clusters_sizes_removed = []

with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Clusters_sizes.tsv", "r") as file:
    for line in file:
        cluster = line.strip().split("\t")[0]
        if cluster in clusters_kept:
            clusters_sizes_kept.append(line)
        elif cluster in clusters_removed:
            clusters_sizes_removed.append(line)

# Écrire les clusters gardés dans le fichier Clusters_sizes_kept.tsv
with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Clusters_sizes_kept.tsv", "w") as file:
    for line in clusters_sizes_kept:
        file.write(line)

# Écrire les clusters retirés dans le fichier Clusters_sizes_removed.tsv
with open("/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Clusters_sizes_removed.tsv", "w") as file:
    for line in clusters_sizes_removed:
        file.write(line)

# Afficher un message de confirmation
print("Fichiers Candidate_clusters_kept.m8, Candidate_clusters_removed.m8, Clusters_sizes_kept.tsv et Clusters_sizes_removed.tsv générés avec succès.")

