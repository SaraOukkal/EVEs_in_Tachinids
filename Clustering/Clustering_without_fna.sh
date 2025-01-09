#!/bin/bash
#SBATCH -J Clustering
#SBATCH --partition=normal
#SBATCH --cpus-per-task=30
#SBATCH --time=10:00:00
#SBATCH -o /beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Clustering/Clustering.out
#SBATCH -e/beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Clustering/Clustering.error
#SBATCH --constraint="skylake|haswell|broadwell"
#SBATCH --exclude=pbil-deb27

mmseqs="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"
Output_dir="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/"
Scripts_dir="/beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Clustering/"

##Create DB :
$mmseqs createdb "$Output_dir"Concatenated_Candidates_Schizophora_Viruses.faa "$Output_dir"Concatenated_Candidates_Schizophora_Viruses
##Clustering :
$mmseqs cluster "$Output_dir"Concatenated_Candidates_Schizophora_Viruses "$Output_dir"Candidate_clusters "$Output_dir"Candidate_clusters_tmp --threads 30 -s 7.5 --cluster-mode 1 --cov-mode 0 -c 0.30 -e 0.001
##Create tsv :
$mmseqs createtsv "$Output_dir"Concatenated_Candidates_Schizophora_Viruses "$Output_dir"Concatenated_Candidates_Schizophora_Viruses  "$Output_dir"Candidate_clusters "$Output_dir"Candidate_clusters.tsv
