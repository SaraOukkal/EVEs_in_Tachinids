#!/bin/bash

# Chemins des fichiers et dossiers
candidates_file="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Candidates_info_EVEs.tsv"
te_folder="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/TE/TE_positions"
output_eves_file="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Events/EVEs_info.tsv"
tmp_dir="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/tmp"

# Créer le répertoire temporaire s'il n'existe pas
mkdir -p "$tmp_dir"

# Pour chaque fichier TE par espèce
for te_file in "$te_folder"/*_TE_no_overlaps.bed; do
    species=$(basename "$te_file" | cut -d '_' -f 1,2)
    
    # Convertir candidates_file en fichier BED
    awk -v OFS='\t' -v species="$species" 'NR>1 && $4 == species {print $5, $6, $7}' "$candidates_file" > "$tmp_dir/${species}_candidates.bed"
    
    # Utiliser bedtools pour trouver les régions de candidates_file qui n'intersectent pas avec les TE
    output_non_overlap="$tmp_dir/${species}_non_overlap.bed"
    bedtools intersect -a "$tmp_dir/${species}_candidates.bed" -b "$te_file" -v > "$output_non_overlap"
    
    # Filtrer les lignes correspondantes dans candidates_file
    grep -Fwf "$output_non_overlap" "$candidates_file" > "$tmp_dir/${species}_filtered_candidates.tsv"
done

# Concaténer tous les fichiers filtrés pour obtenir le fichier final
cat "$tmp_dir"/*_filtered_candidates.tsv > "$output_eves_file"

# Ajouter le header du fichier candidates_file au fichier final
head -n 1 "$candidates_file" > "$output_eves_file.tmp" && cat "$output_eves_file" >> "$output_eves_file.tmp" && mv "$output_eves_file.tmp" "$output_eves_file"

echo "Le traitement est terminé."

