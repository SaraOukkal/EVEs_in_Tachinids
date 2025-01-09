#!/bin/bash

# Chemins vers les outils et scripts
TRIMAL="/beegfs/home/soukkal/trimAl_Linux_x86-64/trimal"
READAL="/beegfs/home/soukkal/trimAl_Linux_x86-64/readal"
RESTORE_HEADERS_SCRIPT="/beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Domestication/Restore_headers_modified.py"

# Chemins des répertoires
CLUSTER_ALIGNMENT_PATH="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Domestication/Ali_codon/"

# Fichier contenant la liste des clusters
CLUSTER_LIST="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Domestication/Ali_codon/Failed_trimming_clusters_updated.txt"

# Lire la liste des clusters et exécuter les commandes pour chaque cluster
while IFS= read -r cluster; do
    echo "Processing cluster: $cluster"

    # Définir les chemins des fichiers d'entrée et de sortie
    INPUT_ALIGNMENT_FILE="${CLUSTER_ALIGNMENT_PATH}${cluster}_NT.fna.aln"
    WITHOUT_SHIFT_FILE="${INPUT_ALIGNMENT_FILE}.without_shift"
    WITHOUT_SHIFT_FILE_MODIFIED="${INPUT_ALIGNMENT_FILE}.without_shift.modified"
    UNALIGNED_FILE="${INPUT_ALIGNMENT_FILE}.unaligned"
    OUTPUT_TRIMMED_FILE="${CLUSTER_ALIGNMENT_PATH}${cluster}_NT.fna.aln.trimmed"
    OUTPUT_TRIMMED_RESTORED_FILE="${CLUSTER_ALIGNMENT_PATH}${cluster}_NT.fna.aln.trimmed.restored"

    # Modifier les ":" en "-" 
    echo "Replacing ':' with '-' in headers for $WITHOUT_SHIFT_FILE"
    awk '{ if(substr($0,1,1)==">") { gsub(":", "-", $0); } print $0; }' "$WITHOUT_SHIFT_FILE" > "$WITHOUT_SHIFT_FILE_MODIFIED"

    # Commande Readal 
    echo "Running Readal on $WITHOUT_SHIFT_FILE_MODIFIED"
    $READAL -in "$WITHOUT_SHIFT_FILE_MODIFIED" -out "$UNALIGNED_FILE" -onlyseqs

    # Commande Trimal
    echo "Running Trimal on $WITHOUT_SHIFT_FILE_MODIFIED"
    $TRIMAL -in "$WITHOUT_SHIFT_FILE_MODIFIED" -out "$OUTPUT_TRIMMED_FILE" -automated1 -backtrans "$UNALIGNED_FILE" -resoverlap 0.30 -seqoverlap 30 -fasta -ignorestopcodon

    # Restorer les headers
    echo "Restoring headers for $OUTPUT_TRIMMED_FILE"
    python3 "$RESTORE_HEADERS_SCRIPT" "$WITHOUT_SHIFT_FILE" "$OUTPUT_TRIMMED_FILE" "$OUTPUT_TRIMMED_RESTORED_FILE"

    echo "Finished processing cluster: $cluster"
    echo "--------------------------------------"

done < "$CLUSTER_LIST"

echo "All clusters processed!"
