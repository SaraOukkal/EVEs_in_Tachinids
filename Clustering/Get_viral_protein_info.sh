#!/bin/bash

fichier_num_accession="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/mmseqs2_candidates_Refseq_virus/Accession_viral_list.txt"
fichier_ncbi="/beegfs/data/soukkal/Thesis/Databases/Refseq_viral_db_information_09_09_23.gpff"
fichier_sortie="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/mmseqs2_candidates_Refseq_virus/Viral_protein_name_species.tsv"

grep "LOCUS" -A 1 $fichier_ncbi > $fichier_ncbi.tmp

# Lecture de chaque numéro d'accession et recherche du nom de la protéine correspondante dans le fichier NCBI
while IFS= read -r num_accession; do
    num_accession_no_suffix=$(echo "$num_accession" | sed 's/\.1$//')  # Retire le suffixe ".1"
    nom_prot=$(grep -A 1 "$num_accession_no_suffix" "$fichier_ncbi.tmp" | grep "DEFINITION" | sed 's/^DEFINITION  //; s/ \[.*//')
    nom_espece=$(grep -A 1 "$num_accession_no_suffix" "$fichier_ncbi.tmp" | grep -o '\[.*\]' | sed 's/[][]//g')
    echo -e "$num_accession\t$nom_prot\t$nom_espece" >> "$fichier_sortie"
done < "$fichier_num_accession"

rm $fichier_ncbi.tmp

echo "Résultats enregistrés dans $fichier_sortie"
