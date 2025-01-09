#!/bin/bash

# Répertoire contenant les dossiers
repertoire_de_base="/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/"
fichier_de_sortie="/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/resultats.csv"

# Créer un fichier CSV avec un en-tête
echo "Dossier Complete Single Duplicated Fragmented Missing Total" > "$fichier_de_sortie"

# Parcourir les dossiers dans le répertoire de base
for dossier in "$repertoire_de_base"/*; do
    if [ -d "$dossier" ] && [ "$(basename "$dossier")" != "busco_downloads" ]; then
        # Rechercher le fichier short_summary correspondant
        fichiers_summary=("$dossier"/short_summary.specific.diptera_odb10.*.txt)
        
        if [ ${#fichiers_summary[@]} -gt 0 ]; then
            # Prendre le premier fichier trouvé (vous pouvez ajuster cela si nécessaire)
            fichier_summary="${fichiers_summary[0]}"
            
            # Extraire les informations avec sed
            info=$(sed -nE '/C:/ s/[^0-9.]+/ /gp' "$fichier_summary")
            
            # Extraire le nom du dossier
            dossier_name=$(basename "$dossier")
            
            # Écrire les données dans le fichier de sortie
            echo "$dossier_name $info" >> "$fichier_de_sortie"
        else
            echo "Aucun fichier summary trouvé dans $dossier"
        fi
    fi
done

echo "Terminé. Les résultats ont été enregistrés dans $fichier_de_sortie"
