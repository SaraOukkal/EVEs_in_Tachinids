import os
import pandas as pd


###### Partie 1 : Distribution des profondeur de séquençage des scaffolds contenant des gènes BUSCO :

# Chemin vers le répertoire contenant les résultats BUSCO
busco_results_dir = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/"

# Chemin vers le répertoire contenant les profondeurs de séquençage
depth_dir = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Depth/Diptera_mappings/COV/"

# Chemin vers le fichier Candidates_info.tsv
candidates_info_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Candidates_info.tsv"

# Dictionnaire pour stocker les profondeurs de séquençage
depth_dict = {}

# Lire le fichier Candidates_info.tsv pour obtenir la liste des espèces
candidates_info_df = pd.read_csv(candidates_info_file, sep="\t")
species_list = candidates_info_df['Species'].unique()

#Créer le dictionnaire qui contiendra les quantiles de distribution BUSCO
quantiles_dict = {}

# Parcourir chaque espèce
for species in species_list:
    #print(f"Traitement de l'espèce {species}")

    # Charger le fichier BUSCO correspondant à l'espèce
    busco_file = os.path.join(busco_results_dir, f"{species}/run_diptera_odb10/full_table.tsv")
    if not os.path.exists(busco_file):
        print(f"Le fichier BUSCO pour l'espèce {species} n'existe pas.")
        continue

    # Lire le fichier BUSCO en ignorant les lignes commençant par '#' et ne contenant pas "Missing"
    with open(busco_file, 'r') as f:
        lines = f.readlines()

    # Récupérer les noms des scaffolds contenant des gènes BUSCO
    unique_scaffolds = []
    for line in lines:
        if not line.startswith('#') and "Missing" not in line:
            columns = line.strip().split('\t')
            scaffold = columns[2]
            if scaffold:  # Vérifier si la valeur de la colonne est non vide
                unique_scaffolds.append(scaffold)

    # Supprimer les doublons
    unique_scaffolds = list(set(unique_scaffolds))
    
    # Vérifier si le fichier de profondeur de séquençage existe pour cette espèce
    depth_file = os.path.join(depth_dir, f"{species}.cov")
    if not os.path.exists(depth_file):
        print(f"Pas de fichier de profondeur de séquençage pour l'espèce {species}.")
        continue
    
    # Charger le fichier de profondeur de séquençage
    depth_df = pd.read_csv(depth_file, sep="\t", skiprows=1, names=['rname', 'startpos', 'endpos', 'numreads', 'covbases', 'coverage', 'meandepth', 'meanbaseq', 'meanmapq'])
    
    # Filtrer les profondeurs de séquençage pour les scaffolds contenant des gènes BUSCO
    busco_depths = depth_df[depth_df['rname'].isin(unique_scaffolds)][['rname', 'meandepth']]
    
    # Calculer les quantiles
    quantiles = {
        'quantile_5': busco_depths['meandepth'].quantile(0.05),
        'quantile_95': busco_depths['meandepth'].quantile(0.95),
        'quantile_20': busco_depths['meandepth'].quantile(0.20),
        'quantile_80': busco_depths['meandepth'].quantile(0.80)
    }
    
    # Ajouter les quantiles au dictionnaire par espèce
    quantiles_dict[species] = quantiles

    # Afficher les quantiles pour cette espèce
    print("Quantiles pour l'espèce", species, ":", quantiles)
    print("\n")

print("Traitement terminé.")


###### Partie 2 : Test d'appartenance des candidats au génome de l'insecte :

# Fonction pour déterminer si un candidat est endogène ou non
def is_endogenized(row):
    # Vérifier si la valeur de la colonne Scaffold_depth est NA
    if pd.isna(row['Scaffold_depth']):
        return "Yes"  # Endogène si la valeur est NA
    
    # Récupérer les quantiles 5%, 20%, 80% et 95% pour cette espèce
    species = row['Species']
    species_quantiles = quantiles_dict.get(species, {})
    quantile_5 = species_quantiles.get('quantile_5', None)
    quantile_20 = species_quantiles.get('quantile_20', None)
    quantile_80 = species_quantiles.get('quantile_80', None)
    quantile_95 = species_quantiles.get('quantile_95', None)
    
    # Récupérer les valeurs de la colonne Scaffold_depth et des colonnes BUSCO_count et TE_count
    scaffold_depth = row['Scaffold_depth']
    busco_count = row['BUSCO_count']
    te_count = row['TE_count']
    
    # Vérifier les conditions pour déterminer si le candidat est endogène ou non
    if busco_count > 0:
        return "Yes"  # Endogène si au moins un gène BUSCO est présent
    
    if quantile_20 is not None and quantile_80 is not None:
        if quantile_20 <= scaffold_depth <= quantile_80:
            return "Yes"  # Endogène si la profondeur est entre les quantiles 20% et 80%
    
    if quantile_5 is not None and quantile_95 is not None:
        if quantile_5 <= scaffold_depth <= quantile_95:
            if te_count > 0:
                return "Yes"  # Endogène si la profondeur est entre les quantiles 5% et 95% et TE_count est supérieur à 0
    
    # Si aucune des conditions ci-dessus n'est remplie, le candidat est exogène
    return "No"

# Ajouter une colonne "Endogenized" au DataFrame en appliquant la fonction is_endogenized à chaque ligne
candidates_info_df['Endogenized'] = candidates_info_df.apply(is_endogenized, axis=1)
print(candidates_info_df['Endogenized'])

# Sauvegarder le DataFrame mis à jour dans un nouveau fichier Candidates_info_with_endogenized.tsv
candidates_info_with_endogenized_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Candidates_info_EVEs.tsv"
candidates_info_df.to_csv(candidates_info_with_endogenized_file, sep="\t", index=False, na_rep='NA')

print("Tableau mis à jour avec la colonne 'Endogenized' ajoutée.")

