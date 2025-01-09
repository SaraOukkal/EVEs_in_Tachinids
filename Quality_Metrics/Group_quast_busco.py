import pandas as pd

# Chemins d'accès aux fichiers
quast_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/Quast/Assembly_stats.txt'
busco_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/resultats.tsv'
horizon_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_horizon.txt'
dtol_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_ncbi.txt'
schizophora_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_diptera.txt'

# Charger les résultats QUAST
quast_results = pd.read_csv(quast_results_path, sep='\s+', header=0)

# Charger les résultats BUSCO avec un délimiteur d'espace régulier
busco_results = pd.read_csv(busco_results_path, sep='\s+', header=0)

# Vérifier les colonnes dans quast_results
print("Colonnes du fichier QUAST :", quast_results.columns.tolist())

# Vérifier les colonnes dans busco_results
print("Colonnes du fichier BUSCO :", busco_results.columns.tolist())

# Charger les listes d'espèces pour chaque groupe
horizon_species = pd.read_csv(horizon_species_path, header=None, names=['Species'])
dtol_species = pd.read_csv(dtol_species_path, header=None, names=['Species'])
schizophora_species = pd.read_csv(schizophora_species_path, header=None, names=['Species'])

# Ajouter une colonne de groupe dans les résultats QUAST
quast_results['Group'] = None
quast_results.loc[quast_results['species'].isin(horizon_species['Species']), 'Group'] = 'Tachinidae_Horizon'
quast_results.loc[quast_results['species'].isin(dtol_species['Species']), 'Group'] = 'Tachinidae_DTOL'
quast_results.loc[quast_results['species'].isin(schizophora_species['Species']), 'Group'] = 'Schizophora_NCBI'

# Convertir les valeurs de longueur en Mb
quast_results['Total_Length_Mb'] = quast_results['Total_Length'] / 1_000_000

# Ajouter une colonne de groupe dans les résultats BUSCO
busco_results['Group'] = None
busco_results.loc[busco_results['Dossier'].isin(horizon_species['Species']), 'Group'] = 'Tachinidae_Horizon'
busco_results.loc[busco_results['Dossier'].isin(dtol_species['Species']), 'Group'] = 'Tachinidae_DTOL'
busco_results.loc[busco_results['Dossier'].isin(schizophora_species['Species']), 'Group'] = 'Schizophora_NCBI'

# Calculer la colonne Non_Missing
busco_results['Non_Missing'] = busco_results['Single'] + busco_results['Duplicated'] + busco_results['Fragmented']

# Fusionner les deux DataFrames sur la colonne Group
merged_results = pd.merge(
    quast_results[['species', 'Total_Length_Mb', 'Group']],
    busco_results[['Dossier', 'Single', 'Duplicated', 'Fragmented', 'Missing', 'Group']],
    left_on='Group',
    right_on='Group',
    how='inner'
)

# Renommer les colonnes pour la lisibilité
merged_results.rename(columns={'species': 'Species', 'Total_Length_Mb': 'Assembly_Size(Mb)', 'Dossier': 'Dossier'}, inplace=True)

# Sélectionner les colonnes pertinentes
final_results = merged_results[['Species', 'Group', 'Assembly_Size(Mb)', 'Single', 'Duplicated', 'Fragmented', 'Missing']]

# Trier par taille d'assemblage
final_results.sort_values(by='Assembly_Size(Mb)', inplace=True)

# Enregistrer les résultats fusionnés dans un fichier TSV
output_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/Merged_busco_quast_results.tsv'
final_results.to_csv(output_results_path, sep='\t', index=False)

print(f"Les résultats fusionnés ont été enregistrés sous : {output_results_path}")

