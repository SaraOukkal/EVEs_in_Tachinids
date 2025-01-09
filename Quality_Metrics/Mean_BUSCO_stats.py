import pandas as pd

# Chemin d'accès aux fichiers
busco_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/resultats.tsv'
horizon_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_horizon.txt'
dtol_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_ncbi.txt'
schizophora_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_diptera.txt'

# Charger les résultats BUSCO avec un délimiteur d'espace régulier
busco_results = pd.read_csv(busco_results_path, sep='\s+', header=0)

# Vérifier les colonnes dans busco_results
print("Colonnes du fichier BUSCO :", busco_results.columns.tolist())

# Charger les listes d'espèces pour chaque groupe
horizon_species = pd.read_csv(horizon_species_path, header=None, names=['Species'])
dtol_species = pd.read_csv(dtol_species_path, header=None, names=['Species'])
schizophora_species = pd.read_csv(schizophora_species_path, header=None, names=['Species'])

# Ajouter une colonne de groupe dans les résultats BUSCO
busco_results['Group'] = None
busco_results.loc[busco_results['Dossier'].isin(horizon_species['Species']), 'Group'] = 'Tachinidae_Horizon'
busco_results.loc[busco_results['Dossier'].isin(dtol_species['Species']), 'Group'] = 'Tachinidae_DTOL'
busco_results.loc[busco_results['Dossier'].isin(schizophora_species['Species']), 'Group'] = 'Schizophora_NCBI'

# Calculer la colonne Non_Missing
busco_results['Non_Missing'] = busco_results['Single'] + busco_results['Duplicated'] + busco_results['Fragmented']

# Calculer les moyennes pour chaque groupe
mean_results = busco_results.groupby('Group').agg(
    Mean_Single=('Single', 'mean'),
    Mean_Duplicated=('Duplicated', 'mean'),
    Mean_Fragmented=('Fragmented', 'mean'),
    Mean_Missing=('Missing', 'mean'),
    Mean_Non_Missing=('Non_Missing', 'mean')
).reset_index()

# Enregistrer les résultats moyens dans un fichier TSV
mean_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/mean_results_by_group.tsv'
mean_results.to_csv(mean_results_path, sep='\t', index=False)

print(f"Les résultats moyens ont été enregistrés sous : {mean_results_path}")

