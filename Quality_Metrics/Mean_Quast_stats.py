import pandas as pd

# Chemins d'accès aux fichiers
quast_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/Quast/Assembly_stats.txt'
horizon_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_horizon.txt'
dtol_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_ncbi.txt'
schizophora_species_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_diptera.txt'

# Charger les résultats QUAST
quast_results = pd.read_csv(quast_results_path, sep='\s+', header=0)

# Vérifier les colonnes dans quast_results
print("Colonnes du fichier QUAST :", quast_results.columns.tolist())

# Charger les listes d'espèces pour chaque groupe
horizon_species = pd.read_csv(horizon_species_path, header=None, names=['Species'])
dtol_species = pd.read_csv(dtol_species_path, header=None, names=['Species'])
schizophora_species = pd.read_csv(schizophora_species_path, header=None, names=['Species'])

# Ajouter une colonne de groupe dans les résultats QUAST
quast_results['Group'] = None
quast_results.loc[quast_results['species'].isin(horizon_species['Species']), 'Group'] = 'Tachinidae_Horizon'
quast_results.loc[quast_results['species'].isin(dtol_species['Species']), 'Group'] = 'Tachinidae_DTOL'
quast_results.loc[quast_results['species'].isin(schizophora_species['Species']), 'Group'] = 'Schizophora_NCBI'

# Convertir les valeurs en Mb pour Total_Length
quast_results['Total_Length_Mb'] = quast_results['Total_Length'] / 1_000_000

# Calculer les intervalles et les moyennes
stats_results = quast_results.groupby('Group').agg(
    Length_Min=('Total_Length_Mb', 'min'),
    Length_Max=('Total_Length_Mb', 'max'),
    Mean_Length=('Total_Length_Mb', 'mean'),
    GC_Min=('GC', 'min'),
    GC_Max=('GC', 'max'),
    Mean_GC=('GC', 'mean'),
    N50_Min=('N50', 'min'),
    N50_Max=('N50', 'max'),
    Mean_N50=('N50', 'mean')
).reset_index()

# Enregistrer les résultats statistiques dans un fichier TSV
stats_results_path = '/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/Quast/stats_results_by_group.tsv'
stats_results.to_csv(stats_results_path, sep='\t', index=False)

print(f"Les résultats statistiques ont été enregistrés sous : {stats_results_path}")

