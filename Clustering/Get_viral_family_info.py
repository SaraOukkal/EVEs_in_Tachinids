from ete3 import NCBITaxa
import csv

ncbi = NCBITaxa()

def get_viral_family(virus_species):
    try:
        name2taxid = ncbi.get_name_translator([virus_species])
        taxid = name2taxid[virus_species][0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        for taxid in lineage:
            if ranks[taxid] == 'family':
                return names[taxid]
    except Exception as e:
        # Si la famille n'est pas trouvée avec ETE3, chercher dans le fichier manuel
        with open("/beegfs/data/soukkal/Thesis/Databases/Viral_families_added_manually.tsv", mode='r') as manual_file:
            manual_reader = csv.reader(manual_file, delimiter='\t')
            for row in manual_reader:
                if len(row) == 2 and row[0] == virus_species:
                    return row[1]
        print(f"Pas de famille trouvée pour {virus_species}: {e}")
    
    return "NA"

# Chemin vers le fichier CSV des structures génomiques
structure_path = "/beegfs/data/soukkal/Thesis/Databases/Viral_families_structure.csv"

# Lire le fichier CSV et créer un dictionnaire pour les familles et leur structure génomique
family_to_structure = {}
with open(structure_path, mode='r') as infile:
    reader = csv.DictReader(infile, delimiter=';', fieldnames=['Col1', 'Col2', 'Col3', 'Col4', 'Col5', 'Col6', 'Col7', 'Col8', 'Col9', 'Col10', 'Famille_virale', 'Structure_genomique'])
    for row in reader:
        family_name = row.get('Famille_virale', 'NA')
        genomic_structure = row.get('Structure_genomique', 'NA')
        family_to_structure[family_name] = genomic_structure  # Ajouter la famille et sa structure génomique correspondante

# Chemin vers votre fichier TSV
tsv_path = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/mmseqs2_candidates_Refseq_virus/Viral_protein_name_species.tsv"

output_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidates_summary_table/Viral_information.tsv"
with open(output_file, 'w') as output:
    with open(tsv_path, 'r') as file:
        tsv_reader = csv.reader(file, delimiter='\t')
        for row in tsv_reader:
            if len(row) == 3:  # Assurez-vous que la ligne a exactement 3 colonnes
                num_accession, protein_name, species_name = row
                family = get_viral_family(species_name)
                genomic_structure = family_to_structure.get(family, "NA")
                
                # Si la structure génomique n'est pas trouvée dans le fichier, rechercher manuellement
                if genomic_structure == "NA":
                    with open("/beegfs/data/soukkal/Thesis/Databases/Viral_structure_added_manually.tsv", mode='r') as manual_structure_file:
                        manual_structure_reader = csv.reader(manual_structure_file, delimiter='\t')
                        for row in manual_structure_reader:
                            if len(row) == 2 and row[0] == family:
                                genomic_structure = row[1]
                                break
                    if genomic_structure == "NA":
                        print(f"Pas de structure génomique trouvée pour {family}")
                
                output.write(f"{num_accession}\t{protein_name}\t{species_name}\t{family}\t{genomic_structure}\n")
