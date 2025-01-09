##Récupérer les séquences de Schizophora : 

from ete3 import NCBITaxa

# Charger la base de données taxonomique de NCBI
ncbi = NCBITaxa()

# TaxID de Schizophora
schizophora_taxid = 43738

# Récupérer tous les descendants de Schizophora
schizophora_descendants = ncbi.get_descendant_taxa(schizophora_taxid)

# Créer un dictionnaire pour associer les numéros d'accession aux taxID
accession_to_taxid = {}

# Ouvrir le fichier qui contient les numéros d'accession et les Taxid et remplir le dictionnaire
with open("/beegfs/data/soukkal/Thesis/Databases/CreateTaxDB/Accession_Taxid.tsv", "r") as access_taxid_file:
    for line in access_taxid_file:
        accession, taxid = line.strip().split("\t")
        accession_to_taxid[accession] = taxid

# Ouvrir le fichier FASTA d'entrée et créer un fichier de sortie
with open("/beegfs/data/soukkal/Thesis/Databases/Refseq.faa", "r") as input_file, open("/beegfs/data/soukkal/Thesis/Databases/Refseq_Schizophora.faa", "w") as output_file:
    write_sequence = False  # Nous n'écrirons pas de séquence tant que nous n'aurons pas vérifié le taxid
    for line in input_file:
        if line.startswith(">"):
            accession = line[1:].split()[0]
            taxid = accession_to_taxid.get(accession)
            if taxid and int(taxid) in schizophora_descendants: #Int(taxid) parce que les taxid dans les headers sont des chaines de caractères alors que dans le tableau c'est des chiffres
                write_sequence = True  # Commencez à écrire la séquence
                output_file.write(line)  # Écrire l'en-tête
            else:
                write_sequence = False  # Ne pas écrire la séquence
        elif write_sequence:
            output_file.write(line)  # Écrire la séquence
