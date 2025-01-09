#!/bin/bash
#SBATCH -J mmseqs_taxdb   # Nom de la tâche
#SBATCH -o /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Homology_search/Logs/mmseqs_taxdb.out  # Fichier de sortie standard
#SBATCH -e /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Homology_search/Logs/mmseqs_taxdb.err  # Fichier de sortie d'erreur
#SBATCH --cpus-per-task=8    # Nombre de cœurs CPU par tâche
#SBATCH --mem=8G             # Mémoire par tâche
#SBATCH --time=12:00:00      # Temps maximal d'exécution (hh:mm:ss)

mmseqs="/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs"

##Créer une BDD de mmseqs :  
$mmseqs createdb /beegfs/data/soukkal/Thesis/Databases/Refseq.faa /beegfs/data/soukkal/Thesis/Databases/Refseq
##Créer la base de données taxonomique (utiliser le tableau fait juste avant pour l'option --tax-mapping-file): 
$mmseqs createtaxdb /beegfs/data/soukkal/Thesis/Databases/Refseq tmp --ncbi-tax-dump /beegfs/data/soukkal/Thesis/Databases/CreateTaxDB/ncbi-taxdump --tax-mapping-file /beegfs/data/soukkal/Thesis/Databases/CreateTaxDB/Accession_Taxid.tsv
#Retirer les schizophora et apocrita : 
$mmseqs filtertaxseqdb /beegfs/data/soukkal/Thesis/Databases/Refseq /beegfs/data/soukkal/Thesis/Databases/Refseq_noschizophora_noapocrita --taxon-list '!43738&&!7400'
#Refaire un createtaxdb sur la BDD filtrée : 
$mmseqs createtaxdb /beegfs/data/soukkal/Thesis/Databases/Refseq_noschizophora_noapocrita tmp
