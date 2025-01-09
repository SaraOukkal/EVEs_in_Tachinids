#Récup la liste de tous les gènes BUSCO complets pour chacune des espèces
for file in /beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/*/run_diptera_odb10/full_table.tsv
    do
    grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids.txt;
    done

#Retirer ceux qui sont présents chez moins de 75/80 espèces :
##Compter les occurences de chaque gène : 
sort /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids.txt |uniq -c > /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids_with_counts.txt
##Retirer les espaces des débuts de ligne que génère uniq -c
sed 's/^[ \t]*//' /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids_with_counts.txt > /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids_with_counts_without_spaces.txt
awk '{if($1 > 77) print $2}' /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/complete_busco_ids_with_counts_without_spaces.txt > /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/total_ok_busco_ids.txt
#Choisir 300 gènes parmi les gènes restants :
shuf -n 300 /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/total_ok_busco_ids.txt > /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/final_busco_ids.txt


#Copier les fichiers fasta des nt dans un dossier
#Renommer les headers des séquences avec le nom de l'espèce et le nom du gène
#Faire un fichier par Gène par espèce
mkdir -p /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt
species_list=`cat /beegfs/data/soukkal/Thesis/Tachinid_Project/All_species_list.txt`
Genes_list=`cat /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/final_busco_ids.txt`

for sp in $species_list
    do
    for Gene in $Genes_list
        do
        if test -f "/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/$sp/run_diptera_odb10/busco_sequences/single_copy_busco_sequences/$Gene.fna"
            then
            cp /beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/$sp/run_diptera_odb10/busco_sequences/single_copy_busco_sequences/$Gene.fna /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/$sp.$Gene.fna
            sed -i "/>/c\>$sp" /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/$sp.$Gene.fna
            fi
        done
    done

#Regrouper les fichier pour avoir un fichier par gène BUSCO (mettre dans un sous dossier)
mkdir -p /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/Genes_files/

for sp in $species_list
    do
    for Gene in $Genes_list
        do
        if test -f "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/$sp.$Gene.fna"
            then
            cat /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/$sp.$Gene.fna >> /beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Species_Phylogeny/busco_nt/Genes_files/$Gene.fasta
            fi
        done
    done
