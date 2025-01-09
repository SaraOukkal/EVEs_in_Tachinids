#!/bin/bash

#Ouvrir une liste de noms d’espèces :
species_list=`cat /beegfs/data/soukkal/Thesis/Tachinid_Project/All_species_list.txt`

#On se place dans le dossier qui contient les résultats
cd /beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/Quast/Quast_sum

#On crée le fichier final et on y insère un header :
echo species Total_Length GC N50 > ../Assembly_stats.txt

#On crée une boucle qui va parcourir la liste de noms et appliquer les commandes précédentes à chaque fichier :
for i in $species_list
do
        echo $i > tmp_file.txt
        grep "Total length\|GC\|N50" report_$i.tsv | tail -n 3 | cut -f2 |paste -sd ' ' >> tmp_file.txt
        paste -sd ' '  tmp_file.txt >> ../Assembly_stats.txt
done

rm tmp_file.txt
