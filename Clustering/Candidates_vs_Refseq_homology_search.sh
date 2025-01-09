#!/bin/bash
#SBATCH -J Mmseqs_search
#SBATCH --partition=normal
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH -o /beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Clustering/mmseqs_candidates_vs_virus.out
#SBATCH -e /beegfs/data/soukkal/Thesis/Tachinid_Project/Scripts/Clustering/mmseqs_candidates_vs_virus.error
#SBATCH --constraint="skylake|haswell|broadwell"
#SBATCH --exclude=pbil-deb27

Output_dir="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/mmseqs2_candidates_Refseq_virus"
Candidates="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidates_in_clusters.faa"
Candidates_db="/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidates_in_clusters"
Viral_db="/beegfs/data/soukkal/Thesis/Databases/Refseq_viral_db_nophages_noPDV_dEVE_FV"


#Create a directory where the BlastX results will be written
mkdir -p $Output_dir
#Create the query database (for each query genome/sequence)
/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs createdb $Candidates $Candidates_db
#Run Mmseqs2 search (BlastX equivalent, homology search between query and db)
/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs search $Candidates_db $Viral_db $Output_dir/result_mmseqs2 $Output_dir/tpm -a --start-sens 1 --sens-steps 3 -s 7 --max-accept 1 --threads 4 --remove-tmp-files
#From the previous step you get one result file per thread, the next step will format the result column and also put together all the results in one file:
/beegfs/data/soukkal/TOOLS/mmseqs/bin/mmseqs convertalis --format-output 'query,qlen,tlen,target,pident,qstart,qend,tstart,tend,evalue,tcov' $Candidates_db $Viral_db $Output_dir/result_mmseqs2 $Output_dir/result_mmseqs2.m8