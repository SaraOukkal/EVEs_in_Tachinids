import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import pandas as pd 
import numpy as np 

# Print out a message when the program is initiated.
#print('------------------------------------------------------------------------------------------------------\n')
#print('           ###Extract and Translate loci candidates to DNA and Amino Acide files#######.             \n')
#print('------------------------------------------------------------------------------------------------------\n')

#----------------------------------- PARSE SOME ARGUMENTS ----------------------------------------
parser = argparse.ArgumentParser(description='Extract and Translate loci candidates to DNA and Amino Acide files')
parser.add_argument("-g", "--genome_file", help="the genome assembly file")
parser.add_argument("-b", "--blast_file", help="result_mmseqs2_strand_summary_V.bed file")
parser.add_argument("-aa", "--out_aa", help="the aa output file")
parser.add_argument("-dna", "--out_dna", help="the dna output file")
parser.add_argument("-sp", "--sp_name", help="the species' name ")
args = parser.parse_args()


#Usage example 

# python3  /beegfs/data/soukkal/StageM2/Scripts/Homology_search/Extract_and_translate_loci.py -g /beegfs/data/soukkal/StageM2/Control_genomes/genomes/GCA_001014415.1.fna -b /beegfs/data/soukkal/StageM2/Control_genomes/results/Viral_Homology_results/GCA_001014415.1/run_mmseqs2/result_mmseqs2_strand_summary_V.bed -aa /beegfs/data/soukkal/StageM2/Control_genomes/results/Viral_Homology_results/GCA_001014415.1/Fasta_viral_loci.faa -dna /beegfs/data/soukkal/StageM2/Control_genomes/results/Viral_Homology_results/GCA_001014415.1/Fasta_viral_loci.fna -sp GCA_001014415.1 

#Define the variables 
Genome_assembly=args.genome_file
Blast_loci= args.blast_file
AA_Out_file=args.out_aa
DNA_Out_file=args.out_dna
species_name=args.sp_name


#Open records :
print("Opening the genome file...")
Genome_dict= SeqIO.to_dict(SeqIO.parse(Genome_assembly, "fasta"))

#Open the blast tab file : 
print("Openin the blast table...")
Blast_tab = pd.read_csv(Blast_loci,sep=" ")
#Add plus strand coordinates 

#m = Blast_tab['strand'].eq('-')
#Blast_tab['Newqstart'] = np.where(m, Blast_tab['qlen'].sub(Blast_tab['end']), Blast_tab['start'])
#Blast_tab['Newqend'] = np.where(m, Blast_tab['qlen'].sub(Blast_tab['start']), Blast_tab['end'])

# Genome_assembly="/beegfs/data/soukkal/StageM2/Control_genomes/genomes/GCA_001014415.1.fna"
# Blast_loci="/beegfs/data/soukkal/StageM2/Control_genomes/results/Viral_Homology_results/GCA_001014415.1/run_mmseqs2/result_mmseqs2_strand_summary_V.bed"

#Write new fasta loci aa file 
with open(AA_Out_file,"w") as output_aa:
    for index, row in Blast_tab.iterrows():
        seq_name= row['query']+":"+str(row['start'])+'-'+str(row['end'])+"("+row['strand']+")"+":"+species_name
        if row['strand'] =="+":
            print('>',seq_name,sep="",file=output_aa)
            print(str(Genome_dict[row['query']].seq[row['start']-1:row['end']].translate()),file=output_aa) # since python uses 0-based 
        elif row['strand'] =="-":
            print('>',seq_name,sep="",file=output_aa)
            print((Genome_dict[row['query']][row['start']-1: row['end']].seq.reverse_complement().translate()),file=output_aa)

#Write new fasta loci dna file 
with open(DNA_Out_file,"w") as output_dna:
    for index, row in Blast_tab.iterrows():
        seq_name= row['query']+":"+str(row['start'])+'-'+str(row['end'])+"("+row['strand']+")"+":"+species_name
        if row['strand'] =="+":
            print('>',seq_name,sep="",file=output_dna)
            print(str(Genome_dict[row['query']].seq[row['start']-1:row['end']]),file=output_dna) # since python uses 0-based 
        elif row['strand'] =="-":
            print('>',seq_name,sep="",file=output_dna)
            print((Genome_dict[row['query']][row['start']-1: row['end']].seq.reverse_complement()),file=output_dna)


print("The translated file has been saved to :", AA_Out_file)
print("Nucleotide file has been saved to :", DNA_Out_file)

