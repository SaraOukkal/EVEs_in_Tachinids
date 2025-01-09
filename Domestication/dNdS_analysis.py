#!/usr/bin/python
import sys 
import pandas as pd
#First take the ete3 package modified to get the Muse and Gant model:
#sys.path.insert(0, "/beegfs/home/soukkal/miniconda3/lib/python3.9/site-packages/ete3/")
sys.path.insert(0, "/beegfs/home/bguinet/.local/lib/python3.11/site-packages/")
from ete3 import EvolTree
from ete3.treeview.layouts import evol_clean_layout
from ete3 import faces
#import pandas as pd
import numpy as np
import sys
import argparse
import re
import os
import subprocess
from Bio import SeqIO
from Bio.Align.Applications import MuscleCommandline
from Bio import Phylo
from Bio.Seq import Seq
import statistics
from io import StringIO


# Arguments provenant de Snakemake via argparse
parser = argparse.ArgumentParser(description='dNdS analysis for cluster of genes.')
parser.add_argument('--tree', type=str, required=True, help='Path to the tree file')
parser.add_argument('--aln', type=str, required=True, help='Path to the alignment file')
parser.add_argument('--output_dir', type=str, required=True, help='Output directory')
parser.add_argument('--cluster_name', type=str, required=True, help='Cluster name')
parser.add_argument('--EVEs_list', type=str, required=True, help='Path to the EVEs list file')
parser.add_argument('--event_number', type=str, required=True, help='Event number')

args = parser.parse_args()

tree_file = args.tree
alignment_file = args.aln
output_dir = args.output_dir
cluster_name = args.cluster_name
EVEs_list_path = args.EVEs_list
event_number = args.event_number

print('------------------------------------------------------------------------------------------------------------------------------\n')
print('#dNdS_analyzer.\n')
print('\n')
print('This script will create bash scripts in order to run dNdS analysis for cluster of genes.\n')
print('Each analysis will fill a dataframe with \n')
print('\n')
print(' Cluster_name  Mean_dNdS   Pvalue_dNdS\n')
print(' cluster1      0.34        0.20           <------- foreground branches are not significantly different from 1.    \n') 
print(' cluster2      0.002       0.000003       <------- foreground branches are significantly different from 1.\n')
print('\n')
print('foreground branches being the hymenoptera candidate loci for a viral domestication. \n')
print('If foreground branches are significantly different from 1, then the dN/dS values are under purifying or adaptative selection. \n')
print('------------------------------------------------------------------------------------------------------------------------------\n')

# Charger la liste des EVEs
with open(EVEs_list_path) as f:
    list_of_EVEs_to_test = [item for line in f for item in line.strip().split()]

# if output_dir does not exists, create it 
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

print("list_of_EVEs_to_test;", list_of_EVEs_to_test)

print("----------------------------------------------------------------")
print("Beginning of the dNdS analysis for the cluster :" + cluster_name)
print("----------------------------------------------------------------")
print("\n")

print("- The tree and the alignment codon file of the cluster of interest have been charged -")
#tree = EvolTree(tree_file, binpath="/beegfs/data/soukkal/TOOLS/paml4.9j/bin")
#tree = EvolTree(tree_file, binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")
#print(tree)
# We need to root it otherwise the previous clades are not respected.


print("Opening tree file to root it.")
tree = Phylo.read(tree_file, 'newick')
print(tree)
tree.root_at_midpoint()
output_tree_file = tree_file+ '.rooted'
Phylo.write(tree, output_tree_file, 'newick')
print("Rooted tree saved to : ", output_tree_file)

#R = tree.get_midpoint_outgroup()
#tree.set_outgroup(R)

#rooted_tree_path = tree_file+".Rooted"
#from Bio import Phylo
#Phylo.write(tree, rooted_tree_path, "newick")

tree = EvolTree(output_tree_file, binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")

#tree = EvolTree(tree_file, binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")

print(tree)

List_leaf_node = []
for leaf in tree:
    List_leaf_node.append(leaf.name)
List_node_to_keep = []
for record in SeqIO.parse(alignment_file, "fasta"):
    List_node_to_keep.append(record.id)

if len(List_leaf_node) != len(List_node_to_keep):
    print("TREE PRUNED")
    tree.prune(List_node_to_keep)
    tree.write(format=1, outfile=alignment_file + ".Pruned")
    #tree = EvolTree(alignment_file + ".Pruned", binpath="/beegfs/data/soukkal/TOOLS/paml4.9j/bin")
    tree = EvolTree(alignment_file + ".Pruned", binpath="/beegfs/data/bguinet/TOOLS/paml4.9i/bin")

print("\n")
print('Tree of the cluster : ', cluster_name)
print(tree)
print("\n")

# Modify RecordID to fit tree ID 

output_file = alignment_file + ".newID"
with open(alignment_file, "r") as input_handle, open(output_file, "w") as output_handle:
    for record in SeqIO.parse(input_handle, "fasta"):
        # Replace "(-):" with "_-__", "(+):" with "_+__", and ":" with "_"
        record.id = record.id.replace("(-):", "_-__").replace("(+):", "____").replace(":", "_")
        # Save the modified records to a new file
        record.description = ""
        SeqIO.write(record, output_handle, "fasta")

print("Alignment with new RecordIDs saved to : ", output_file)

df_marks = pd.DataFrame(columns=['Node_name', 'Node id'])
list_of_EVEs_to_test_newID = [record.replace("(-):", "_-__").replace("(+):", "____").replace(":", "_")
 for record in list_of_EVEs_to_test]

ancestor=tree.get_common_ancestor(list_of_EVEs_to_test_newID)
print(list_of_EVEs_to_test)
print(ancestor)
print("Type4")
print(type(ancestor))

for node in ancestor.traverse("levelorder"):	
	try:
		df_marks=df_marks.append({"Node_name":node.name,"Node_id":node.node_id}, ignore_index=True)
	except:
		continue

print(df_marks)
print(df_marks.shape[0]) 

if df_marks.shape[0] < 5:
	df_marks=df_marks.loc[df_marks["Node_name"].str.contains("__")]
else:
	df_marks = df_marks.iloc[1:]

tree.workdir = output_dir
tree.link_to_alignment(output_file)

print("\n")
print("-------------------------------------------------------------------------")
print(".  Adding of node marks...                                               ")
print(".  the marks will allow to force the loci to be set to a dN/dS value =1  ")
print("-------------------------------------------------------------------------")

# mark a group of branches of interest

marks=list(df_marks['Node_id'].astype(int))
print("Marks : ", marks)
tree.mark_tree(marks, ['#1'] * len(marks))

print (tree.write ())


print("--------------------------------------------------")
print("      Running of the the neutral model         ")
print(" all marked sequence are forced to have a dNdS =1 ")
print("--------------------------------------------------")
import shutil
#remove previous run
if os.path.isdir(output_dir +"/"+ '/b_neut'+'.'+str(cluster_name)+'_'+str(event_number)):
 shutil.rmtree(output_dir +"/"+ 'b_neut'+'.'+str(cluster_name)+'_'+str(event_number)+'/')

if os.path.isdir(output_dir +"/"+ 'b_free'+'.'+str(cluster_name)+'_'+str(event_number)):
 shutil.rmtree(output_dir +"/"+ 'b_free'+'.'+str(cluster_name)+'_'+str(event_number)+'/')

print("\n")
print("--------------------------------------")
print("Running of the the dNdS mode")
print("--------------------------------------")
print("\n")
tree.run_model('b_free' + '.' + str(cluster_name) + '_' + str(event_number), CodonFreq=4, estFreq=1,getSE=1,noisy=3,verbose=1)
print("Model run done")

# Load the model
#tree.link_to_evol_model(output_dir +"/"+ 'fb' + '.' + str(cluster_name) + '_' + str(event_number) + '/out', 'fb')
#for model in tree._models:
#    model = tree.get_evol_model(model)
#print(model)

# Here we will get the model output in a dataframe:
#df = pd.read_csv(StringIO(model.__str__()), skiprows=6, names=['marks', 'omega', 'node_ids', 'name'])
#df = df.applymap(lambda x: x.split(":")[1])
# Remove white space
#df = df.applymap(str.strip).rename(columns=str.strip)
#print(df)
#marks = []
#omega = []
#for i in list_of_EVEs_to_test:
#    print(i)
#    s = df.loc[df['name'] == i, 'node_ids']
#    print("s:", s)
    # dNdS = df.loc[df['name'] == i,'omega']
    # print(dNdS)
    # dNdS = float(dNdS)
    # s=int(s)
    # print("dNdS of",i," equal: ", dNdS)
    # print(dNdS)
    # omega.append(dNdS)
#    marks.append(s)
#print("Node marked :")
#for i in list_of_EVEs_to_test:
#    print("- ", i)
# dNdS_average=statistics.mean(omega)
# print("dN/dS average is :",statistics.mean(omega))

#print("\n")
#print("-------------------------------------------------------------------------")
#print(".  Adding of node marks...                                               ")
#print(".  the marks will allow to force the loci to be set to a dN/dS value =1  ")
#print("-------------------------------------------------------------------------")
# mark a group of branches of interest
#tree.mark_tree(marks, ['#1'])
#print("\n")

print("--------------------------------------------------")
print("      Running of the the neutral model         ")
print(" all marked sequence are forced to have a dNdS =1 ")
print("--------------------------------------------------")


tree.run_model('b_neut' + '.' + str(cluster_name) + '_' + str(event_number), CodonFreq=4, estFreq=1, getSE=1,noisy=3,verbose=1)
print("neutral model analysis done.")
print("\n")


print("--------------------------------------------------")
print("      Running of the the free-ratio model            ")
print(" all marked sequence are allowed to be free       ")
print("--------------------------------------------------")
# tree.run_model('b_free'+'.'+str(cluster_name)+'_'+str(event_number),CodonFreq=4,estFreq=1)
print("free-ratio model analysis done.")
print("\n")
tree.link_to_evol_model(output_dir + "/" + 'b_neut' + '.' + str(cluster_name) + '_' + str(event_number) + '/out', 'b_neut')
for model in tree._models:
    fb_model_neut = tree.get_evol_model(model)
# fb_model_neut = tree.get_evol_model('b_neut'+"."+cluster_name+"_"+event_number)
# fb_model_free = tree.get_evol_model("b_free"+"."+cluster_name+"_"+str(event_number))
tree.link_to_evol_model(output_dir + "/" + 'b_free' + '.' + str(cluster_name) + '_' + str(event_number) + '/out', 'b_free')
for model in tree._models:
    fb_model_free = tree.get_evol_model(model)

# Charging the free model output
df_free = pd.read_csv(StringIO(fb_model_free.__str__()), skiprows=6, names=['marks', 'omega', 'node_ids', 'name'])
df_free = df_free.applymap(lambda x: x.split(":")[1])
# Remove white space
df_free = df_free.applymap(str.strip).rename(columns=str.strip)
print(df_free)

print("\n")
omega = []
for i in list_of_EVEs_to_test_newID:
    s = df_free.loc[df_free['name'] == i, 'node_ids']
    # print(s)
    dNdS = df_free.loc[df_free['name'] == i, 'omega']
    print("dNdS : ", dNdS)
    dNdS = float(dNdS)
    s = int(s)
    # print("dNdS of",i," equal: ", dNdS)
    # print(dNdS)
    omega.append(dNdS)

print("Free model output: ")
print(df_free)
print("\n")
print("dN/dS average is :", omega[0])  # Take the first omega value
print("\n")




print("Free model output: ")
print(df_free)
print("\n")
print("dN/dS average is :",omega[0]) #Take the first omega value
print("\n")

#Run of the likelihhod ratio test...
print("------------------------------------------------------------------")
print("           Running of the likelihood ratio test                   ")
print("                                                                  ")
print("         if pvalue 'b_free' vs 'b_neut' < 0.05:                   ")
print(" the foreground branches are significantly different from 1 '     ")
print("------------------------------------------------------------------")

#pvalue=tree.get_most_likely ('b_free'+"."+Cluster_name+"_"+str(event_number), 'b_neut'+"."+Cluster_name+"_"+str(event_number))	
try:
    		from scipy.stats import chi2
    		chi_high = lambda x, y: 1 - chi2.cdf(x, y)
except ImportError:
    		from .utils import chi_high
def get_most_likely2(altn, null):
 		if hasattr(altn, 'lnL') and hasattr(null, 'lnL'):
  			if  null.lnL - altn.lnL < 0:
    				return chi_high(2 * abs(altn.lnL - null.lnL),float(altn.np - null.np))
  			else:
  				warn("\nWARNING: Likelihood of the alternative model is ""smaller than null's (%f - %f = %f)" % (null.lnL, altn.lnL, null.lnL - altn.lnL) + "\nLarge differences (> 0.1) may indicate mistaken ""assigantion of null and alternative models")
  				return 1
pvalue=get_most_likely2(fb_model_free,fb_model_neut)
print("pvalue estimated :", pvalue)
print("\n")

#get the Standar Error of the dN/dS estimation : 
with open(output_dir +"/"+'b_free'+'.'+str(cluster_name)+'_'+str(event_number)+'/out', "r") as ifile:
		for line in ifile:
			if line.startswith("SEs for parameters:"):
				SE=next(ifile, ' ').strip()
				SE=re.split('\s+', SE)
				SE=SE[-1]

#Insert the value into a new dataframe 
df = pd.DataFrame(columns=("Clustername","Event","Mean_dNdS","Pvalue_dNdS","SE_dNdS"))
df=df.append({"Clustername":cluster_name,"Event":str(event_number),"Mean_dNdS":omega[0],"Pvalue_dNdS":pvalue,"SE_dNdS":SE},ignore_index=True)
print("Global result :")
print(df)
df.to_csv(output_dir+ "/"+"dNds_"+cluster_name+"_"+str(event_number)+".out",sep="\t",index=False)


"""
# Run of the likelihhod ratio test...
print("------------------------------------------------------------------")
print("           Running of the likelihood ratio test                   ")
print("                                                                  ")
print("         if pvalue 'b_free' vs 'b_neut' < 0.05:                   ")
print(" the foreground branches are significantly different from 1 '     ")
print("------------------------------------------------------------------")

# pvalue=tree.get_most_likely ('b_free'+"."+cluster_name+"_"+str(event_number), 'b_neut'+"."+cluster_name+"_"+str(event_number))

try:
	from scipy.stats import chi2
	chi_high = lambda x: chi2.sf(x, 1)
	LRT = 2 * (fb_model_free.LnL - fb_model_neut.LnL)
	print(LRT)
	pvalue = chi_high(LRT)
	print("pvalue :", pvalue)
except ImportError:
	print("Module scipy.stats not found")
	pvalue = "NA"

print("\n")
print("-------------------------------------------------------------------")
print(".              Filling of the dataframe with output:              ")
print("-------------------------------------------------------------------")
print("\n")

data = {'Cluster_name': [cluster_name], 'Mean_dNdS': [omega[0]], 'Pvalue_dNdS': [pvalue]}
df = pd.DataFrame(data)
df.to_csv(output_dir + "dNdS_output_" + cluster_name + "_" + event_number + ".csv", index=False)

#import csv

## Define the data
#data = {'Cluster_name': [cluster_name], 'Mean_dNdS': [omega[0]], 'Pvalue_dNdS': [pvalue]}
# Create the output file path
#output_file = output_dir + "dNdS_output_" + cluster_name + "_" + event_number + ".csv"

# Open the file and write the data
#with open(output_file, mode='w', newline='') as file:
#    writer = csv.writer(file)
#    # Write the header
#    writer.writerow(['Cluster_name', 'Mean_dNdS', 'Pvalue_dNdS'])
#    # Write the data row
#    writer.writerow([cluster_name, omega[0], pvalue])

print("Analysis completed!")
"""
