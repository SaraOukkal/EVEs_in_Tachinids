import os
import pandas as pd

import re

import re

def Get_candidates_info(candidates_file):
    candidates_info = {}  # Dictionnaire pour stocker les informations des candidats par espèce
    with open(candidates_file, 'r') as infile:
        for line in infile:
            columns = line.strip().split('\t')
            candidat = columns[0]
            viral_protein_accession = columns[3]
            evalue = float(columns[9])
            candidate_info = {'viral_protein_accession': viral_protein_accession, 'evalue': evalue}
            
            # Utilisation d'une expression régulière pour extraire les informations
            match = re.match(r'([^:]+):(\d+)-(\d+)(\([-+]\)):(.+)$', candidat)
            if match:
                scaffold = match.group(1)
                start = int(match.group(2))
                end = int(match.group(3))
                strand = match.group(4)[1]  # Retirer les parenthèses
                species = match.group(5)
            else:
                # Afficher un avertissement si aucune correspondance n'est trouvée
                print("Avertissement : Aucune correspondance trouvée pour le candidat:", candidat)
                continue
            
            species_info = {'species': species, 'scaffold': scaffold, 'start': start, 'end': end, 'strand': strand}
            
            # Mettre à jour le dictionnaire avec les informations extraites
            if species not in candidates_info:
                candidates_info[species] = {}
            candidates_info[species][candidat] = {**species_info, **candidate_info}
    return candidates_info


def Get_virus_info(virus_info_file, candidates_info):
    for species, species_data in candidates_info.items():
        with open(virus_info_file, 'r') as infile:
            for line in infile:
                accession, protein_name, virus_species, virus_family, virus_genomic_structure = line.strip().split('\t')
                for candidat, candidat_data in species_data.items():
                    if candidat_data['viral_protein_accession'] == accession:
                        candidat_data.update({'viral_protein_name': protein_name,
                                              'virus_species': virus_species,
                                              'virus_family': virus_family,
                                              'virus_genomic_structure': virus_genomic_structure})
                        break

def Get_BUSCO_count(candidates_info, busco_folder):
    for species, species_data in candidates_info.items():
        busco_file = os.path.join(busco_folder, f"{species}_scaffold_BUSCO_count.tsv")
        with open(busco_file, 'r') as infile:
            next(infile)  # Ignorer le header
            busco_counts = {}  
            for line in infile:
                scaffold, busco_count = line.strip().split('\t')
                busco_counts[scaffold] = int(busco_count)

            for candidat, candidat_data in species_data.items():
                scaffold = candidat_data['scaffold']
                busco_count = busco_counts.get(scaffold, 0)
                candidat_data['BUSCO_count'] = busco_count
                if scaffold not in busco_counts:
                    candidat_data['BUSCO_count'] = 0

def Get_TE_count(candidates_info, te_folder):
    for species, species_data in candidates_info.items():
        te_file = os.path.join(te_folder, f"{species}_scaffold_count.tsv")
        with open(te_file, 'r') as infile:
            next(infile)  # Ignorer le header
            te_counts = {}  
            for line in infile:
                scaffold, te_count = line.strip().split('\t')
                te_counts[scaffold] = int(te_count)

            for candidat, candidat_data in species_data.items():
                scaffold = candidat_data['scaffold']
                te_count = te_counts.get(scaffold, 0)
                candidat_data['TE_count'] = te_count
                if scaffold not in te_counts:
                    candidat_data['TE_count'] = 0

def add_scaffold_depth(candidates_info, depth_folder):
    for species, species_data in candidates_info.items():
        depth_file = os.path.join(depth_folder, f"{species}.cov")
        
        # Vérifier si le fichier de l'espèce existe
        if not os.path.exists(depth_file):
            print(f"Fichier de profondeur introuvable pour l'espèce {species}. Ignoré.")
            continue
        
        scaffold_depths = {}
        with open(depth_file, 'r') as infile:
            next(infile)  # Ignorer le header
            for line in infile:
                columns = line.strip().split('\t')
                scaffold = columns[0]
                depth = float(columns[6])  # Mean depth
                scaffold_depths[scaffold] = depth

        for candidat_data in species_data.values():
            scaffold = candidat_data['scaffold']
            if scaffold in scaffold_depths:
                depth = scaffold_depths[scaffold]
            else:
                depth = "NA"  # Mettre "NA" si la profondeur n'est pas disponible
            candidat_data['scaffold_depth'] = depth

def add_species_group(candidates_info, horizon_species_file, dtol_species_file,
                      ncbi_diptera_species_file, ncbi_hymenoptera_species_file):
    # Load species lists
    with open(horizon_species_file, 'r') as infile:
        horizon_species = set(line.strip() for line in infile)
    with open(dtol_species_file, 'r') as infile:
        dtol_species = set(line.strip() for line in infile)
    with open(ncbi_diptera_species_file, 'r') as infile:
        ncbi_diptera_species = set(line.strip() for line in infile)
    with open(ncbi_hymenoptera_species_file, 'r') as infile:
        ncbi_hymenoptera_species = set(line.strip() for line in infile)

    for species, species_data in candidates_info.items():
        if species in horizon_species:
            group = "Tachinidae"
            origin = "Horizon"
        elif species in dtol_species:
            group = "Tachinidae"
            origin = "DTOL"
        elif species in ncbi_diptera_species:
            group = "Free_schizophora"
            origin = "NCBI"
        elif species in ncbi_hymenoptera_species:
            group = "Hymenoptera"
            origin = "NCBI"
        else:
            group = "Unknown"
            origin = "Unknown"

        for candidat_data in species_data.values():
            candidat_data['group'] = group
            candidat_data['origin'] = origin

def add_cluster(candidates_info, cluster_file):
    # Dictionnaire pour stocker les informations sur les clusters
    clusters = {}
    # Ouvrir le fichier de clusters une seule fois
    with open(cluster_file, 'r') as infile:
        for line in infile:
            cluster_name, sequence = line.strip().split('\t')
            clusters[sequence] = cluster_name

    for species, species_data in candidates_info.items():
        for candidate, candidate_data in species_data.items():
            if candidate in clusters:
                candidate_data['cluster'] = clusters[candidate]

def create_candidates_table(candidates_info):
    table = []
    for species, species_data in candidates_info.items():
        for candidate, data in species_data.items():
            viral_protein_accession = data.get('viral_protein_accession', 'NA')
            viral_protein_name = data.get('viral_protein_name', 'NA')
            virus_species = data.get('virus_species', 'NA')
            virus_family = data.get('virus_family', 'NA')
            virus_genomic_structure = data.get('virus_genomic_structure', 'NA')
            scaffold = data.get('scaffold', 'NA')
            scaffold_depth = data.get('scaffold_depth', 'NA')
            BUSCO_count = data.get('BUSCO_count', 'NA')
            TE_count = data.get('TE_count', 'NA')
            evalue = data.get('evalue', 'NA')
            group = data.get('group', 'NA')
            origin = data.get('origin', 'NA')
            start = data.get('start', 'NA')
            end = data.get('end', 'NA') 
            strand = data.get('strand', 'NA')
            cluster = data.get('cluster', 'NA')
            
            row = [candidate, group, origin, species, scaffold, start, end, strand, evalue, viral_protein_accession, viral_protein_name, virus_species, virus_family, virus_genomic_structure, BUSCO_count, TE_count, scaffold_depth, cluster]
            table.append(row)
    return table

def main():
    # Step 1: Process candidates file
    candidates_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidates_summary_table/Mmseqs_candidats_Refseq_virus/result_mmseqs2.m8"
    print("Processing candidates file:", candidates_file)
    candidates_info = Get_candidates_info(candidates_file)

    #print_candidates_info(candidates_info)

    # Step 2: Add virus information
    virus_info_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Candidates_summary_table/Viral_information.tsv"
    print("Adding virus information from file:", virus_info_file)
    Get_virus_info(virus_info_file, candidates_info)
    #print_candidates_info(candidates_info)

    # Step 3: Get BUSCO count
    busco_folder = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Stats/BUSCO/Scaffolds_BUSCO_Count"
    print("Getting BUSCO count from folder:", busco_folder)
    Get_BUSCO_count(candidates_info, busco_folder)
    #print_candidates_info(candidates_info)

    # Step 4: Get TE count
    te_folder = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/TE/TE_count"
    print("Getting TE count from folder:", te_folder)
    Get_TE_count(candidates_info, te_folder)
    #print_candidates_info(candidates_info)

    # Step 5: Add scaffold depth
    depth_folder = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Depth/Diptera_mappings/COV"
    print("Adding scaffold depth from folder:", depth_folder)
    add_scaffold_depth(candidates_info, depth_folder)
    #print_candidates_info(candidates_info)

    # Step 6: Add species group and origin
    horizon_species_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_horizon.txt"
    dtol_species_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_tachinid_ncbi.txt"
    ncbi_diptera_species_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_diptera.txt"
    ncbi_hymenoptera_species_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_hymenoptera.txt"
    print("Adding species group and origin from files:", horizon_species_file, dtol_species_file,
          ncbi_diptera_species_file, ncbi_hymenoptera_species_file)
    add_species_group(candidates_info, horizon_species_file, dtol_species_file,
                      ncbi_diptera_species_file, ncbi_hymenoptera_species_file)
    #print_candidates_info(candidates_info)

    # Step 7: Add cluster information
    cluster_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Clustering/Filtered_Clusters/Candidate_clusters_kept_filtered.m8"
    print("Adding cluster information from file:", cluster_file)
    add_cluster(candidates_info, cluster_file)
    #print_candidates_info(candidates_info)

    # Create candidates table
    print("Creating candidates table...")
    table = create_candidates_table(candidates_info)

    # Convert table to DataFrame
    columns = ['Candidat', 'Group', 'Origin', 'Species', 'Scaffold', 'Start', 'End', 'Strand', 'Evalue', 'Virus_accession', 'Viral_protein', 'Virus_species', 'Virus_family', 'Virus_structure', 'BUSCO_count', 'TE_count', 'Scaffold_depth', 'Cluster']
    df = pd.DataFrame(table, columns=columns)

    # Output DataFrame to a TSV file
    output_file = "/beegfs/data/soukkal/Thesis/Tachinid_Project/Results/Genomic_environment/Candidates_info.tsv"
    print("Saving DataFrame to file:", output_file)
    df.to_csv(output_file, sep='\t', index=False)

if __name__ == "__main__":
    main()