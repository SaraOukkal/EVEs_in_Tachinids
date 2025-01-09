# Charger la liste des espèces à exclure depuis le fichier.txt
with open('/beegfs/data/soukkal/Thesis/Tachinid_Project/Species_list_diptera.txt', 'r') as species_file:
    species_to_exclude = [line.strip().replace('_', ' ') for line in species_file]

# Ouvrir le fichier fasta en lecture
with open('/beegfs/data/soukkal/Thesis/Databases/Refseq_Schizophora.faa', 'r') as fasta_file:
    # Ouvrir le fichier de sortie en écriture
    with open('/beegfs/data/soukkal/Thesis/Databases/Refseq_Schizophora_filtered.faa', 'w') as output_file:
        current_header = None
        exclude_sequence = False
        sequence_lines = []

        # Parcourir chaque ligne du fichier fasta
        for line in fasta_file:
            if line.startswith('>'):
                # Si c'est un en-tête
                if current_header is not None and not exclude_sequence:
                    # Écrire l'en-tête précédent et sa séquence associée dans le fichier de sortie
                    output_file.write(current_header + '\n')  # Ajouter une nouvelle ligne après l'en-tête
                    output_file.write(''.join(sequence_lines) + '\n')  # Écrire la séquence
                current_header = line.strip()
                species_name = current_header.split('[')[-1].split(']')[0].strip().replace('_', ' ')

                # Vérifier si l'espèce doit être exclue
                exclude_sequence = any(species in species_name for species in species_to_exclude)

                # Ajouter la mention "Refseq-Schizophora" à la fin de l'en-tête
                current_header += ' Refseq-Schizophora'

                # Réinitialiser les lignes de séquence
                sequence_lines = []
            else:
                # Construire la séquence
                sequence_lines.append(line.strip())

        # Écrire la dernière séquence si elle ne doit pas être exclue
        if current_header is not None and not exclude_sequence:
            output_file.write(current_header + '\n')
            output_file.write(''.join(sequence_lines) + '\n')

