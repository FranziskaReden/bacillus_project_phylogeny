import pandas as pd
import os
import argparse
from tqdm import tqdm
import numpy as np

from utils import calculate_sequence_identity, read_fasta, get_sequences_by_names

def find_closest_genome(assembly, genome_names):

    # Read in marker gene of assembly
    marker_genes_assembly = read_fasta(f'bacillus_assemblies/prokka/{assembly}_assembly/{assembly}_assembly_marker_genes.faa')
    results = {}

    # Create resulst dict
    for genome in genome_names:
        results[genome] = []

    # Calculate identity of each marker gene to each genome in the relevant genus
    for marker, seq in marker_genes_assembly.items():
        marker_genes = read_fasta(f'marker_genes/all/faa/{marker}.faa')
        for id, seq2 in marker_genes.items():
            if id in genome_names:
                identity = calculate_sequence_identity(seq, seq2)
                results[id].append(identity)

    # Calculate average identity for each genome
    for key, item in results.items():
        if len(item) > 0:
            results[key] = np.mean(item)
        else:
            results[key] = 0.0

    # Return genome with the highest average identity
    closest_hit = max(results, key=results.get)

    return closest_hit
    
def get_all_rep_names():

    rep_names = []
    
    for marker_genes in os.listdir('marker_genes/reps/faa/'):
        rep_names = []
        fasta_file = read_fasta(f'marker_genes/reps/faa/{marker_genes}')
        for key in fasta_file.keys():
            if key not in rep_names:
                rep_names.append(key)

    return set(rep_names)

def get_identity(seq1, seq2):
    '''
    Calculate the sequence identity between two sequences.
    Parameters:
        seq1 (str): First sequence.
        seq2 (str): Second sequence.
    Returns:
        float: Sequence identity as a fraction between 0 and 1.
    '''
    matches = 0
    for i in range (len(seq1)):
        # Ignore gaps
        if seq1[i] != '-' and seq2[i] != '-':
            if seq1[i] == seq2[i]:
                matches += 1
    # Return relative identity (relative to average length of the two sequences)
    return matches/((len(seq1.replace('-', '')) + len(seq2.replace('-', '')))/2)

def calculate_identity_rep_genomes(rep_genomes):

    identity_matrix = {}
    for a in rep_genomes:
        identity_matrix[a] = {}
        for b in rep_genomes:
            identity_matrix[a][b] = []

    markers = os.listdir('marker_genes/reps/faa/alignment')
    markers_dict = {}
    for marker in markers:
        markers_dict[marker] = read_fasta(f'marker_genes/reps/faa/alignment/{marker}')

    for i in range(0, len(rep_genomes)):
        print(i)
        for j in range(i+1, len(rep_genomes)):
            a = rep_genomes[i]
            b = rep_genomes[j]
            for marker in markers:
                fasta_file = markers_dict[marker]
                if a in fasta_file.keys() and b in fasta_file.keys():
                    identity = get_identity(fasta_file[a], fasta_file[b])
                    identity_matrix[a][b].append(identity)
                    identity_matrix[b][a].append(identity)

    print(identity_matrix)

    for key in identity_matrix.keys():
        for key2 in identity_matrix[key].keys():
            if len(identity_matrix[key][key2]) > 0:
                identity_matrix[key][key2] = np.mean(identity_matrix[key][key2])
            else:
                identity_matrix[key][key2] = 1

    return identity_matrix

def cluster_genomes(identity_dict:dict, relatives:list, closest:list, threshold:float=0.95):
    '''
    Cluster genomes based on a sequence identity threshold.
    Parameters:
        df (pd.DataFrame): DataFrame containing pairwise sequence identities.
        threshold (float): Sequence identity threshold for clustering.
    Returns:
        list: List of clusters, where each cluster is a list of genome names.
    '''
    genomes = set(list(identity_dict.keys())+relatives+closest)

    check = {}
    for key in identity_dict.keys():
        check[key] = False

    genomes_list = list(identity_dict.keys())
    for i in range (len(genomes_list)):
        a = genomes_list[i]
        if check[a] is True:
            continue
        check[a] = True
        cluster = [a]
        for j in range (i+1, len(genomes_list)):
            b = genomes_list[j]
            if identity_dict[a][b] >= threshold:
                cluster.append(b)
                check[b] = True
        
        average_distance = {}
        if len(cluster) > 1:
            for clu1 in cluster:
                average_distance[clu1] = []
                for clu2 in cluster:
                    if clu1 != clu2:
                        average_distance[clu1].append(identity_dict[clu1][clu2])
                average_distance[clu1] = np.mean(average_distance[clu1])

            representative = min(average_distance, key=average_distance.get)
            for clu in cluster:
                if clu != representative and clu not in relatives and clu not in closest:
                    genomes.discard(clu)

    return genomes

def prepare_ali_files(genomes:list, assembly_groups:dict):
    '''
    Prepare alignment files for the selected genomes.
    Parameters:
        genomes (list): List of genome names to include in the alignment files.
        assembly_groups (dict): Dictionary with assembly groups to keep or remove.
    '''

    os.mkdir('alignments', exist_ok=True)
    os.mkdir('alignments/faa', exist_ok=True)
    os.mkdir('alignments/fna', exist_ok=True)

    # Add all representative genomes to the alignment files
    for marker in os.listdir('marker_genes/reps/faa/'):
        sequences_aa = get_sequences_by_names(genomes, f'marker_genes/reps/faa/{marker}')
        sequences_dna = get_sequences_by_names(genomes, f'marker_genes/reps/faa/{marker}')

        # Add all assemblies to the alignment files
        for assembly in os.listdir('bacillus_assemblies/prokka/'):
            if assembly not in assembly_groups['remove']:
                marker_strip = marker.replace(f'.{marker.split(".")[-1]}', '')
                assembly_sequences_aa = get_sequences_by_names([marker_strip], 
                                f'bacillus_assemblies/prokka/{assembly}/{assembly}_marker_genes.faa')
                assembly_sequences_dna = get_sequences_by_names([marker_strip], 
                                f'bacillus_assemblies/prokka/{assembly}/{assembly}_marker_genes.fna')
                
                # If marker gene is present in the assembly, add it to the sequences to be aligned
                if len(assembly_sequences_aa) > 0:
                    # Check if the assembly should be renamed to a group name
                    if assembly in assembly_groups['keep'].keys():
                        sequences_aa[assembly_groups['keep'][assembly]] = assembly_sequences_aa[marker_strip]
                        sequences_dna[assembly_groups['keep'][assembly]] = assembly_sequences_dna[marker_strip]
                    # Otherwise, just add it with its original name
                    else:
                        sequences_aa[assembly] = assembly_sequences_aa[marker_strip]
                        sequences_dna[assembly] = assembly_sequences_dna[marker_strip]

        # Write the sequences to the alignment files    
        with open(f'alignments/faa/{marker}', 'w') as w:
            for key, seq in sequences_aa.items():
                w.write(f'>{key}\n{seq}\n')
        with open(f'alignments/fna/{marker}', 'w') as w:
            for key, seq in sequences_dna.items():
                w.write(f'>{key}\n{seq}\n')
    
def main():
    
    # Read in taxonomy
    taxonomy = pd.read_csv("gtdb/bac120_taxonomy_expanded.tsv", sep='\t')

    # Read in relevant genera for downstream analysis
    relevant_genera = pd.read_csv("config_genera.csv", sep=',')

    rep_genomes = get_all_rep_names()
    print(f'Number of representative genomes: {len(rep_genomes)}')

    identity_dict = calculate_identity_rep_genomes(list(rep_genomes))
    df = pd.DataFrame(identity_dict)
    df.to_csv('marker_genes/reps/identity_matrix.tsv', sep='\t', index=True)

    # Get the closest genome (in the relevant genus) for each assembly
    closest_genomes = []
    for row in relevant_genera.itertuples():
        genome_names = taxonomy[taxonomy['genus'] == row.genus]['accession'].to_list()
        print(f'Finding closest genome for assembly {row.assembly}...')
        closest_genome = find_closest_genome(row.assembly, genome_names)
        print(f'Closest genome to assembly {row.assembly} is {closest_genome}.')
        closest_genomes.append(closest_genome)

    # Cluster the representative genomes based on a sequence identity threshold
    clustered_genomes = cluster_genomes(identity_dict, 
                                        list(relevant_genera['closest_putative_relative']), 
                                        closest_genomes, 
                                        threshold=0.95)
    print(f'Number of representative genomes after clustering: {len(clustered_genomes)}')
    
    assembly_groups = {'keep': {'SRL221_assembly': 'SRL221_244_assembly',
                                'SRL662_assembly': 'SRL662_656_658_assembly',
                                'SRL218_assembly': 'SRL218_215_224_assembly', 
                                'SRL398_assembly': 'SRL398_342_assembly'}, 
                       'remove': ['SRL244_assembly', 'SRL656_assembly', 
                                  'SRL658_assembly', 'SRL215_assembly', 
                                  'SRL224_assembly', 'SRL342_assembly']}
    
    prepare_ali_files(clustered_genomes, assembly_groups)
    return

if __name__ == "__main__":
    main()