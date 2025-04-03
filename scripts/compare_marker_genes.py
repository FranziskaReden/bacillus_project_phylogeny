import pandas as pd
import os
import argparse
from tqdm import tqdm
import numpy as np

from utils import calculate_sequence_identity, read_fasta, get_sequences_by_names

def read_marker_genes_assembly(assembly, folder, reference) -> str:
    '''
    Compare the marker gene of the assembly with the marker genes of the representative 
    genomes. This will give an indication for in which genus we can expect the 
    assembly to fall. The results are written into file 
    f'{assembly}/{assembly}_marker_genes_identity.tsv'.

    Input: Name of assembly, name of folder of assembly, path to representative marker genes

    Returns: Path to and name of results file.
    '''

    results_file = os.path.join(folder, f'{assembly}/{assembly}_marker_genes_identity.tsv')
    if os.path.exists(results_file):
        return results_file
    
    # Read in marker genes found in assembly
    asembly_seqs = read_fasta(os.path.join(folder, f'{assembly}/{assembly}_marker_genes.fna'))
    # Create empty dataframe to hold results, write header into file
    identity_matrix = pd.DataFrame(columns = asembly_seqs.keys())
    identity_matrix.to_csv(results_file, sep='\t', index = True)
    # Iterate over each marker gene in assembly
    for seq_id, sequence in asembly_seqs.items():
        # Read in correponding marker genes from the gtdb
        #print(f"Calculating sequence identity of {seq_id}...")
        gtdb_marker_genes = read_fasta(os.path.join(reference, f'{seq_id}.fna'))
        # Compare 
        for seq_id2, sequence2 in tqdm(gtdb_marker_genes.items(), desc=f"Comparing {seq_id}", leave=False):
            identity = calculate_sequence_identity(sequence, sequence2)
            identity_matrix.loc[seq_id2, seq_id] = identity

        # Write results
        identity_matrix.to_csv(results_file, sep='\t', index = True)

    # Calculate the average identity for each row, ignoring NaN values
    identity_matrix['average_identity'] = identity_matrix.mean(axis=1, skipna=True)
    # Write results
    identity_matrix.to_csv(results_file, sep='\t', index = True)

    return results_file

def get_closest_genome(assembly:str, folder:str, genus:str, taxonomy:pd.DataFrame) -> set:
    '''
    Search for closest genomes given the marker genes (all genomes, not only representative) 
    given the most porbable genus of the assembly has already been found.
    Will calculate the identity between assembly marker genes and reference marker genes.
    The reuslts are written into closest_genomes.tsv file.
    
    Input: name of assembly, genus name, taxonomy table

    Returns: The name of the clostest genome
    '''

    target_file = os.path.join(folder, f'{assembly}/{assembly}_closest_genomes.tsv')
    print(f'Retrieving closest genome for assembly {assembly}...')

    if os.path.exists(target_file):
        results = pd.read_csv(target_file, sep='\t')
        return results.at[results['identity'].idxmax(), 'genome']

    genome_names = set(taxonomy[taxonomy['genus'] == genus]['accession'])
    print(f'Number of potential genomes {len(genome_names)}...')
    asembly_seqs = read_fasta(os.path.join(folder, f'{assembly}/{assembly}_marker_genes.fna'))

    identities = {}
    for seq_id, sequence in tqdm(asembly_seqs.items()):
        # Read in correponding marker genes from the gtdb
        #print(f"Calculating sequence identity of {seq_id}...")
        gtdb_marker_genes = get_sequences_by_names(genome_names, os.path.join('marker_genes/all/fna/', f'{seq_id}.fna'))
        # Compare 
        for seq_id2, sequence2 in gtdb_marker_genes.items():
            identities.setdefault(seq_id2, []).append(calculate_sequence_identity(sequence, sequence2))

    with open(target_file, 'w') as w:
        w.write('genome\tidentity\n')
        for key,item in identities.items():
            w.write(f'{key}\t{np.mean(item)}\n')

    max_identity = 0
    best_genome = None
    for key,item in identities.items():
        if np.mean(item) > max_identity:
            max_identity = np.mean(item)
            best_genome = key

    return best_genome

def extract_best(assembly, folder, taxonomy, identity_file, limit=0.90):
    '''
    Will extract the closest reference genomes to the assembly in question according to 
    the gtdb marker genes. All representative genomes with identity limit will be 
    considered. The genus of said genomes are then assumed to be the same for the 
    assembly genome. All representative genomes in said genus will be chosen in next step.
    Given the likely genus of the assembly, the average identity between the assembly 
    marker genes and the marker genes of ALL (not only representative) genomes
    falling into said genus will be calculates. The results will be written into 
    file '{assembly}/{assembly}_closest_genomes.tsv'. The closest genome will be 
    further considered in tthe downstream analysis.

    Input: 
        - Name of assembly
        - folder in which to find the assembly
        - file holding the average identities calculated for representative genomes
        - limit (default 90).
    
    Output: Name of closest genome
    '''

    selected_genera = []
    closest_genomes = []

    identity_matrix = pd.read_csv(identity_file, sep='\t', index_col=0)
    identity_matrix = identity_matrix.sort_values(by = 'average_identity', ascending = False)

    best_genomes = identity_matrix[identity_matrix['average_identity'] >= limit].index
    if len(best_genomes) == 0:
        best_genomes = [identity_matrix['average_identity'].idxmax()]
    
    for genome in best_genomes:
        idx = taxonomy[taxonomy['accession'] == genome].index
        genus = taxonomy.at[idx[0], 'genus']

        if genus not in selected_genera:
            selected_genera.append(genus)

    if len(selected_genera) > 1:
        print(f"Warning! Two selected genera for assembly: {selected_genera}")
    for genus in selected_genera:
        best_genome = get_closest_genome(assembly, folder, genus, taxonomy)
        closest_genomes.append(best_genome)

    return closest_genomes

def main():
    
    parser = argparse.ArgumentParser(description="Compare marker genes")
    parser.add_argument("--assembly", "-a", type=str, help="Name of the assembly.")
    parser.add_argument("--folder", "-f", type=str, default = 'bacillus_assemblies/prokka/', help="Path to folder holding assemblies.")
    parser.add_argument("--reference", '-r', type=str, default = 'marker_genes/reps/fna/', help="Path to folder holding reference marker genes.")
    parser.add_argument("--taxonomy", '-t', type=str, default = "gtdb/bac120_taxonomy_expanded.tsv", help= "Path to file holding taxonomy table.")
    args = parser.parse_args()
    
    assembly = args.assembly
    print(f"Assembly: {assembly}")

    taxonomy = pd.read_csv(args.taxonomy, sep='\t')

    results_file = read_marker_genes_assembly(assembly, args.folder, args.reference)
    closest_genome = extract_best(assembly, args.folder, taxonomy, results_file)

    print(f'Closest genome to assembly {assembly} is {closest_genome}.')

    return 0

if __name__ == "__main__":
    main()