import argparse
import pandas as pd
import os
from Bio import SeqIO
import numpy as np
from tqdm import tqdm

from utils import filter_fasta_by_names, retrieve_all_reps_names, subsample_tree, combine_fasta

def retrieve_chosen_genomes(genomes:set, input_folder = "marker_genes/all/fna/", target_folder = "alignments/fna"):

    os.makedirs(target_folder.split('/')[0], exist_ok=True)
    os.makedirs(target_folder, exist_ok=True)
    
    print(f'Retrieve chosen genomes...')
    for file in os.listdir(input_folder):
        file_path = os.path.join(input_folder, file)
        output_path = os.path.join(target_folder, file)

        length = filter_fasta_by_names(genomes, file_path, output_path)

    return length

def get_closest_genome(assembly:str, folder) -> set:

    target_file = os.path.join(folder, f'{assembly}/{assembly}_closest_genomes.tsv')
    print(f'Retrieving closest genome for assembly {assembly}...')

    if os.path.exists(target_file):
        results = pd.read_csv(target_file, sep='\t')
        return results.at[results['identity'].idxmax(), 'genome']
    else:
        raise FileNotFoundError(f"The file {target_file} does not exist.")

def extract_best(assembly, folder, file_name, taxonomy, limit=0.90):

    closest_genomes = []

    identity_matrix = pd.read_csv(os.path.join(folder, f'{assembly}/{assembly}_marker_genes_identity.tsv'), sep='\t', index_col=0)
    identity_matrix = identity_matrix.sort_values(by = 'average_identity', ascending = False)

    best_genomes = identity_matrix[identity_matrix['average_identity'] >= limit].index
    if len(best_genomes) == 0:
        best_genomes = [identity_matrix['average_identity'].idxmax()]
    
    for genome in best_genomes:
        idx = taxonomy[taxonomy['accession'] == genome].index
        lineage = taxonomy.at[idx[0], 'lineage']
        genus = taxonomy.at[idx[0], 'genus']

        #print(f'Best hit: {genome} with average identity of {identity_matrix.at[genome, 'average_identity']} and lineage {lineage}.')
        
        with open(file_name, 'a') as a:
            a.write(f'{assembly}\t{genome}\t{identity_matrix.at[genome, 'average_identity']}\t{lineage}\t{genus}\n')

    best_genome = get_closest_genome(assembly, folder)
    closest_genomes.append(best_genome)

    return closest_genomes

def choose_genomes_according_to_identity(folder:str, taxonomy:pd.DataFrame, file_name:str):

    all_cloeset_genomes = []
    with open(file_name, 'w') as w:
        w.write(f'assembly\tgenome\thit\tlineage\tgenus\n')
    for file in os.listdir(folder):
        if 'assembly' in file:
            print(f'Getting best hits for assembly {file}...')
            closest_genomes = extract_best(file, folder, file_name, taxonomy, limit=0.90)
            for g in closest_genomes:
                if g not in all_cloeset_genomes:
                    all_cloeset_genomes.append(g)
                    print(f'Closest genome {g}')
                else:
                    print(f'Closest genome {g} (already added)')

    return all_cloeset_genomes

def choose_genomes(closest_genomes:list, taxonomy:pd.DataFrame, file_name:str, rep_genomes:list) -> set:

    closest_genomes = set(closest_genomes)

    best_hits = pd.read_csv(file_name, sep='\t')
    genera = best_hits['genus'].unique()

    genomes = []
    for genus in genera:
        pot_genomes = taxonomy[taxonomy['genus'] == genus]['accession'].to_list()
        for genome in pot_genomes:
            if genome not in genomes:
                genomes.append(genome)

    print(f'List of genera: {genera}')
    print(f'Number of genomes in chosen genera: {len(genomes)}.')
    
    genomes = set(genomes)
    rep_genomes = set(rep_genomes)

    genomes = genomes.intersection(rep_genomes)
    print(f'Number of rep genomes in genera: {len(genomes)}')

    return closest_genomes.union(genomes)

def combine_fasta_files():

    dna_folder = 'alignments/fna'
    aa_folder = 'alignments/faa'
    
    print(f'Write marker genes of chosen genomes and assemblies into folder {dna_folder}...')
    for file in os.listdir(dna_folder):
        target_file = os.path.join(dna_folder, file)
        assembly_marker_genes = f'marker_genes/assemblies/fna/{file}'
        if os.path.exists(assembly_marker_genes):
            combine_fasta(target_file, assembly_marker_genes, target_file)
        else:
            print(f'WARNING: file {assembly_marker_genes} does not exist!')

    print(f'Write marker genes of chosen genomes and assemblies into folder {aa_folder}...')
    for file in os.listdir(aa_folder):
        target_file = os.path.join(aa_folder, file)
        assembly_marker_genes = f'marker_genes/assemblies/faa/{file}'
        if os.path.exists(assembly_marker_genes):
            combine_fasta(target_file, assembly_marker_genes, target_file)
        else:
            print(f'WARNING: file {assembly_marker_genes} does not exist!')

def main():
    
    parser = argparse.ArgumentParser(description="Compare marker genes")
    parser.add_argument("--folder", "-f", type=str, default = 'bacillus_assemblies/prokka/', help="Path to folder holding assemblies.")
    args = parser.parse_args()

    if not os.path.exists('gtdb/bac120_marker_genes_reps_r220/names.txt'):
        retrieve_all_reps_names()

    with open('gtdb/bac120_marker_genes_reps_r220/names.txt') as file:
        rep_genomes = [line.strip() for line in file]

    taxonomy = pd.read_csv(f'gtdb/bac120_taxonomy_expanded.tsv', sep='\t')
    file_name = 'marker_genes/best_hits.tsv'

    # Get the closest genomes to the assemblies
    # Write representative genomes with identity >90 into file for all assemblies
    closest_genomes = choose_genomes_according_to_identity(args.folder, taxonomy, file_name)
    print(f'Number of closest genomes: {len(closest_genomes)}')

    # Get all representative genomes of the genera in question
    chosen_genomes = choose_genomes(closest_genomes, taxonomy, file_name, rep_genomes)
    
    selected = retrieve_chosen_genomes(chosen_genomes)
    _ = retrieve_chosen_genomes(chosen_genomes, input_folder = "marker_genes/all/faa/", target_folder = "alignments/faa")
    print(f'Number of selected genomes {selected}.')

    combine_fasta_files()
    subsample_tree(f'gtdb/bac120.tree', chosen_genomes, taxonomy, target_file=f'gtdb/chosen_genomes.tree')

    return 0

if __name__ == "__main__":
    main()