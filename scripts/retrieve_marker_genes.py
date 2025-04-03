import pandas as pd
import os
import Bio.SeqIO
import io

from utils import read_taxonomy, subsample_tree, read_pgap, read_in_gff, clean_folder

def expand_taxonomy(file_name:str) -> pd.DataFrame:
    """
    Expand the taxonomy of the GTDB genomes into separate columns for each rank.
    """

    target_file = file_name.split(".")[0] + "_expanded.tsv"
    if os.path.exists(target_file):
        print(f"Reading in {target_file}...")
        return pd.read_csv(target_file, sep='\t')
    
    taxonomy = read_taxonomy(file_name)

    print("""Expanding the taxonomy of the GTDB genomes into separate columns for each rank...""")
    # Create two new columns, name of genome (species name) and rank
    taxonomy['name'] = None
    taxonomy['rank'] = None

    # Create a new column for each taxonomic rank
    abbr = {"d": "domain", "p": "phylum", "c": "class", "o": "order", "f": "family", "g": "genus", "s": "species"}
    for _,item in abbr.items():
        taxonomy[item] = None

    # Iterate over the rows of the dataframe
    for idx in taxonomy.index:
        # Split the lineage into ranks
        tmp_ranks = taxonomy.at[idx, 'lineage'].split(';')
        # Remove empty strings
        if "" in tmp_ranks: 
            tmp_ranks.remove("")

        # Iterate over the ranks
        for rank in tmp_ranks:
            tmp_name = rank.split('__')[1]
            tmp_rank = abbr[rank.split('__')[0]]

            # Update the rank and name of the clade
            taxonomy.at[idx, tmp_rank] = tmp_name

            # Keep updating the rank and name of the clade until the species level
            taxonomy.at[idx, "rank"] = tmp_rank
            taxonomy.at[idx, "name"] = tmp_name
    
    print("Writing the expanded taxonomy to file...")
    taxonomy.to_csv(target_file, sep='\t', index=False)

    return taxonomy

def get_relevant_genomes(rank:str, name:str, taxonomy:pd.DataFrame) -> set:
    """
    Get a subset of the genomes based on a specific rank and name using the taxonomy dataframe
    """
    print(f"Getting a subset of the genomes with rank {rank} and clade {name}...")
    genomes = taxonomy[taxonomy[rank].str.contains(name)]['accession'].to_list()

    return set(genomes)

def subsample_marker_genes(genomes:set, type='fna', rep = True) -> set:
    """
    Subsample the marker genes of a specific type from the GTDB genomes
    """

    folder_0 = "marker_genes"

    print(f'''Subsampling marker genes of type {type} from {'reps' if rep is True else 'all'} genomes...''')
    if rep == True:
        folder_1 = "reps"
        target_folder = f"{folder_0}/{folder_1}/{type}"
        input_folder = f"gtdb/bac120_marker_genes_reps_r220/{type}"
    else:
        folder_1 = "all"
        target_folder = f"{folder_0}/{folder_1}/{type}"
        input_folder = f"gtdb/bac120_marker_genes_all_r220/{type}"

    if os.path.exists(input_folder) == False:
        raise FileNotFoundError(f"The folder {input_folder} does not exist.")

    for folder in [f'{folder_0}', f'{folder_0}/{folder_1}', target_folder]:
        if os.path.exists(folder) == False:
            os.mkdir(folder)

    marker_genes = []

    for file in os.listdir(input_folder):
        marker_genes.append(file.replace(f'.{file.split('.')[-1]}', ''))
        if os.path.exists(os.path.join(target_folder, file)) == False:
            with open(os.path.join(target_folder, file), 'w') as out_fasta:
                for record in Bio.SeqIO.parse(os.path.join(input_folder, file), "fasta"):
                    if record.id in genomes:
                        Bio.SeqIO.write(record, out_fasta, "fasta")

    return set(marker_genes)
        
def retireve_gene_families(file:str, marker_genes:set) -> pd.DataFrame:

    gene_families = read_pgap(file)
    gene_families['file_name'] = None 

    print("Comparing the gene families to the marker genes...")
    indeces = []
    for marker in marker_genes:
        marker_tmp = marker.replace(f'.{marker.split('.')[-1]}', '')

        tmp_index = gene_families[gene_families['source_identifier'].str.contains(marker_tmp, na=False)].index

        if len(tmp_index) <= 0:
            print(marker)
        for idx in tmp_index:
            gene_families.at[idx, 'file_name'] = marker
            if marker_tmp == 'PF02576':
                gene_families.at[idx, 'gene_symbol'] = 'rimP'
            elif marker_tmp == 'PF03726':
                gene_families.at[idx, 'gene_symbol'] = 'pnp'

            indeces.append(idx)

    return pd.DataFrame(gene_families[['file_name', 'source_identifier', 'gene_symbol']], index = indeces)

def get_annotated_marker_genes(assembly:str, gene_families:pd.DataFrame, input_folder = 'bacillus_assemblies/prokka/'):

    gff = read_in_gff(os.path.join(input_folder, f'{assembly}/{assembly}.gff'))
    merged_marker_genes = pd.merge(gene_families, gff, on = 'gene_symbol', how='left')

    fna_file = os.path.join(input_folder, f'{assembly}/{assembly}.ffn')
    faa_file = os.path.join(input_folder, f'{assembly}/{assembly}.faa')

    if os.path.exists(fna_file) is False:
        raise FileNotFoundError(f"The file {fna_file} does not exist.")

    if os.path.exists(faa_file) is False:
        raise FileNotFoundError(f"The file {faa_file} does not exist.")

    target_file = os.path.join(input_folder, f'{assembly}/{assembly}_marker_genes.fna')
    #if os.path.exists(target_file) == False:
    with open(target_file, 'w') as out_fasta:
        for record in Bio.SeqIO.parse(os.path.join(fna_file), "fasta"):
            if record.id in merged_marker_genes['ID'].to_list():
                idx = merged_marker_genes[merged_marker_genes['ID'] == record.id].index[0]
                record.id = merged_marker_genes.at[idx, 'file_name']
                record.description = record.id
                Bio.SeqIO.write(record, out_fasta, "fasta")
    
    target_file = os.path.join(input_folder, f'{assembly}/{assembly}_marker_genes.faa')
    #if os.path.exists(target_file) == False:
    with open(target_file, 'w') as out_fasta:
        for record in Bio.SeqIO.parse(os.path.join(faa_file), "fasta"):
            if record.id in merged_marker_genes['ID'].to_list():
                idx = merged_marker_genes[merged_marker_genes['ID'] == record.id].index[0]
                record.id = merged_marker_genes.at[idx, 'file_name']
                record.description = record.id
                Bio.SeqIO.write(record, out_fasta, "fasta")

    merged_marker_genes.to_csv(os.path.join(input_folder, f'{assembly}/{assembly}_marker_genes.tsv'),
              sep='\t', index = False)

def get_marker_genes_all_assemblies(folder:str, gene_families:pd.DataFrame):

    for subfolder in os.listdir(folder):
        if subfolder not in ['hmm_PGAP.tsv', 'hmm_PGAP_marker_genes.tsv']:
            if not os.path.exists(os.path.join(folder, f'{subfolder}/{subfolder}_marker_genes.fna')):
                print(f'Retrieving marker genes for assembly {subfolder}...')
                get_annotated_marker_genes(subfolder, gene_families, input_folder=folder)

def write_marker_genes_all_assemblies(assembly_folder):

    # Create all folder (if they do not exist yet)
    folder = 'marker_genes'
    subfolder = f'{folder}/assemblies'
    target_fna = f'{subfolder}/fna'
    target_faa = f'{subfolder}/faa'
    if not os.path.exists(folder):
        os.mkdir(folder)
    if not os.path.exists(subfolder):
        os.mkdir(subfolder)

    print(f'Writing assembly marker genes into folder {subfolder}...')
    
    # Ensure that target folder are empty.
    # Delete files if need be.
    clean_folder(target_fna)
    clean_folder(target_faa)

    # For each marker gene, write marker gene of assembly in target folder and file
    for file in os.listdir(assembly_folder):
        if 'assembly' in file:
            marker_genes_fna = os.path.join(assembly_folder, f'{file}/{file}_marker_genes.fna')
            marker_genes_faa = os.path.join(assembly_folder, f'{file}/{file}_marker_genes.faa')

            # For DNA sequences
            for record in Bio.SeqIO.parse(marker_genes_fna, "fasta"):
                with open(os.path.join(target_fna, f'{record.id}.fna'), 'a') as out_fasta:
                    record.id = file
                    record.description = record.id
                    Bio.SeqIO.write(record, out_fasta, "fasta")

            # For AA sequences
            for record in Bio.SeqIO.parse(marker_genes_faa, "fasta"):
                with open(os.path.join(target_faa, f'{record.id}.faa'), 'a') as out_fasta:
                    record.id = file
                    record.description = record.id
                    Bio.SeqIO.write(record, out_fasta, "fasta")

def main():
    taxonomy = expand_taxonomy("gtdb/bac120_taxonomy.tsv")
    genomes = get_relevant_genomes('class', 'Bacilli', taxonomy)

    subsample_tree("gtdb/bac120.tree", genomes, taxonomy)

    marker_genes = subsample_marker_genes(genomes, 'faa')
    _ = subsample_marker_genes(genomes, 'fna')
    _ = subsample_marker_genes(genomes, 'faa', rep=False)
    _ = subsample_marker_genes(genomes, 'fna', rep=False) 

    gene_families = retireve_gene_families('bacillus_assemblies/prokka/hmm_PGAP.tsv', marker_genes)
    gene_families.to_csv('hmm_PGAP_marker_genes.tsv', sep='\t', index = False)
    get_marker_genes_all_assemblies('bacillus_assemblies/prokka', gene_families)

    write_marker_genes_all_assemblies('bacillus_assemblies/prokka')
    
if __name__ == "__main__":
    main()