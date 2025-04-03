import os
from Bio import SeqIO
from Bio.Align import PairwiseAligner
from Bio import Phylo
import pandas as pd
import shutil

def read_partitions(file:str) -> list:

    partitions = []
    with open(file) as f:
        lines = f.readlines()

    for line in lines:
        limit = line.strip().split(' = ')[-1].split('-')
        partitions.append([int(limit[0]-1), int(limit[1])])

    return partitions

def clean_folder(folder_path):
    """
    Ensure the folder is empty by deleting all files and subdirectories within it.
    If the folder does not exist, it will be created.
    """
    if os.path.exists(folder_path):
        # Iterate through all files and subdirectories in the folder
        for filename in os.listdir(folder_path):
            file_path = os.path.join(folder_path, filename)
            try:
                # Remove files
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                # Remove directories
                elif os.path.isdir(file_path):
                    shutil.rmtree(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")
    else:
        # Create the folder if it doesn't exist
        os.makedirs(folder_path)

def read_taxonomy(file_name:str) -> pd.DataFrame:
    """
    Read the taxonomy of the GTDB genomes
    """
    # Read the taxonomy of the GTDB genomes
    print('Reading in the taxonomy of the GTDB genomes...')
    if os.path.exists(file_name) == False:
        raise FileNotFoundError(f"The file {file_name} does not exist.")
    taxonomy = pd.read_csv(file_name, sep='\t', header=0, names=['accession', 'lineage'])
    return taxonomy

def check_column(attributes:str):

    name = None 
    id = None

    attributes = attributes.split(';')
    for att in attributes:
        att_2 = att.split('=')
        if att_2[0] == 'Name':
            name = att_2[1]
        elif att_2[0] == 'ID':
            id = att_2[1]
    
    return name, id

def read_in_gff(gff_file:str) -> pd.DataFrame:
    """
    Read in Annotation file.
    """

    if os.path.exists(gff_file) == False:
        raise FileNotFoundError(f"The file {gff_file} does not exist.")
    
    # Read in file
    col_names = ["SeqID", "Source", "Type", "Start", "End", "Score", "Strand", "Phase", "Attributes"]
    gff = read_until_comment(gff_file, col_names, delimiter='\t')

    gff['gene_symbol'] = None
    gff['ID'] = None

    # Retrieve gene symbol from attributes column
    for idx in gff.index:
        gff.at[idx, 'gene_symbol'], gff.at[idx,'ID'] = check_column(gff.at[idx, 'Attributes'])

    return gff

def read_until_comment(file, names, delimiter="\t"):

    with open(file, "r") as f:
        lines = f.readlines()
        count = 0
        for line in lines:
            if line[0] == '#' and count > 0:
                break
            elif line[0] != '#': 
                count += 1

    return pd.read_csv(file, sep=delimiter, header = None, names = names, comment = '#', nrows = count)

def read_pgap(file_name:str) -> pd.DataFrame:
    """
    Read the gene families from the PGAP output
    """
    print(f"Reading in the gene families from {file_name}...")
    if os.path.exists(file_name) == False:
        raise FileNotFoundError(f"The file {file_name} does not exist.")
    gene_families = pd.read_csv(file_name, sep = "\t")
    return gene_families

def subsample_tree(tree_file:str, genomes:set, taxonomy:pd.DataFrame, target_file = None):
    """Subsample the tree to include only the genomes of interest"""
    
    target_file=f'{tree_file}.subsample.tree' if not target_file else target_file
    
    if os.path.exists(target_file) == True:
        return
    
    print(f"Subsampling the tree {tree_file}...")
    tree = Phylo.read(tree_file, "newick")
    common_ancestor = tree.common_ancestor(genomes)
    Phylo.write(common_ancestor, target_file, format="newick")

    with open() as t: 
        newick = t.readlines(target_file)[0]
    newick.strip("\n")

    for acc in genomes: 
        newick = newick.replace(acc, taxonomy[taxonomy["accession"] == acc]["name"].values[0])

    with open(f'{target_file}_rename.tree', "w") as w:
        w.write(newick)

def retrieve_all_reps_names():

    seq_set = set()

    path = 'gtdb/bac120_marker_genes_reps_r220/fna/'
    for file in os.listdir(path):
        sequences = read_fasta(os.path.join(path, file))
        seq_set = seq_set.union(sequences.keys())

    with open('gtdb/bac120_marker_genes_reps_r220/names.txt', 'w') as w:
        for name in seq_set: 
            w.write(f'{name}\n')

def calculate_sequence_identity(seq1, seq2):
    aligner = PairwiseAligner()
    aligner.mode = "local"  # Use local alignment instead of global
    
    # Perform local pairwise alignment
    alignments = aligner.align(seq1, seq2)
    
    # Get the best local alignment
    best_alignment = alignments[0]
    
    # Extract the aligned sequences
    aligned_seq1, aligned_seq2 = best_alignment.sequences
    
    # Compute the number of identical matches
    matches = sum(res1 == res2 for res1, res2 in zip(aligned_seq1, aligned_seq2))
    
    # Use the length of the aligned region (excluding leading/trailing gaps)
    aligned_length = max(len(aligned_seq1.replace("-", "")), len(aligned_seq2.replace("-", "")))
    
    # Calculate sequence identity
    identity = matches / aligned_length if aligned_length > 0 else 0
    
    return identity

def read_fasta(file_path):
    """
    Reads a FASTA file and returns a dictionary of sequences as strings. Key is
    the name of the sequence.
    """
    sequences = {}
    for record in SeqIO.parse(file_path, "fasta"):
        sequences[record.id] = str(record.seq)
    return sequences

def read_fasta_object(file_path):
    """Reads a FASTA file and returns a dictionary of SeqRecord objects."""
    sequences = {record.id: record for record in SeqIO.parse(file_path, "fasta")}
    return sequences

def filter_fasta_by_names(seq_ids: set, input_fasta: str, output_fasta: str):
    """
    Extracts sequences from an input FASTA file whose IDs match those in the given set 
    and writes them to an output FASTA file.

    Parameters:
        seq_ids (set): A set of sequence IDs to filter.
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file.

    Returns:
        int: The number of sequences written to the output file.
    """
    # Read and filter sequences
    selected_records = [record for record in SeqIO.parse(input_fasta, "fasta") if record.id in seq_ids]

    # Write to output FASTA file if matches were found
    if selected_records:
        SeqIO.write(selected_records, output_fasta, "fasta")

    return len(selected_records)

def get_sequences_by_names(seq_ids: set, input_fasta: str) -> dict:
    """
    Extracts sequences from an input FASTA file whose IDs match those in the given set
    and returns them as a dictionary.

    Parameters:
        seq_ids (set): A set of sequence IDs to filter.
        input_fasta (str): Path to the input FASTA file.

    Returns:
        dict: A dictionary where keys are sequence IDs and values are sequences as strings.
    """
    sequences = {
        record.id: str(record.seq)
        for record in SeqIO.parse(input_fasta, "fasta")
        if record.id in seq_ids
    }
    return sequences

def combine_fasta(file1, file2, output_file):
    """
    Combine two FASTA files into a single output file.

    Args:
        file1 (str): Path to the first FASTA file.
        file2 (str): Path to the second FASTA file.
        output_file (str): Path to the output FASTA file.
    """
    # Read sequences from both files
    sequences = []
    for fasta_file in [file1, file2]:
        with open(fasta_file, "r") as handle:
            sequences.extend(list(SeqIO.parse(handle, "fasta")))

    # Write all sequences to the output file
    with open(output_file, "w") as output_handle:
        SeqIO.write(sequences, output_handle, "fasta")
