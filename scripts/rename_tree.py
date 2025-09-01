import os
from Bio import Phylo
import pandas as pd
from io import StringIO
import re
import math

def extract_support(label):
    match = re.search(r"support=([^\];]+)", label)
    if match:
        try:
            val = float(match.group(1))
            if math.isnan(val):
                return 0.0
            return val
        except ValueError:
            return 0.0
    return None

def read_taxonomy(file_name:str) -> pd.DataFrame:
    """
    Read the taxonomy of the GTDB genomes
    """
    # Read the taxonomy of the GTDB genomes
    print('Reading in the taxonomy of the GTDB genomes...')
    if os.path.exists(file_name) == False:
        raise FileNotFoundError(f"The file {file_name} does not exist.")
    taxonomy = pd.read_csv(file_name, sep='\t')
    return taxonomy

def read_tree(file_name:str, taxonomy:pd.DataFrame):

    colors = [
    "#377eb8",  # Blue
    "#ff7f00",  # Orange
    "#4daf4a",  # Green
    "#f781bf",  # Pink
    "#a65628",  # Brown
    "#984ea3",  # Purple
    "#ffff33",  # Yellow
    "#00ced1",  # Cyan
    "#b2df8a"   # Olive
]
    
    colors_dict = {}

    with open(file_name) as t: 
        newick = t.readlines()[0]
    newick.strip("\n")

    tree = Phylo.read(StringIO(newick), "newick")

    for clade in tree.get_nonterminals():
        if clade.name and "support=" in clade.name:
            support = extract_support(clade.name)
            clade.confidence = support
            clade.name = None

    with open('color_code.tsv', 'w') as w:
        w.write('taxon\tgroup\tcolor\n')
        for clade in tree.get_terminals():
            if 'assembly' not in clade.name:
                species = taxonomy[taxonomy['accession'] == clade.name]['species'].values[0]
                genus = taxonomy[taxonomy['accession'] == clade.name]['genus'].values[0]

                if genus in colors_dict.keys():
                    w.write(f'{clade.name}\t{genus}\t{colors_dict[genus]}\n')
                elif not not colors:
                    color = colors.pop(0)
                    colors_dict[genus] = color
                    w.write(f'{clade.name}\t{genus}\t{color}\n')
            else:
                w.write(f'{clade.name}\tassembly\t#e41a1c\n')

    Phylo.write(tree, "species_tree_support.nwk", "newick")

    with open('color_code_rename.tsv', 'w') as w:
        w.write('taxon\tgroup\tcolor\n')
        for clade in tree.get_terminals():
            if 'assembly' not in clade.name:
                species = taxonomy[taxonomy['accession'] == clade.name]['species'].values[0]
                genus = taxonomy[taxonomy['accession'] == clade.name]['genus'].values[0]
                clade.name = species

                if genus in colors_dict.keys():
                    w.write(f'{clade.name}\t{genus}\t{colors_dict[genus]}\n')
                elif not not colors:
                    color = colors.pop(0)
                    colors_dict[genus] = color
                    w.write(f'{clade.name}\t{genus}\t{color}\n')
            else:
                w.write(f'{clade.name}\tassembly\t#e41a1c\n')

    Phylo.write(tree, "species_tree_support_rename.nwk", "newick")

def main():
    taxonomy = read_taxonomy("gtdb/bac120_taxonomy_expanded.tsv")
    treefile = "species_tree.nwk"

    read_tree(treefile, taxonomy)
    return 0 

if __name__ == "__main__": 
    main()