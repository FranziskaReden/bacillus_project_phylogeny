from Bio import Phylo
import argparse
from io import StringIO
import pandas as pd

def get_distance(tree, leave1, leave2):
     return tree.distance(leave1, leave2)

def calculate_balistic_distances(tree, target):

    results_dict = {}
    minimal_dict = {}

    for clade in tree.get_terminals():
        if 'assembly' in clade.name:
            results_dict[clade.name] = {}
            tmp_min = 10000
            minimal_dict[clade.name] = []
            for clade2 in tree.get_terminals():
                if clade2.name != clade.name:
                    distance = float(get_distance(tree, clade.name, clade2.name))
                    results_dict[clade.name][clade2.name] = distance
                    if distance < tmp_min and 'assembly' not in clade2.name:
                        tmp_min = distance
                        minimal_dict[clade.name] = [[clade2.name, distance]]
                    elif distance == tmp_min and 'assembly' not in clade2.name:
                        minimal_dict[clade.name].append([clade2.name, distance])

    df = pd.DataFrame.from_dict(results_dict, orient='index')
    df.to_csv(target, sep='\t', index = True)

    with open(target.replace('.tsv', '_minimal.tsv'), 'w') as w:
        w.write('assembly\tclosest_putative_relative\tdistance\n')
        for key, items in minimal_dict.items():
            w.write(f'{key}\t')
            for item in items:
                w.write(f'{item[0]}\t{item[1]}\t')
            w.write('\n')

def read_tree(treefile):
    with open(treefile) as t: 
        newick = t.readlines()[0]
    newick.strip("\n")

    tree = Phylo.read(StringIO(newick), "newick")

    return tree

def main():
    file_name = 'species_tree/species_tree_support_rename.nwk'
    tree = read_tree(file_name)
    calculate_balistic_distances(tree, 'species_tree/ballistic_distances.tsv')

if __name__ == "__main__":
    main()