import os
import numpy as np
from Bio import AlignIO

def calculate_average_alignment_length(folder_path):
    # List all files in the folder that end with .afa
    alignment_files = [f for f in os.listdir(folder_path) if f.endswith('.afa')]
    
    lengths = []
    alignments = len(alignment_files)

    print(f'Number of alignments deteced: {alignments}')
    
    for alignment_file in alignment_files:
        file_path = os.path.join(folder_path, alignment_file)
        
        # Parse the alignment file
        alignment = AlignIO.read(file_path, "fasta")
        
        # Sum the length of the alignments
        for record in alignment:
            lengths.append(len(record.seq))
            break

    print(f'Number of alignments deteced: {len(lengths)}')

    # Calculate the average length
    return np.mean(lengths)

# Provide the path to your folder
folder_path = "alignments/reduced/fna/"
average_length = calculate_average_alignment_length(folder_path)

print(f"Average alignment length: {average_length:.2f} bp")