from Bio import SeqIO
import argparse
import os

def create_dict(records, translate=False):

    seq_dict = {}
    length = 0

    for seq_record in records:

        if translate is True: 
            seq_dict[str(seq_record.id).strip('\n')] = str(seq_record.seq.translate())
        
        else: 
            seq_dict[str(seq_record.id).strip('\n')] = str(seq_record.seq)
            length = len(str(seq_record.seq))

    return seq_dict, length

def overlay_ali(aa_aligned, dna_seq): 

    dna_aligned = ""

    tmp = 0
    for aa in aa_aligned: 
        if aa == "-": 
            dna_aligned += "---"
        else: 
            dna_aligned += dna_seq[tmp:tmp+3]
            tmp += 3
    
    if tmp == len(dna_seq)-3:
        dna_aligned += dna_seq[tmp:tmp+3]
    elif tmp == len(dna_seq):
         dna_aligned += "---"
    else:
        print('WARNING! No stop codon (or lack thereof)! Different length!')

    return dna_aligned                

def write_dna_alignment(aa_file_name, dna_file_name, file_name_aligned):

    fna_dict, length_dna = create_dict(SeqIO.parse(dna_file_name, 'fasta'))
    faa_dict, length_aa = create_dict(SeqIO.parse(aa_file_name, 'fasta'))

    empty_seq = "-"*(length_aa*3+3)
    with open(file_name_aligned, "w") as w:  
        for genome in fna_dict.keys():
            w.write(">"+genome+"\n")
            dna_aligned = overlay_ali(faa_dict[genome], fna_dict[genome])
            if len(dna_aligned) != len(empty_seq): 
                print("Warning: Sequences of different length!")
            if dna_aligned.replace('-', '') != fna_dict[genome]:
                print(f'Aligned\t\t{dna_aligned.replace('-', '')}\nOriginal\t{fna_dict[genome]}')
            w.write(dna_aligned+"\n")

def main():

    parser = argparse.ArgumentParser(description="Overlay AA alignment onto DNA sequences.")
    parser.add_argument("--marker", "-m", type=str, help="Name of the assembly.")
    args = parser.parse_args()

    print(f'Overlay AA alignment onto DNA sequences for marker {args.marker}...')

    dna_folder = 'alignments/fna/muscle/'
    if not os.path.exists(dna_folder):
        print(f'Create folder {dna_folder}...')
        os.mkdir(dna_folder)

    write_dna_alignment(f'alignments/faa/muscle/{args.marker}.faa.afa', 
                        f'alignments/fna/{args.marker}.fna',
                        f'alignments/fna/muscle/{args.marker}.fna.afa')
    
    print(f'DNA alignment was written into alignments/fna/muscle/{args.marker}.fna.afa')
    
    return 0

if __name__ == "__main__":
    main()