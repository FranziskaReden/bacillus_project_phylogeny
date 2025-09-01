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
    for i in range (len(aa_aligned)): 
        aa = aa_aligned[i]
        if aa == "-": 
            dna_aligned += "---"
        else:
            dna_aligned += dna_seq[tmp:tmp+3]
            tmp += 3
    
    if tmp == len(dna_seq)-3 and aa_aligned[-1] != "*":
        dna_aligned += dna_seq[tmp:tmp+3]
        tmp += 3
    elif tmp == len(dna_seq) and aa_aligned[-1] != "*":
        dna_aligned += "---"
        tmp += 3
    elif tmp != len(dna_seq):
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
                print(dna_aligned)
                print(empty_seq)
            if dna_aligned.replace('-', '') != fna_dict[genome]:
                print(f'Aligned\t\t{dna_aligned.replace('-', '')}\nOriginal\t{fna_dict[genome]}')
            w.write(dna_aligned+"\n")

def main():

    parser = argparse.ArgumentParser(description="Overlay AA alignment onto DNA sequences.")
    parser.add_argument("--marker", "-m", type=str, help="Name of the assembly.")
    parser.add_argument("--folder", '-f', type=str, default='alignments', help="Folder with the alignments.")
    args = parser.parse_args()

    print(f'Overlay AA alignment onto DNA sequences for marker {args.marker}...')

    dna_folder = os.path.join(args.folder, 'fna/')
    if not os.path.exists(dna_folder):
        print(f'Create folder {dna_folder}...')
        os.mkdir(dna_folder)

    write_dna_alignment(os.path.join(args.folder, f'faa/{args.marker}.faa.afa'), 
                        os.path.join(args.folder, f'fna/{args.marker}.fna'),
                        os.path.join(args.folder, f'fna/{args.marker}.fna.afa'))
    
    print(f'DNA alignment was written into {args.folder}/fna/{args.marker}.fna.afa')
    
    return 0

if __name__ == "__main__":
    main()