from Bio import SeqIO

def split_sequence_by_gap(input_fasta, output_fasta):
    sequences = SeqIO.parse(input_fasta, "fasta")
    
    with open(output_fasta, "w") as output_handle:
        for seq_record in sequences:
            seq = str(seq_record.seq)
            gap_positions = [i for i, letter in enumerate(seq) if letter == 'N']
            split_positions = [0] + [pos + 1 for pos in gap_positions] + [len(seq)]
            
            for i in range(1, len(split_positions)):
                start = split_positions[i-1]
                end = split_positions[i]
                sub_seq = seq[start:end].replace('N', '')
                
                if sub_seq:
                    sub_record = seq_record[start:end]
                    sub_record.seq = sub_seq
                    sub_record.id = f"{seq_record.id}_part{i}"
                    sub_record.description = f"part {i} from {start} to {end}"
                    SeqIO.write(sub_record, output_handle, "fasta")

if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Split a FASTA sequence by GAPs (Ns).")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    
    args = parser.parse_args()
    
    split_sequence_by_gap(args.input, args.output)
