def standardize_fasta(input_fasta, output_fasta, line_length=100):
    with open(input_fasta, 'r') as infile, open(output_fasta, 'w') as outfile:
        sequence = ""
        header = ""
        for line in infile:
            line = line.strip()
            if line.startswith(">"):
                if sequence:
                    # Write the sequence in chunks of line_length
                    for i in range(0, len(sequence), line_length):
                        outfile.write(sequence[i:i+line_length] + '\n')
                header = line
                outfile.write(header + '\n')
                sequence = ""
            else:
                sequence += line
        
        # Write the last sequence
        if sequence:
            for i in range(0, len(sequence), line_length):
                outfile.write(sequence[i:i+line_length] + '\n')

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Standardize FASTA file to have sequences with fixed line length")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")
    parser.add_argument("-l", "--length", type=int, default=100, help="Line length for sequences (default: 100)")

    args = parser.parse_args()

    standardize_fasta(args.input, args.output, args.length)
