import argparse

def fasta_statistics(input_file, output_file):
    total_scaffolds = 0
    total_bases = 0
    total_Ns = 0
    scaffold_lengths = []
    gc_count = 0

    with open(input_file, "r") as file:
        sequence = ""
        for line in file:
            if line.startswith(">"):
                if sequence:
                    # Process the previous sequence
                    length = len(sequence)
                    scaffold_lengths.append(length)
                    total_bases += length
                    total_Ns += sequence.count('N')
                    gc_count += sequence.count('G') + sequence.count('C')
                    sequence = ""
                total_scaffolds += 1
            else:
                sequence += line.strip()

        # Process the last sequence
        if sequence:
            length = len(sequence)
            scaffold_lengths.append(length)
            total_bases += length
            total_Ns += sequence.count('N')
            gc_count += sequence.count('G') + sequence.count('C')

    scaffold_lengths = sorted(scaffold_lengths, reverse=True)
    longest_scaffold = scaffold_lengths[0]
    shortest_scaffold = scaffold_lengths[-1]
    avg_length = total_bases / total_scaffolds
    gc_percentage = (gc_count / total_bases) * 100

    def calculate_Nxx(scaffold_lengths, threshold):
        cumulative_length = 0
        total_length = sum(scaffold_lengths)
        for i, length in enumerate(scaffold_lengths):
            cumulative_length += length
            if cumulative_length >= total_length * threshold:
                return length, i + 1

    N50, L50 = calculate_Nxx(scaffold_lengths, 0.5)
    N75, L75 = calculate_Nxx(scaffold_lengths, 0.75)
    N90, L90 = calculate_Nxx(scaffold_lengths, 0.9)

    scaffolds_1M = sum(1 for x in scaffold_lengths if x >= 1000000)
    scaffolds_100K = sum(1 for x in scaffold_lengths if x >= 100000)
    scaffolds_10K = sum(1 for x in scaffold_lengths if x >= 10000)
    scaffolds_1K = sum(1 for x in scaffold_lengths if x >= 1000)

    length_1M = sum(x for x in scaffold_lengths if x >= 1000000)
    length_100K = sum(x for x in scaffold_lengths if x >= 100000)
    length_10K = sum(x for x in scaffold_lengths if x >= 10000)
    length_1K = sum(x for x in scaffold_lengths if x >= 1000)

    with open(output_file, "w") as out:
        out.write(f"Total scaffolds: {total_scaffolds}\n")
        out.write(f"Total base (bp): {total_bases}\n")
        out.write(f"Total N (bp): {total_Ns}\n")
        out.write(f"Average length (bp): {avg_length:.2f}\n")
        out.write(f"Longest scaffold (bp): {longest_scaffold}\n")
        out.write(f"Shortest scaffold (bp): {shortest_scaffold}\n")
        out.write(f"L50: {L50}\n")
        out.write(f"N50: {N50}\n")
        out.write(f"L75: {L75}\n")
        out.write(f"N75: {N75}\n")
        out.write(f"L90: {L90}\n")
        out.write(f"N90: {N90}\n")
        out.write(f"GC (%): {gc_percentage:.2f}\n")
        out.write(f"Total scaffolds (>= 1000000 bp): {scaffolds_1M}\n")
        out.write(f"Total scaffolds (>= 100000 bp): {scaffolds_100K}\n")
        out.write(f"Total scaffolds (>= 10000 bp): {scaffolds_10K}\n")
        out.write(f"Total scaffolds (>= 1000 bp): {scaffolds_1K}\n")
        out.write(f"Total length (>= 1000000 bp): {length_1M}\n")
        out.write(f"Total length (>= 100000 bp): {length_100K}\n")
        out.write(f"Total length (>= 10000 bp): {length_10K}\n")
        out.write(f"Total length (>= 1000 bp): {length_1K}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate statistics for a FASTA file.")
    parser.add_argument('-a', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file to save statistics')
    
    args = parser.parse_args()

    # Call the fasta_statistics function with the input and output files
    fasta_statistics(args.input, args.output)

if __name__ == "__main__":
    main()
