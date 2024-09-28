import argparse

def fasta_statistics(input_file, output_file):
    total_scaffolds = 0
    total_bases = 0
    total_Ns = 0
    scaffold_lengths = []
    gc_count = 0
    current_scaffold_length = 0

    with open(input_file, "r") as file:
        for line in file:
            line = line.strip()
            if not line:
                continue  # 跳过空行
            if line.startswith(">"):
                if total_scaffolds > 0:
                    # 处理前一个scaffold
                    scaffold_lengths.append(current_scaffold_length)
                    current_scaffold_length = 0
                total_scaffolds += 1
            else:
                seq_line = line.upper()
                length = len(seq_line)
                current_scaffold_length += length
                total_bases += length
                total_Ns += seq_line.count('N')
                gc_count += seq_line.count('G') + seq_line.count('C')

        # 处理最后一个scaffold
        if current_scaffold_length > 0:
            scaffold_lengths.append(current_scaffold_length)

    scaffold_lengths.sort(reverse=True)
    longest_scaffold = scaffold_lengths[0]
    shortest_scaffold = scaffold_lengths[-1]
    avg_length = total_bases / total_scaffolds
    gc_percentage = (gc_count / total_bases) * 100

    def calculate_Nxx(scaffold_lengths, threshold):
        cumulative_length = 0
        total_length = sum(scaffold_lengths)
        target = total_length * threshold
        for i, length in enumerate(scaffold_lengths):
            cumulative_length += length
            if cumulative_length >= target:
                return length, i + 1

    N50, L50 = calculate_Nxx(scaffold_lengths, 0.5)
    N75, L75 = calculate_Nxx(scaffold_lengths, 0.75)
    N90, L90 = calculate_Nxx(scaffold_lengths, 0.9)

    thresholds = [1000000, 100000, 10000, 1000]
    scaffold_counts = {}
    length_sums = {}
    for t in thresholds:
        scaffold_counts[t] = 0
        length_sums[t] = 0
    for length in scaffold_lengths:
        for t in thresholds:
            if length >= t:
                scaffold_counts[t] += 1
                length_sums[t] += length

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
        for t in thresholds:
            out.write(f"Total scaffolds (>= {t} bp): {scaffold_counts[t]}\n")
        for t in thresholds:
            out.write(f"Total length (>= {t} bp): {length_sums[t]}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate statistics for a FASTA file.")
    parser.add_argument('-a', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file to save statistics')
    
    args = parser.parse_args()

    fasta_statistics(args.input, args.output)

if __name__ == "__main__":
    main()
