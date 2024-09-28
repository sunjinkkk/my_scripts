import argparse

def fasta_sequence_lengths(input_file, output_file):
    with open(input_file, "r") as file, open(output_file, "w") as out:
        current_length = 0
        for line in file:
            line = line.strip()
            if line.startswith(">"):
                if current_length > 0:
                    out.write(f"{current_length}\n")
                    current_length = 0  # 重置长度
            else:
                current_length += len(line)
        
        # 处理最后一个序列
        if current_length > 0:
            out.write(f"{current_length}\n")

def main():
    parser = argparse.ArgumentParser(description="Calculate sequence lengths in a FASTA file.")
    parser.add_argument('-a', '--input', required=True, help='Input FASTA file')
    parser.add_argument('-o', '--output', required=True, help='Output file to save sequence lengths')
    
    args = parser.parse_args()

    fasta_sequence_lengths(args.input, args.output)

if __name__ == "__main__":
    main()

