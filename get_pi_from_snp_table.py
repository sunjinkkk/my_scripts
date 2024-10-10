import pandas as pd

def calculate_pi_in_bins(chrom_data, samples, bin_size=100000, step_size=10000):
    """Calculate number of variants and mean pi value in specified bins for a chromosome."""
    results = []
    start = chrom_data['POS'].min()
    end = chrom_data['POS'].max()
    
    # Slide the bin across the chromosome
    for bin_start in range(start, end, step_size):
        bin_end = bin_start + bin_size
        bin_data = chrom_data[(chrom_data['POS'] >= bin_start) & (chrom_data['POS'] < bin_end)]
        
        if not bin_data.empty:
            # Filter the genotypes based on the provided samples
            genotype_data = bin_data[samples]
            # Calculate allele counts
            allele_counts = genotype_data.apply(lambda x: x.value_counts(), axis=1).fillna(0)
            # Calculate pi for each SNP
            pi_values = []
            for index, counts in allele_counts.iterrows():
                if len(counts) >= 2:
                    major_allele_count = counts.iloc[0]
                    minor_allele_count = counts.iloc[1]
                else:
                    major_allele_count = counts.iloc[0]
                    minor_allele_count = 0
                pi_value = 2 * (major_allele_count / (major_allele_count + minor_allele_count)) * \
                           (minor_allele_count / (major_allele_count + minor_allele_count))
                pi_values.append(pi_value)
            
            # Compute average pi for the window
            mean_pi = pd.Series(pi_values).mean()
            n_variants = len(bin_data)
            results.append({
                'CHROM': chrom_data['chrom'].iloc[0],
                'BIN_START': bin_start,
                'BIN_END': bin_end,
                'N_VARIANTS': n_variants,
                'PI': mean_pi
            })
    
    return results

def main(snp_file, group_file, bin_size, step_size, output_file):
    # Read the SNP data
    snp_data = pd.read_csv(snp_file, sep='\t')
    # Read the group file
    groups = pd.read_csv(group_file, header=None)
    samples = groups[0].tolist()
    
    # Filter SNP data based on the samples
    genotype_columns = ['chrom', 'POS'] + samples
    snp_data = snp_data[genotype_columns]
    
    # Calculate pi for each chromosome
    bin_results = []
    for chrom, group in snp_data.groupby('chrom'):
        chrom_bin_results = calculate_pi_in_bins(group, samples, bin_size, step_size)
        bin_results.extend(chrom_bin_results)
    
    # Save the results to a file
    pi_bin_df = pd.DataFrame(bin_results)
    pi_bin_df.to_csv(output_file, sep='\t', index=False)
    print(f"Results saved to {output_file}")

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Sliding window pi calculation for specified samples.")
    parser.add_argument('-s', '--snp_file', required=True, help='Input SNP file in TSV format')
    parser.add_argument('-g', '--group_file', required=True, help='File with sample names to keep')
    parser.add_argument('-b', '--bin_size', type=int, default=100000, help='Window size (default: 100000)')
    parser.add_argument('-t', '--step_size', type=int, default=10000, help='Step size (default: 10000)')
    parser.add_argument('-o', '--output_file', required=True, help='Output file for results')

    args = parser.parse_args()

    main(args.snp_file, args.group_file, args.bin_size, args.step_size, args.output_file)
