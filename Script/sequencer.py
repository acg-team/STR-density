import argparse
import pandas as pd
import sys
import re

def interactive_prompt():
    """Provides an interactive command line interface to gather inputs from the user."""
    print("Welcome to the STR Density Calculator.")
    print("This program calculates STR density for each gene in the provided GTF and BED files.")
    gtf = input("Enter the path to your GTF file: ")
    bed = input("Enter the path to your BED file: ")
    return gtf, bed

def parse_arguments():
    """Parse command-line arguments or use interactive prompts if not provided."""
    parser = argparse.ArgumentParser(description='Calculate STR densities in genes using GTF and BED files.')
    parser.add_argument('--gtf', type=str, help='Path to the GTF file with gene annotations.')
    parser.add_argument('--bed', type=str, help='Path to the BED file with STR regions.')
    args = parser.parse_args()

    if not args.gtf or not args.bed:
        print("Not all command-line arguments provided, switching to interactive mode.")
        args.gtf, args.bed = interactive_prompt()

    return args

def load_data(filepath, file_type):
    """Load data files safely with error handling."""
    try:
        if file_type == "BED":
            return pd.read_csv(filepath, sep='\t', header=None, names=['chrom', 'start', 'end', 'name', 'sequence'])
        elif file_type == "GTF":
            return pd.read_csv(filepath, sep='\t', header=None, comment='#', names=['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute'])
    except FileNotFoundError:
        sys.exit(f"Error: The file {filepath} does not exist.")
    except pd.errors.EmptyDataError:
        sys.exit(f"Error: The file {filepath} is empty or invalid.")
    except Exception as e:
        sys.exit(f"Error loading file {filepath}: {str(e)}")

def get_gene_name(attribute):
    """Extract the gene name from the attribute string using regex to handle different formatting."""
    match = re.search(r'gene_name\s*"([^"]+)"', attribute)
    if match:
        return match.group(1)
    return "Unknown_Gene"  # Return a placeholder if the gene name isn't found

def calculate_density(gtf_data, bed_data):
    """Calculate STR density for each gene in the GTF file."""
    results = []
    genes = gtf_data[gtf_data['feature'] == 'gene']
    for _, gene in genes.iterrows():
        gene_name = get_gene_name(gene['attribute'])
        gene_interval = (gene['start'], gene['end'])
        chromosome = gene['seqname']
        str_data = bed_data[bed_data['chrom'] == chromosome]
        total_overlap = 0
        gene_length = gene['end'] - gene['start'] + 1
        for _, str_row in str_data.iterrows():
            str_interval = (str_row['start'], str_row['end'])
            overlap = max(0, min(gene_interval[1], str_interval[1]) - max(gene_interval[0], str_interval[0]))
            total_overlap += overlap
        if gene_length > 0:
            density = total_overlap / gene_length
        else:
            density = 0
        results.append((gene_name, total_overlap, density))
    return results

def main():
    args = parse_arguments()
    gtf_data = load_data(args.gtf, "GTF")
    bed_data = load_data(args.bed, "BED")
    densities = calculate_density(gtf_data, bed_data)
    for gene_name, total_overlap, density in densities:
        print(f"Gene: {gene_name}, Total Overlap: {total_overlap}, STR Density: {density}")

if __name__ == '__main__':
    main()
