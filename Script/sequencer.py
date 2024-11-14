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
    output = input("Enter the path for your output file: ")
    return gtf, bed, output

def parse_arguments():
    """Parse command-line arguments or use interactive prompts if not provided."""
    parser = argparse.ArgumentParser(description='Calculate STR densities in genes using GTF and BED files.')
    parser.add_argument('--gtf', type=str, help='Path to the GTF file with gene annotations.')
    parser.add_argument('--bed', type=str, help='Path to the BED file with STR regions.')
    parser.add_argument('--output', type=str, help='Path to the output CSV file.')
    args = parser.parse_args()

    if not args.gtf or not args.bed or not args.output:
        print("Not all command-line arguments provided, switching to interactive mode.")
        args.gtf, args.bed, args.output = interactive_prompt()

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

def get_gene_info(attribute):
    """Extract the gene name and gene ID from the attribute string using regex to handle different formatting."""
    gene_name_match = re.search(r'gene_name\s*"([^"]+)"', attribute)
    gene_id_match = re.search(r'gene_id\s*"([^"]+)"', attribute)
    gene_name = gene_name_match.group(1) if gene_name_match else "Unknown_Gene"
    gene_id = gene_id_match.group(1) if gene_id_match else "Unknown_ID"
    return gene_name, gene_id

def calculate_density(gtf_data, bed_data):
    """Calculate STR density for each gene in the GTF file."""
    results = []
    genes = gtf_data[gtf_data['feature'] == 'gene']
    for _, gene in genes.iterrows():
        gene_name, gene_id = get_gene_info(gene['attribute'])
        gene_interval = (gene['start'], gene['end'])
        
        str_data = bed_data[bed_data['chrom'] == gene['seqname']]
        str_data = str_data[(str_data["start"] > gene_interval[0]) & (str_data["end"] <= gene_interval[1])] # select only STRs that are located fully within gene
        
        total_overlap = (str_data["end"] - str_data["start"]).sum()
        gene_length = gene['end'] - gene['start'] + 1
        if gene_length > 0:
            density = total_overlap / gene_length
        else:
            density = 0

        results.append({"gene_id": gene_id, "gene_name": gene_name, "total_overlap": total_overlap, "density": density})
        print(f"Working on gene number {len(results)} (out of {len(genes)})", end='\r', file=sys.stderr)
    print()
    return pd.DataFrame(results)

def main():
    args = parse_arguments()
    gtf_data = load_data(args.gtf, "GTF")
    bed_data = load_data(args.bed, "BED")
    densities_df = calculate_density(gtf_data, bed_data)
    densities_df.to_csv(args.output, index=False)
    print(f"Results saved to {args.output}")

if __name__ == '__main__':
    main()
