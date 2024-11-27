import argparse
import pandas as pd
import sys
import re
import json
from pandas import IntervalIndex

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
    """Optimized calculation of STR density for genes, exons, and introns."""


    bed_data['str_type'] = bed_data['sequence'].apply(lambda seq: f"{len(seq)}-mer")

    # Pre-filter gene and exon data
    genes = gtf_data[gtf_data['feature'] == 'gene'].copy()
    exons = gtf_data[gtf_data['feature'] == 'exon'].copy()
    results = []

    for chrom in gtf_data['seqname'].unique():
        chrom_strs = bed_data[bed_data['chrom'] == chrom]
        chrom_genes = genes[genes['seqname'] == chrom]
        chrom_exons = exons[exons['seqname'] == chrom]

        for idx, gene in chrom_genes.iterrows():
            gene_name, gene_id = get_gene_info(gene['attribute'])
            gene_start, gene_end = gene['start'], gene['end']

            # Gene STR density
            gene_strs = chrom_strs[(chrom_strs['start'] >= gene_start) & (chrom_strs['end'] <= gene_end)]
            total_gene_overlap = (gene_strs['end'] - gene_strs['start']).sum()
            gene_length = gene_end - gene_start + 1
            gene_density = total_gene_overlap / gene_length if gene_length > 0 else 0

            # STR type aggregation
            str_type_counts = gene_strs['str_type'].value_counts().to_dict()

            # Exon STR density
            overlapping_exons = chrom_exons[
                (chrom_exons['start'] <= gene_end) & (chrom_exons['end'] >= gene_start)
            ]
            exon_intervals_gene = IntervalIndex.from_arrays(overlapping_exons['start'], overlapping_exons['end'], closed='both')
            total_exon_overlap = sum(
                (chrom_strs[(chrom_strs['start'] >= exon.left) & (chrom_strs['end'] <= exon.right)]['end'] -
                 chrom_strs[(chrom_strs['start'] >= exon.left) & (chrom_strs['end'] <= exon.right)]['start']).sum()
                for exon in exon_intervals_gene
            )
            total_exon_length = sum(exon.length for exon in exon_intervals_gene)
            exon_density = total_exon_overlap / total_exon_length if total_exon_length > 0 else 0

            # Intron STR density
            intron_intervals = [
                (exon_intervals_gene[i].right + 1, exon_intervals_gene[i + 1].left - 1)
                for i in range(len(exon_intervals_gene) - 1)
                if exon_intervals_gene[i].right + 1 < exon_intervals_gene[i + 1].left
            ]
            total_intron_overlap = sum(
                (chrom_strs[(chrom_strs['start'] >= start) & (chrom_strs['end'] <= end)]['end'] -
                 chrom_strs[(chrom_strs['start'] >= start) & (chrom_strs['end'] <= end)]['start']).sum()
                for start, end in intron_intervals
            )
            total_intron_length = sum(end - start + 1 for start, end in intron_intervals)
            intron_density = total_intron_overlap / total_intron_length if total_intron_length > 0 else 0

            # Append results
            results.append({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "total_gene_overlap": total_gene_overlap,
                "gene_density": round(gene_density, 6),
                "total_exon_overlap": total_exon_overlap,
                "exon_density": round(exon_density, 6),
                "total_intron_overlap": total_intron_overlap,
                "intron_density": round(intron_density, 6),
                "str_types": json.dumps(str_type_counts)
            })
        print(f"Processed chromosome {chrom}", file=sys.stderr)

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

