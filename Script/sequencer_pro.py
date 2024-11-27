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
    """Calculate STR density for genes, exons, and introns, including STR types."""
    # Precompute STR types in BED data
    bed_data['str_type'] = bed_data['sequence'].apply(lambda seq: f"{len(seq)}-mer")

    # Filter genes and exons from GTF
    genes = gtf_data[gtf_data['feature'] == 'gene'].copy()
    exons = gtf_data[gtf_data['feature'] == 'exon'].copy()

    # Initialize results list
    results = []

    # Process by chromosome to reduce filtering overhead
    for chrom in gtf_data['seqname'].unique():
        chrom_strs = bed_data[bed_data['chrom'] == chrom]
        chrom_genes = genes[genes['seqname'] == chrom]
        chrom_exons = exons[exons['seqname'] == chrom]

        for _, gene in chrom_genes.iterrows():
            gene_name, gene_id = get_gene_info(gene['attribute'])
            gene_start, gene_end = gene['start'], gene['end']

            # Filter STRs within the gene region
            gene_strs = chrom_strs[(chrom_strs['start'] >= gene_start) & (chrom_strs['end'] <= gene_end)]
            total_gene_overlap = (gene_strs['end'] - gene_strs['start']).sum()
            gene_length = gene_end - gene_start + 1
            gene_density = total_gene_overlap / gene_length if gene_length > 0 else 0

            # Aggregate STR types
            str_type_counts = gene_strs['str_type'].value_counts().to_dict()

            # Filter exons within the gene region
            gene_exons = chrom_exons[(chrom_exons['start'] >= gene_start) & (chrom_exons['end'] <= gene_end)]
            exon_overlaps = chrom_strs.apply(
                lambda row: (gene_exons['start'] <= row['end']) & (gene_exons['end'] >= row['start']),
                axis=1
            )
            total_exon_overlap = exon_overlaps.sum()
            total_exon_length = (gene_exons['end'] - gene_exons['start'] + 1).sum()
            exon_density = total_exon_overlap / total_exon_length if total_exon_length > 0 else 0

            # Calculate intron regions
            intron_intervals = [
                (gene_exons.iloc[i]['end'] + 1, gene_exons.iloc[i + 1]['start'] - 1)
                for i in range(len(gene_exons) - 1)
                if gene_exons.iloc[i]['end'] + 1 < gene_exons.iloc[i + 1]['start']
            ]
            total_intron_overlap = sum(
                chrom_strs[(chrom_strs['start'] >= start) & (chrom_strs['end'] <= end)]['end'] -
                chrom_strs[(chrom_strs['start'] >= start) & (chrom_strs['end'] <= end)]['start']
                for start, end in intron_intervals
            )
            total_intron_length = sum(end - start + 1 for start, end in intron_intervals)
            intron_density = total_intron_overlap / total_intron_length if total_intron_length > 0 else 0

            # Append results
            results.append({
                "gene_id": gene_id,
                "gene_name": gene_name,
                "total_gene_overlap": total_gene_overlap,
                "gene_density": gene_density,
                "total_exon_overlap": total_exon_overlap,
                "exon_density": exon_density,
                "total_intron_overlap": total_intron_overlap,
                "intron_density": intron_density,
                "str_types": str_type_counts
            })
        print(f"Processed chromosome {chrom}", file=sys.stderr)

    return pd.DataFrame(results)


def main():
    args = parse_arguments()
    gtf_data = load_data(args.gtf, "GTF")
    bed_data = load_data(args.bed, "BED")
    densities = calculate_density(gtf_data, bed_data)
    for gene_name, total_overlap, density in densities:
        print(f"Gene: {gene_name}, Total Overlap: {total_overlap}, STR Density: {density}")

if __name__ == '__main__':
    main()
