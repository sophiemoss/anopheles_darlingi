import re
import csv
import argparse
import os

def parse_gff_to_gene_product_map(gff_file):
    """Parse the GFF file and create a mapping of gene names to product information."""
    gene_product_map = {}
    product_pattern = re.compile(r"\b(LOC[a-zA-Z0-9_.-]+)\b.*?product=([^\";]+)")

    with open(gff_file, 'r') as gff:
        for line in gff:
            match = product_pattern.search(line)
            if match:
                gene_name = match.group(1)
                product_info = match.group(2)
                gene_product_map[gene_name] = product_info
    return gene_product_map

def extract_unique_genes(annotation_line, gene_pattern):
    """Extract unique gene names from an annotation line."""
    return set(gene_pattern.findall(annotation_line))

def process_tsv_file(input_tsv, gff_file, output_file):
    """Process the input TSV file and extract gene products from the GFF file."""
    # Compile regex patterns
    gene_pattern = re.compile(r"LOC[a-zA-Z0-9_.-]+")

    # Parse GFF file to create a gene-to-product map
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")
    gene_product_map = parse_gff_to_gene_product_map(gff_file)

    # Open output file for writing
    with open(output_file, 'w', newline='') as out_tsv:
        writer = csv.writer(out_tsv, delimiter='\t')
        writer.writerow(["Chromosome", "Position", "Gene", "Product"])  # Write headers

        # Process input TSV
        with open(input_tsv, 'r') as in_tsv:
            reader = csv.DictReader(in_tsv, delimiter='\t')  # Specify tab delimiter
            for row in reader:
                chromosome = row["Chromosome"]
                position = row["Position"]
                annotation = row["Gff_Annotation"]

                # Extract unique gene names
                gene_names = extract_unique_genes(annotation, gene_pattern)
                for gene_name in gene_names:
                    # Retrieve product info from the map
                    product_info = gene_product_map.get(gene_name, "Product not found")
                    writer.writerow([chromosome, position, gene_name, product_info])

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process TSV file and extract gene products.")
    parser.add_argument("input_tsv", help="Path to the input TSV file.")
    parser.add_argument("gff_file", help="Path to the GFF file.")
    parser.add_argument("output_file", help="Path to the output TSV file.")
    args = parser.parse_args()

    process_tsv_file(args.input_tsv, args.gff_file, args.output_file)
