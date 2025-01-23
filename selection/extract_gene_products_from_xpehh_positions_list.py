import re
import csv
import argparse
import os

def extract_unique_genes(annotation_line, gene_pattern):
    """Extract unique gene names from an annotation line."""
    return set(gene_pattern.findall(annotation_line))

def find_gene_product_in_gff(gene_name, gff_content, product_pattern):
    """Find the first occurrence of a gene name in the GFF file and extract product information."""
    for line in gff_content:
        match = product_pattern.search(line)
        if match:
            return match.group(1)
    return "Product not found"

def process_csv_file(input_csv, gff_file, output_file):
    """Process the input CSV file and extract gene products from the GFF file."""
    # Compile regex patterns
    gene_pattern = re.compile(r"LOC[a-zA-Z0-9_.-]+")
    product_pattern = re.compile(rf"\b(LOC[a-zA-Z0-9_.-]+)\b.*?product=([^\";]+)")

    # Validate GFF file
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF file not found: {gff_file}")

    # Read GFF content into memory
    with open(gff_file, 'r') as gff:
        gff_content = gff.readlines()

    # Open output file for writing
    with open(output_file, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Chromosome", "Position", "Gene", "Product"])  # Write headers

        # Process input CSV
        with open(input_csv, 'r') as in_csv:
            reader = csv.DictReader(in_csv)
            for row in reader:
                chromosome = row["Chromosome"]
                position = row["Position"]
                annotation = row["Gff_Annotation"]

                # Extract unique gene names
                gene_names = extract_unique_genes(annotation, gene_pattern)
                for gene_name in gene_names:
                    # Find the first occurrence of the gene in the GFF file
                    product_info = find_gene_product_in_gff(gene_name, gff_content, product_pattern)
                    writer.writerow([chromosome, position, gene_name, product_info])

    print(f"Output written to {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process CSV file and extract gene products.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    parser.add_argument("gff_file", help="Path to the GFF file.")
    parser.add_argument("output_file", help="Path to the output CSV file.")
    args = parser.parse_args()

    process_csv_file(args.input_csv, args.gff_file, args.output_file)
