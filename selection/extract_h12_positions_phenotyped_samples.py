import csv
import re
import argparse

def parse_gff_to_gene_ranges(gff_file):
    """Parse the GFF file to extract gene ranges, associated products, and descriptions."""
    gene_ranges = []
    gene_pattern = re.compile(r"\b(LOC[a-zA-Z0-9_.-]+)\b")
    product_pattern = re.compile(r"product=([^;]+)")
    description_pattern = re.compile(r"description=([^;]+)")
    
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith("#"):
                continue  # Skip comments
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] != "gene":
                continue  # Only process lines with "gene" type
            
            # Extract start, end, gene name, product, and description
            start = int(fields[3])
            end = int(fields[4])
            attributes = fields[8]
            
            gene_match = gene_pattern.search(attributes)
            product_match = product_pattern.search(attributes)
            description_match = description_pattern.search(attributes)
            
            gene_name = gene_match.group(1) if gene_match else "Unknown"
            product = product_match.group(1) if product_match else "Unknown"
            description = description_match.group(1) if description_match else "Unknown"
            
            gene_ranges.append((fields[0], start, end, gene_name, product, description))  # (Chromosome, Start, End, Gene, Product, Description)
    return gene_ranges

def find_gene_for_position(chromosome, position, gene_ranges):
    """Find the gene, product, and description for a given chromosome and position."""
    for gene_data in gene_ranges:
        gene_chr, start, end, gene_name, product, description = gene_data
        if gene_chr == chromosome and start <= position <= end:
            return gene_name, product, description
    return "No gene", "No product", "No description"

def process_csv_and_match_genes(input_csv, gff_file, output_csv):
    """Match genes, products, and descriptions to positions from the CSV file."""
    # Parse the GFF file
    gene_ranges = parse_gff_to_gene_ranges(gff_file)
    
    # Open output CSV for writing
    with open(output_csv, 'w', newline='') as out_csv:
        writer = csv.writer(out_csv)
        writer.writerow(["Chromosome", "Position", "H12", "Gene", "Product", "Description"])  # Write headers
        
        # Process input CSV
        with open(input_csv, 'r') as in_csv:
            reader = csv.DictReader(in_csv)
            for row in reader:
                chromosome = row["Chromosome"]
                position = int(float(row["Position"]))  # Convert to integer for range comparison
                h12 = row["H12"]
                
                # Find matching gene, product, and description
                gene_name, product, description = find_gene_for_position(chromosome, position, gene_ranges)
                writer.writerow([chromosome, position, h12, gene_name, product, description])

    print(f"Output written to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Match genes, products, and descriptions to positions in a CSV file.")
    parser.add_argument("input_csv", help="Path to the input CSV file.")
    parser.add_argument("gff_file", help="Path to the GFF file.")
    parser.add_argument("output_csv", help="Path to the output CSV file.")
    args = parser.parse_args()

    process_csv_and_match_genes(args.input_csv, args.gff_file, args.output_csv)
