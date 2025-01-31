import re

# Define file paths
dups_file = "phenotype_dups_to_investigate.txt"
gff_file = "GCF_943734745.1_idAnoDarlMG_H_01_genomic.gff"
output_file = "gene_products_output_2.txt"

def extract_unique_genes(annotation_line):
    """Extract unique gene names from a SNPEFF annotation line."""
    gene_pattern = re.compile(r"LOC\d+|XM_[0-9]+\.[0-9]+")
    return set(gene_pattern.findall(annotation_line))

def find_gene_product_in_gff(gene_name, gff_content):
    """Find the first occurrence of a gene name in the GFF file and extract product information."""
    gene_pattern = re.compile(fr"\b{gene_name}\b.*?product=([^\";]+)")
    for line in gff_content:
        match = gene_pattern.search(line)
        if match:
            return match.group(1)
    return "Product not found"

def process_files(dups_file, gff_file, output_file):
    """Process the dups file and the GFF file to extract gene names and products."""
    # Read GFF file content into a list
    with open(gff_file, 'r') as gff:
        gff_content = gff.readlines()

    # Open output file for writing
    with open(output_file, 'w') as output:
        
        # Read the duplications file
        with open(dups_file, 'r') as dups:
            lines = dups.readlines()

        chr_pos = None
        annotation_buffer = []

        for line in lines:
            line = line.strip()

            # Check if the line contains Chr and Pos
            chr_pos_match = re.match(r"(\S+)\s+(\d+)", line)
            if chr_pos_match:
                # Process buffered annotations for the previous Chr and Pos
                if chr_pos and annotation_buffer:
                    full_annotation = " ".join(annotation_buffer)
                    gene_names = extract_unique_genes(full_annotation)
                    if gene_names:
                        output.write(f"{chr_pos}\n")
                        for gene_name in gene_names:
                            product_info = find_gene_product_in_gff(gene_name, gff_content)
                            output.write(f"{gene_name}\t{product_info}\n")

                # Update Chr and Pos
                chr_pos = f"{chr_pos_match.group(1)}\t{chr_pos_match.group(2)}"
                annotation_buffer = []  # Clear buffer for new annotations
                continue

            # If it's an annotation line, add it to the buffer
            if chr_pos:
                annotation_buffer.append(line)

        # Process any remaining buffered annotations after the loop
        if chr_pos and annotation_buffer:
            full_annotation = " ".join(annotation_buffer)
            gene_names = extract_unique_genes(full_annotation)
            if gene_names:
                output.write(f"{chr_pos}\n")
                for gene_name in gene_names:
                    product_info = find_gene_product_in_gff(gene_name, gff_content)
                    output.write(f"{gene_name}\t{product_info}\n")

# Run the script
process_files(dups_file, gff_file, output_file)
print(f"Output written to {output_file}")

