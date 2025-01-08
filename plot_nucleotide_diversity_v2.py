# Run this script using python plot_nucleotide_diversity.py <path_to_windowed.pi_file>

import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Plot nucleotide diversity from a .windowed.pi file.")
parser.add_argument("file", help="Path to the .windowed.pi file to be processed")
args = parser.parse_args()

file_path = args.file

# Extract gene name from the file name
gene_name = file_path.split('_')[0]

# Load the data
data = pd.read_csv(file_path, sep="\t", comment='#', header=0)

# Rename columns for consistency
data.columns = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI"]

# Debugging: Print the first few rows to verify data integrity
print("Loaded Data Preview:")
print(data.head())

# Check for unique chromosomes and bins
print("Unique Chromosomes:")
print(data["CHROM"].unique())
print(f"Total bins: {len(data)}")

# Calculate the genomic position as the center of each bin
data["GENOMIC_POSITION"] = data["BIN_START"] + (data["BIN_END"] - data["BIN_START"]) // 2

# Debugging: Ensure all bins have GENOMIC_POSITION calculated
print("Data with GENOMIC_POSITION:")
print(data[["CHROM", "BIN_START", "BIN_END", "GENOMIC_POSITION", "PI"]].head())

# Create the plot
plt.figure(figsize=(12, 6))

# Plot nucleotide diversity for each chromosome
for chrom in data["CHROM"].unique():
    chrom_data = data[data["CHROM"] == chrom]
    plt.plot(
        chrom_data["GENOMIC_POSITION"],
        chrom_data["PI"],
        label=gene_name,
        alpha=0.8,
        marker="o",
        linestyle="-",
    )

# Add labels and title
plt.title(f"Nucleotide Diversity for {gene_name}")
plt.xlabel("Genomic Position (bp)")
plt.ylabel("Nucleotide Diversity (PI)")
plt.legend(title="Gene", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save or show the plot
output_file = f"{gene_name}_nucleotide_diversity_plot.png"
plt.savefig(output_file, dpi=300)
plt.show()

print(f"Plot saved as {output_file}")
