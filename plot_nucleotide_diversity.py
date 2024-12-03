import pandas as pd
import matplotlib.pyplot as plt

# Load the data
# Replace 'colony_old_nuc_div_window_10kb.windowed.pi' with your file name
file_path = "colony_old_nuc_div_window_10kb.windowed.pi"
data = pd.read_csv(file_path, sep="\t")

# Ensure the data is in the correct format
# Rename columns for easier access (if needed)
data.columns = ["CHROM", "BIN_START", "BIN_END", "N_VARIANTS", "PI"]

# Convert BIN_START to a genomic coordinate for plotting
data["GENOMIC_POSITION"] = data["BIN_START"] + (data["BIN_END"] - data["BIN_START"]) // 2

# Plot nucleotide diversity (PI)
plt.figure(figsize=(12, 6))
for chrom in data["CHROM"].unique():
    chrom_data = data[data["CHROM"] == chrom]
    plt.plot(
        chrom_data["GENOMIC_POSITION"],
        chrom_data["PI"],
        label=chrom,
        alpha=0.8,
        marker="o",
        linestyle="-",
    )

# Add labels and title
plt.title("Nucleotide Diversity Across the Genome")
plt.xlabel("Genomic Position (bp)")
plt.ylabel("Nucleotide Diversity (PI)")
plt.legend(title="Chromosome", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

# Save or show the plot
plt.savefig("nucleotide_diversity_plot.png", dpi=300)
plt.show()
