import pandas as pd
import matplotlib.pyplot as plt

# Load the Tajima's D data
file_path = "rondonia_samples_darlingi_filtered_phased.Tajima.D"
data = pd.read_csv(file_path, sep="\t")

# Ensure the column names are correct
data.columns = ["CHROM", "BIN_START", "N_SNPS", "TajimaD"]

# Filter out rows with "nan" in TajimaD
data = data.dropna(subset=["TajimaD"])

# Convert BIN_START to a genomic coordinate for plotting
data["GENOMIC_POSITION"] = data["BIN_START"]

# Get the unique chromosomes
chromosomes = data["CHROM"].unique()

# Create a plot for each chromosome
for chrom in chromosomes:
    chrom_data = data[data["CHROM"] == chrom]
    
    plt.figure(figsize=(12, 6))
    plt.plot(
        chrom_data["GENOMIC_POSITION"],
        chrom_data["TajimaD"],
        alpha=0.8,
        marker="o",
        linestyle="-",
    )
    
    # Add labels, title, and horizontal line
    plt.title(f"Rondonia State Mosquitoes - Tajima's D Across Chromosome {chrom}")
    plt.xlabel("Genomic Position (bp)")
    plt.ylabel("Tajima's D")
    plt.axhline(y=0, color='gray', linestyle='--', linewidth=0.8, label="Neutral Expectation")
    plt.legend(["Tajima's D"], loc='upper left')
    
    # Save the plot
    plt.tight_layout()
    plt.savefig(f"rondonia_mosquitoes_tajimas_d_{chrom}.tiff", dpi=300)
    plt.show()
