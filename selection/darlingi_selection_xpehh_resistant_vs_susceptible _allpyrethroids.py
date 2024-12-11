# ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS

# %%
import os
import numpy as np
import allel
import zarr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import gffutils

# %% set wd
os.chdir('/mnt/storage11/sophie/darlingi/darlingi_database/phenotyped_colony')
os.getcwd()

# %%
# convert phased, filtered, VCF file to zarr file
# already converted to zarr
#allel.vcf_to_zarr('2019melasglobal_finalfiltered_gambiaealigned_phased.vcf.gz', '2019melasglobal_finalfiltered_gambiaealigned_phased.zarr', fields='*', overwrite=True)

# %%
callset = zarr.open('phenotyped_darlingi_filtered_phased.zarr', mode='r')
#callset.tree(expand=True)

# %%
## convert zarr file to genotype array
gt = allel.GenotypeDaskArray(callset['calldata/GT'])
print(gt.shape)

# %%
## import metadata
df_samples=pd.read_csv('/mnt/storage11/sophie/darlingi/darlingi_resistance_metadata.csv',sep=',',usecols=['sample','population','region','country','site','resistance_status'])
df_samples.head()
df_samples.groupby(by=['resistance_status']).count

# %%
## VCF is phased so we can convert genotype arrays made earlier to haplotype array

sample_ids = callset['samples'][:]

# %% Create arrays needed for pyrethroid resistant samples
# Get sample identifiers for pyrethroid_res samples from df_samples
pyrethroid_res_sample_ids = df_samples[
    (df_samples['resistance_status'] == 'alpha-cypermethrin-resistant') | 
    (df_samples['resistance_status'] == 'etofenprox-resistant')
]['sample'].values

# Find indices of these samples in the genotype array
pyrethroid_res_indices = np.array([np.where(sample_ids == id)[0][0] for id in pyrethroid_res_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", pyrethroid_res_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for pyrethroid_res samples using the indices
gt_pyrethroid_res_samples = gt.take(pyrethroid_res_indices, axis=1)

# %% Create arrays needed for pyrethroid susceptible samples
# Get sample identifiers for pyrethroid_sus samples from df_samples
pyrethroid_sus_sample_ids = df_samples[
    (df_samples['resistance_status'] == 'alpha-cypermethrin-susceptible') | 
    (df_samples['resistance_status'] == 'etofenprox-susceptible')
]['sample'].values

# Find indices of these samples in the genotype array
pyrethroid_sus_indices = np.array([np.where(sample_ids == id)[0][0] for id in pyrethroid_sus_sample_ids if id in sample_ids])
# Verify the indices are within the correct range
print("Max index:", pyrethroid_sus_indices.max(), "Sample array size:", len(sample_ids))
# Select genotypes for pyrethroid_sus samples using the indices
gt_pyrethroid_sus_samples = gt.take(pyrethroid_sus_indices, axis=1)

# %% these are from a phased VCF so we can convert the genotype arrays to haplotype arrays

h_array_pyrethroid_res = gt_pyrethroid_res_samples.to_haplotypes().compute()
h_array_pyrethroid_res

h_array_pyrethroid_sus = gt_pyrethroid_sus_samples.to_haplotypes().compute()
h_array_pyrethroid_sus

# %% we need variant positions
pos = callset['variants/POS'][:]
chrom = callset['variants/CHROM'][:]

# %% compute xpehh
# xpehh_raw = allel.xpehh(h_sus, h_res, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=20000, is_accessible=None, use_threads=True)
xpehh_raw = allel.xpehh(h_array_pyrethroid_sus, h_array_pyrethroid_res, pos, use_threads=True)
xpehh_raw

# %% look for where the biggest signal is
xpehh_hit_max = np.nanargmax(xpehh_raw)
xpehh_hit_max

# %% genomic position of top hit
pos[xpehh_hit_max]

# %% Plot the raw xp-ehh values

fig, ax = plt.subplots()
ax.hist(xpehh_raw[~np.isnan(xpehh_raw)], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% standardise xpehh-raw

allele_counts_array = gt.count_alleles(max_allele=3).compute()
xpehh_std = allel.standardize_by_allele_count(xpehh_raw, allele_counts_array[:, 1])

# %% plot standardised xp-ehh values

fig, ax = plt.subplots()
ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
ax.set_xlabel('Raw XP-EHH')
ax.set_ylabel('Frequency (no. variants)');

# %% Plot standardized xp-ehh

# define chromosome lengths and colours 
chromosome_lengths = {
    '2': 94951917,
    '3': 71270736,
    'anop_mito': 15395,
    'anop_X': 13401267
}

# Calculate cumulative offsets for each chromosome
cumulative_lengths = {}
cumulative_length = 0
for chrom, length in chromosome_lengths.items():
    cumulative_lengths[chrom] = cumulative_length
    cumulative_length += length

# %% Plot XP-EHH

# Set threshold
susceptible_threshold = 5
resistant_threshold = -5

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Ensure that pos, chrom, and xpehh_std are all numpy arrays to support advanced indexing
pos = np.array(callset['variants/POS'][:])
chrom = np.array(callset['variants/CHROM'][:])
xpehh_standardised_values = np.array(xpehh_std[0])

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Define colors for each chromosome (for illustration)
chromosome_colours = {
    '2': '#3d348b', '3': '#f7b801', 'anop_mito': '#f35b04', 'anop_X': '#119DA4'
}

# Create a list to hold the legend patches
legend_patches = []

# Filtered chromosomes list, assuming cumulative_lengths are defined for these
filtered_chroms = ['2', '3', 'anop_X', 'anop_mito']

# Iterate through each chromosome to plot its variants
for unique_chrom in filtered_chroms:
    chrom_mask = chrom == unique_chrom
    
    chrom_positions = pos[chrom_mask]
    chrom_xpehh_values = xpehh_standardised_values[chrom_mask]
    
    non_nan_mask = ~np.isnan(chrom_xpehh_values)
    chrom_positions_no_nan = chrom_positions[non_nan_mask]
    chrom_xpehh_values_no_nan = chrom_xpehh_values[non_nan_mask]
    
    adjusted_positions = chrom_positions_no_nan + cumulative_lengths.get(unique_chrom, 0)

    # Conditions for plotting
    solid_mask = (chrom_xpehh_values_no_nan >= susceptible_threshold) | (chrom_xpehh_values_no_nan <= resistant_threshold)
    faded_mask = ~solid_mask
    
    # Plot solid points for values above 5 or below -5
    ax.scatter(adjusted_positions[solid_mask], 
               chrom_xpehh_values_no_nan[solid_mask], 
               color=chromosome_colours[unique_chrom], alpha=1.0, s=10)
    
    # Plot faded points for other values
    ax.scatter(adjusted_positions[faded_mask], 
               chrom_xpehh_values_no_nan[faded_mask], 
               color=chromosome_colours[unique_chrom], alpha=0.1, s=10)
    
    # Add patch for the legend
    patch = mpatches.Patch(color=chromosome_colours[unique_chrom], label=unique_chrom)
    legend_patches.append(patch)

# Add significance threshold lines and legend
ax.axhline(y=susceptible_threshold, color='black', linestyle='--', label='Susceptible Threshold')
ax.axhline(y=resistant_threshold, color='black', linestyle='--', label='Resistant Threshold')
ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

# Set labels
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('XP-EHH')
plt.tight_layout()
plt.savefig('xpehh_plot_600ppi.png', dpi=600)  # Save at 600 PPI

# %% list all positions with xpehh value over or below a certain threshold

susceptible_threshold_mask = xpehh_standardised_values >= susceptible_threshold
resistant_threshold_mask = xpehh_standardised_values <= resistant_threshold

# %% Apply the mask to filter the data
sus_significant_chrom = chrom[susceptible_threshold_mask]
sus_significant_pos = pos[susceptible_threshold_mask]
sus_significant_xpehh = xpehh_standardised_values[susceptible_threshold_mask]

res_significant_chrom = chrom[resistant_threshold_mask]
res_significant_pos = pos[resistant_threshold_mask]
res_significant_xpehh = xpehh_standardised_values[resistant_threshold_mask]

# %% Combine the filtered data into a structured array
sus_significant_xpehh_data = np.column_stack((sus_significant_chrom, sus_significant_pos, sus_significant_xpehh))
res_significant_xpehh_data = np.column_stack((res_significant_chrom, res_significant_pos, res_significant_xpehh))

# %% Convert the structured array to a pandas DataFrame for easier handling
df_significant_sus_xpehh = pd.DataFrame(sus_significant_xpehh_data, columns=['Chromosome', 'Position', 'XPEHH'])
df_significant_res_xpehh = pd.DataFrame(res_significant_xpehh_data, columns = ['Chromosome', 'Position', 'XPEHH'])

# %% Save to csv
df_significant_sus_xpehh.to_csv(f'df_allpyrethroids_significant_sus_xpehh_susceptible_threshold_{susceptible_threshold}.csv', index=False)
df_significant_res_xpehh.to_csv(f'df_allpyrethroids_significant_res_xpehh_resistant_threshold_{resistant_threshold}.csv', index=False)


# %% Annotate the Susceptible XPEHH file

print("Using GFF file to bring in annotations for these positions")

# Parameters
input_file_name = f"df_allpyrethroids_significant_sus_xpehh_susceptible_threshold_{susceptible_threshold}.csv"
output_file_name = f"df_allpyrethroids_significant_sus_xpehh_susceptible_threshold_{susceptible_threshold}_annotated.csv"
gff_file = '/mnt/storage11/sophie/reference_genomes/An_darlingi_ncbi/GCF_943734745.1_idAnoDarlMG_H_01_genomic.gff'

# Function to find and format the GFF line(s) that overlap a given position
def find_overlapping_gff_lines(chromosome, position, gff_file):
    chromosome_mapping = {
        "NC_064873.1": "anop_X",
        "NC_064874.1": "2",
        "NC_064875.1": "3",
        "NC_064876.1": "anop_mito"
    }
    reverse_mapping = {v: k for k, v in chromosome_mapping.items()}
    gff_chromosome = reverse_mapping.get(chromosome)

    if not gff_chromosome:
        print(f"No mapping found for chromosome {chromosome}")
        return []

    overlapping_lines = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == "":
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            if parts[0] == gff_chromosome:
                start = int(parts[3])
                end = int(parts[4])
                if start <= position <= end:
                    print(f"Match found: {chromosome}:{position} in {parts[0]}:{start}-{end}")
                    annotation = parts[8]
                    formatted_line = f"{parts[0]}:{start}-{end} | {annotation}"
                    overlapping_lines.append(formatted_line)
    if not overlapping_lines:
        print(f"No overlap found for {chromosome}:{position}")
    return overlapping_lines

# Writing output file
with open(output_file_name, "w") as outfile:
    # Write the header line
    outfile.write("Chromosome\tPosition\tXPEHH Value\tGff_Annotation\n")

    with open(input_file_name, "r") as infile:
        next(infile)  # Skip the header line
        for line in infile:
            # Split the line using comma as the delimiter
            parts = line.strip().split(",")
            if len(parts) < 3:  # Ensure the line has at least three columns
                print(f"Skipping malformed line: {line.strip()}")
                continue

            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]

            # Find overlapping annotations in the GFF file
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            if overlapping_gff_lines:
                gff_annotation = "; ".join(overlapping_gff_lines)
            else:
                gff_annotation = "no overlapping line found in the GFF file for this position"

            # Write the annotated result to the output file
            outfile.write(f"{chromosome}\t{position}\t{XPEHH}\t{gff_annotation}\n")

print(f"XP-EHH significant values for susceptible samples identified and GFF annotations written to: {output_file_name}")

# %% Annotate the Resistant XPEHH file
print("Using GFF file to annotate resistant XPEHH positions")

input_file_name = f"df_allpyrethroids_significant_res_xpehh_resistant_threshold_{resistant_threshold}.csv"
output_file_name = f"df_allpyrethroids_significant_res_xpehh_resistant_threshold_{resistant_threshold}_annotated.csv"
gff_file = '/mnt/storage11/sophie/reference_genomes/An_darlingi_ncbi/GCF_943734745.1_idAnoDarlMG_H_01_genomic.gff'

# Function to find and format the GFF line(s) that overlap a given position
def find_overlapping_gff_lines(chromosome, position, gff_file):
    chromosome_mapping = {
        "NC_064873.1": "anop_X",
        "NC_064874.1": "2",
        "NC_064875.1": "3",
        "NC_064876.1": "anop_mito"
    }
    reverse_mapping = {v: k for k, v in chromosome_mapping.items()}
    gff_chromosome = reverse_mapping.get(chromosome)

    if not gff_chromosome:
        print(f"No mapping found for chromosome {chromosome}")
        return []

    overlapping_lines = []
    with open(gff_file, 'r') as gff:
        for line in gff:
            if line.startswith('#') or line.strip() == "":
                continue
            
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            if parts[0] == gff_chromosome:
                start = int(parts[3])
                end = int(parts[4])
                if start <= position <= end:
                    print(f"Match found: {chromosome}:{position} in {parts[0]}:{start}-{end}")
                    annotation = parts[8]
                    formatted_line = f"{parts[0]}:{start}-{end} | {annotation}"
                    overlapping_lines.append(formatted_line)
    if not overlapping_lines:
        print(f"No overlap found for {chromosome}:{position}")
    return overlapping_lines

# Writing output file
with open(output_file_name, "w") as outfile:
    # Write the header line
    outfile.write("Chromosome\tPosition\tXPEHH Value\tGff_Annotation\n")

    with open(input_file_name, "r") as infile:
        next(infile)  # Skip the header line
        for line in infile:
            # Split the line using comma as the delimiter
            parts = line.strip().split(",")
            if len(parts) < 3:  # Ensure the line has at least three columns
                print(f"Skipping malformed line: {line.strip()}")
                continue

            chromosome, position, XPEHH = parts[0], int(parts[1]), parts[2]

            # Find overlapping annotations in the GFF file
            overlapping_gff_lines = find_overlapping_gff_lines(chromosome, position, gff_file)
            
            if overlapping_gff_lines:
                gff_annotation = "; ".join(overlapping_gff_lines)
            else:
                gff_annotation = "no overlapping line found in the GFF file for this position"

            # Write the annotated result to the output file
            outfile.write(f"{chromosome}\t{position}\t{XPEHH}\t{gff_annotation}\n")

print(f"XP-EHH significant values for resistant samples identified and GFF annotations written to: {output_file_name}")

# %%