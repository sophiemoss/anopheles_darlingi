# ########### Cross-population extended haplotype homozygosity (XPEHH) ###########

## Compute the unstandardized cross-population extended haplotype homozygosity score (XPEHH) for each variant.
## allel.xpehh(h1, h2, pos, map_pos=None, min_ehh=0.05, include_edges=False, gap_scale=20000, max_gap=200000, is_accessible=None, use_threads=True)
# create h1 and h2, selecting all variants instead of segregating variants only, which is what we did in iHS
# xpehh is calculated for every single variant

# %%
import os
import numpy as np
import allel
import zarr
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import sys
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
df_samples=pd.read_csv('/mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/darlingi_resistance_metadata.csv',sep=',',usecols=['sample','population','region','country','site','resistance_status', 'pyrethroid_resistance_status'])
df_samples.head()
df_samples.groupby(by=['pyrethroid_resistance_status']).count

# %%
## VCF is phased so we can convert genotype arrays made earlier to haplotype array

sample_ids = callset['samples'][:]

# %% Create arrays needed for pyrethroid resistant samples
# Get sample identifiers for pyrethroid_res samples from df_samples
pyrethroid_res_sample_ids = df_samples[
    (df_samples['pyrethroid_resistance_status'] == 'resistant')
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
    (df_samples['pyrethroid_resistance_status'] == 'susceptible')
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

#allele_counts_array = gt.count_alleles(max_allele=3).compute()
#xpehh_std = allel.standardize_by_allele_count(xpehh_raw, allele_counts_array[:, 1])

# %% plot standardised xp-ehh values

# fig, ax = plt.subplots()
# ax.hist(xpehh_std[0][~np.isnan(xpehh_std[0])], bins=20)
# ax.set_xlabel('Raw XP-EHH')
# ax.set_ylabel('Frequency (no. variants)');


#%% Plot raw XP-EHH, define chromosome lengths and colours 
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

# %% Plot raw XP-EHH

# Set thresholds
susceptible_threshold = 5
resistant_threshold = -5

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Ensure that pos, chrom, and xpehh_values are all numpy arrays
pos = np.array(callset['variants/POS'][:])
chrom = np.array(callset['variants/CHROM'][:])
xpehh_values = np.array(xpehh_raw)

# Define colors for each chromosome
chromosome_colours = {
    '2': '#3d348b', '3': '#f7b801', 'anop_mito': '#f35b04', 'anop_X': '#119DA4'
}

# Create a list to hold the legend patches
legend_patches = []

# List of chromosomes to plot
filtered_chroms = ['2', '3', 'anop_X', 'anop_mito']

# Iterate through each chromosome to plot its variants
for unique_chrom in filtered_chroms:
    chrom_mask = chrom == unique_chrom
    
    chrom_positions = pos[chrom_mask]  # Keep all positions
    chrom_xpehh_values = xpehh_values[chrom_mask]  # Keep all XP-EHH values
    
    # Adjust positions to be cumulative (if necessary)
    adjusted_positions = chrom_positions + cumulative_lengths.get(unique_chrom, 0)

    # Conditions for plotting
    solid_mask = (~np.isnan(chrom_xpehh_values)) & ((chrom_xpehh_values >= susceptible_threshold) | (chrom_xpehh_values <= resistant_threshold))
    faded_mask = (~np.isnan(chrom_xpehh_values)) & ~solid_mask
    nan_mask = np.isnan(chrom_xpehh_values)  # Find NaN values

    # Plot solid points (highlighted XP-EHH values)
    ax.scatter(adjusted_positions[solid_mask], 
               chrom_xpehh_values[solid_mask], 
               color=chromosome_colours[unique_chrom], alpha=1.0, s=10)

    # Plot faded points (regular XP-EHH values)
    ax.scatter(adjusted_positions[faded_mask], 
               chrom_xpehh_values[faded_mask], 
               color=chromosome_colours[unique_chrom], alpha=0.1, s=10)

    # Plot NaN values as empty/transparent points (preserving positions)
    ax.scatter(adjusted_positions[nan_mask], 
               np.zeros_like(adjusted_positions[nan_mask]), 
               color='white', alpha=0.0, s=1)

    # Add chromosome to legend
    patch = mpatches.Patch(color=chromosome_colours[unique_chrom], label=unique_chrom)
    legend_patches.append(patch)

# Add threshold lines and legend
ax.axhline(y=susceptible_threshold, color='black', linestyle='--', label='Susceptible Threshold')
ax.axhline(y=resistant_threshold, color='black', linestyle='--', label='Resistant Threshold')
ax.legend(handles=legend_patches, title='Chromosome', bbox_to_anchor=(1.05, 1), loc='upper left')

# Set labels
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('XP-EHH')
plt.tight_layout()
plt.savefig('xpehh_plot_600ppi.png', dpi=600)  # Save at 600 PPI
plt.show()

####################### PERMUTATIONS ######################
# %% Caluclate permutations where the resistant and susceptible labels are switched between different sample IDs

permuted_xpehh_values = []
for i in range(5):
    # Get the indices from df_samples
    df_res_sus_samples = df_samples[df_samples.pyrethroid_resistance_status.isin(['resistant', 'susceptible'])]
    indices = df_res_sus_samples.index.tolist()
    # Shuffle the indices
    np.random.shuffle(indices)
    # Split the indices into two groups
    half_length = len(indices) // 2
    pop1_indices = indices[:half_length]
    pop2_indices = indices[half_length:]
   
    # make new genotype arrays from permutations
    gt_pop1 = gt.take(pop1_indices, axis =1)
    gt_pop2 = gt.take(pop2_indices, axis =1)

    # create the new haplotype arrays
    h_pop1 = gt_pop1.to_haplotypes().compute()
    h_pop2 = gt_pop2.to_haplotypes().compute()

    # check that pos and genotype are the same size
    if len(pos)==len(gt_pop1):
        print("Length of positions and genotypes in the genotype array are the same, script continuing")
    else:
        print("Something is wrong with the genotype_all array as the length of pos_all and genotype_all are different. Stopping script.")
        sys.exit()  # This will stop the script. If you want the script to continue anyway, # out this line
    
    # Calcalate and plot xpehh for iteration
    xpehh_raw_iterated = allel.xpehh(h_pop1, h_pop2, pos, use_threads=True)
    print(f"Calculated xpehh using allel.xpehh for permutation {i}")

    # Store the xp-ehh values for this iteration
    permuted_xpehh_values.append((xpehh_raw_iterated))

# Notify of finishing permutations
print("Permutations calculated")

#%% Save the permuted Xpehh values as a dataframe and as a csv, so that you do not need to calcualte them again if taking a break in analysis
permuted_xpehh_values_df = pd.DataFrame(permuted_xpehh_values)
# Save permuted_xpehh_values as a csv so that you do not need to calculate again if taking a break in analysis
permuted_xpehh_values_df.to_csv(f'permuted_xpehh_values.csv', index=False)

#%% Plot the iterated XP-EHH values

# Set up the plot
fig, ax = plt.subplots(figsize=(10, 6))

# Ensure that pos, chrom, and xpehh_std are all numpy arrays to support advanced indexing
pos = np.array(callset['variants/POS'][:])
chrom = np.array(callset['variants/CHROM'][:])
permuted_xpehh_values_all = np.array(permuted_xpehh_values_df)

# Define colors for each chromosome (for illustration)
chromosome_colours = {
    '2': 'grey', '3': 'lightgrey', 'anop_mito': 'lightgrey', 'anop_X': 'grey'
}

# Create a list to hold the legend patches
legend_patches = []

# Filtered chromosomes list, assuming cumulative_lengths are defined for these
filtered_chroms = ['2', '3', 'anop_X', 'anop_mito']

# Iterate through each chromosome to plot its variants
# Iterate over each of the permuted arrays separately to plot all of the points
for i, permuted_xpehh_values in enumerate(permuted_xpehh_values_all):
    for unique_chrom in filtered_chroms:
        chrom_mask = chrom == unique_chrom
        
        chrom_positions = pos[chrom_mask]  # Positions remain fixed
        chrom_xpehh_values = permuted_xpehh_values[chrom_mask]  # Extract values from current sub-array
        
        non_nan_mask = ~np.isnan(chrom_xpehh_values)
        chrom_positions_no_nan = chrom_positions[non_nan_mask]
        chrom_xpehh_values_no_nan = chrom_xpehh_values[non_nan_mask]
        
        adjusted_positions = chrom_positions_no_nan + cumulative_lengths.get(unique_chrom, 0)

        # Apply threshold filters
        solid_mask = (chrom_xpehh_values_no_nan >= susceptible_threshold) | (chrom_xpehh_values_no_nan <= resistant_threshold)
        faded_mask = ~solid_mask
        
        # Scatter plot for each permutation separately
        ax.scatter(adjusted_positions[solid_mask], chrom_xpehh_values_no_nan[solid_mask], 
                   color=chromosome_colours[unique_chrom], alpha=1.0, s=10, label=f"Perm {i+1}" if i == 0 else "")
        ax.scatter(adjusted_positions[faded_mask], chrom_xpehh_values_no_nan[faded_mask], 
                   color=chromosome_colours[unique_chrom], alpha=0.1, s=10)

# Set labels
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('XP-EHH value')
plt.tight_layout()
plt.savefig('permuted_xpehh_plot_600ppi.png', dpi=600)  # Save at 600 PPI

#%% ### NEED TO ENSURE THAT SIGNIFICANCE IS CORRECT FOR POSITIVES AND NEGATIVES. THIS DOES NOT WORK AS IT JUST GIVES POSITIVES.

# Convert list of arrays into a 2D NumPy array (shape: num_permutations x num_variants)
permuted_xpehh_values_all_array = np.vstack(permuted_xpehh_values_all)

# Compute the 99th percentile for each variant (ignoring NaNs)
percentile_99_values = np.nanpercentile(permuted_xpehh_values_all_array, 99, axis=0)

# percentile_99_values now contains the 99th percentile XP-EHH for each variant
# you can look at it using np.nanmax(percentile_99_values). Need to use np.nanmax to ignore nans, otherwise np.max() would return nan if NaNs are present.

# Now compare the actual standardized XpEHH value to the computed 99th percentile

# Identify significant positions
significant_positions = pos[xpehh_raw > percentile_99_values]

# Extract corresponding XP-EHH values
significant_xpehh_values = xpehh_raw[xpehh_raw > percentile_99_values]

# Create DataFrame of significant results
significant_variants_df = pd.DataFrame({
    'Position': significant_positions,
    'XP-EHH': significant_xpehh_values,
    'Threshold_99th': percentile_99_values[xpehh_raw > percentile_99_values]
})

# Save results
significant_variants_df.to_csv('significant_xpehh_variants.csv', index=False)

# %% Visualise significant variants

# Plot the XP-EHH values with the 99th percentile threshold
fig, ax = plt.subplots(figsize=(10, 6))

# Scatter plot of all XP-EHH values
ax.scatter(pos, xpehh_std[0], alpha=0.3, label="All XP-EHH values", s=5)

# Highlight significant variants
ax.scatter(significant_positions, significant_xpehh_values, color='red', label="Significant (99th percentile)", s=10)

# Add the 99th percentile threshold line
ax.axhline(y=np.nanmax(percentile_99th_values), color='black', linestyle='--', label="99th Percentile Threshold")

# Labels and legend
ax.set_xlabel('Genomic Position (bp)')
ax.set_ylabel('Standardized XP-EHH')
ax.legend()
plt.title('XP-EHH with Significant Variants')

# Save plot
plt.savefig('xpehh_significant_variants_plot.png', dpi=600)
plt.show()
















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



























### Annotations ###
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
