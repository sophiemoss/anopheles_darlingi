#!/bin/bash

# cut the reference fasta to only contain chromosomes not contigs
# samtools faidx AnoDarl_H01.genomic.fasta
# samtools faidx AnoDarl_H01.genomic.fasta NC_064873.1 NC_064874.1 NC_064875.1 NC_064612.1 >  AnoDarl_H01.genomic.chromosomesonly.fasta

# if want only mito genome, subset vcf and subset the reference fasta file
# samtools faidx AnoDarl_H01.genomic.fasta NC_064612.1 > AnoDarl_H01_NC_064612.1.fasta


# Input multi-sample VCF file
multi_sample_vcf="renamedchr_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz"

# Output directory for gVCF files
output_dir="gvcf_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Ensure the VCF file is indexed
if [ ! -f "${multi_sample_vcf}.tbi" ]; then
    echo "Indexing VCF file..."
    bcftools index "$multi_sample_vcf"
fi

# Iterate over samples in the multi-sample VCF
bcftools query -l "$multi_sample_vcf" |
while read -r sample; do
    echo "Processing sample: $sample"

    # Extract sample-specific variants
    bcftools view -s "$sample" "$multi_sample_vcf" -Oz -o "$output_dir/${sample}.vcf.gz"
    bcftools index "$output_dir/${sample}.vcf.gz"

    # Create consensus sequence
    bcftools consensus -f AnoDarl_H01.genomic.chromosomesonly.fasta "$output_dir/${sample}.vcf.gz" |
    awk '/^>/ {next} {printf $0}' |   # Remove headers and concatenate all sequences
    awk -v header=">$sample" 'BEGIN {print header} {print}' > "$output_dir/${sample}_consensus.fasta"
done

echo "Processing completed."

# Concatenate the fasta files into one file
cat "$output_dir"/*_consensus.fasta > allsamples_consensus_wholegenome.fasta

echo "Concatenated fasta files for each sample to produce allsamples_consensus_wholegenome.fasta"
