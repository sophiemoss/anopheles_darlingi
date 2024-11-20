#!/bin/bash

# Input multi-sample VCF file
multi_sample_vcf="mito_only_filtered_darlingi_genotyped_NC_064612.1.vcf.gz"

# Output directory for gVCF files
output_dir="gvcf_output"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Extract SNPs from the multi-sample VCF
bcftools view -V indels -c 1 -a "$multi_sample_vcf" -Oz -o "${multi_sample_vcf%.vcf.gz}.snps.vcf.gz"
bcftools index "${multi_sample_vcf%.vcf.gz}.snps.vcf.gz"

# Iterate over samples in the multi-sample VCF
bcftools query -l "$multi_sample_vcf" |
while read -r sample; do
    echo "Processing sample: $sample"

    # Extract sample-specific variants
    bcftools view -s "$sample" "${multi_sample_vcf%.vcf.gz}.snps.vcf.gz" -Oz -o "$output_dir/${sample}.vcf.gz"
    bcftools index "$output_dir/${sample}.vcf.gz"

    # Create consensus sequence
    bcftools consensus -f AnoDarl_H01_NC_064612.1.fasta -m "$output_dir/${sample}_mask.bed" "$output_dir/${sample}.vcf.gz" |
    tr '*' 'N' |
    sed 's/>.*/>'"$sample"'/' > "$output_dir/${sample}_consensus.fasta"

done

echo "Processing completed."
