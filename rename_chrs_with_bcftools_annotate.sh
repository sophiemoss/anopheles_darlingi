## Rename chromosomes back again
F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## make chr_map.txt
anop_X    NC_064873.1
anop_mito    NC_064612.1
2    NC_064874.1
3    NC_064875.1

bcftools annotate --rename-chrs chr_map.txt F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -Oz -o renamedchr_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz
bcftools query -f '%CHROM\n' F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt
bcftools csq -p a -f AnoDarl_H01.genomic.fasta -g GCF_943734745.1_idAnoDarlMG_H_01_genomic.gff renamedchr_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -Oz -o darlingi_snps.vcf.gz


snpEff build -gff3 -c /mnt/storage11/sophie/miniconda3/envs/snpeff/share/snpeff-5.1-2/snpEff.config -noCheckCds -noCheckProtein -v Anopheles_darlingi_2

