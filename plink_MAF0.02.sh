## Create a distance matrix for PCA
plink --vcf F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz --distance square --double-id --out wholgenome_darlingi --threads 6 --allow-extra-chr

## now create for every chromosome
## using just the mitochondria
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
bcftools view -r anop_mito F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz -Oz -o mito_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
tabix -f -p vcf mito_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
plink --vcf mito_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz --distance square --double-id --out mito_only_darlingi --threads 6 --allow-extra-chr

## just chromosome X
bcftools view -r anop_X F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz -Oz -o X_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
tabix -f -p vcf X_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
plink --vcf X_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz --distance square --double-id --out X_only_darlingi --threads 6 --allow-extra-chr

## just chromosome 2
bcftools view -r 2 F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz -Oz -o 2_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
tabix -f -p vcf 2_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
plink --vcf 2_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz --distance square --double-id --out 2_only_darlingi --threads 6 --allow-extra-chr

## just chromosome 3
bcftools view -r 3 F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz -Oz -o 3_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
tabix -f -p vcf 3_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz
plink --vcf 3_only_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz --distance square --double-id --out 3_only_darlingi --threads 6 --allow-extra-chr
