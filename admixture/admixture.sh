## Admixture

# Identify which samples are in your vcf:
bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING
# MAF > 0.01 filter has already been applied

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | cut -f1 | sort | uniq

zcat F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^anop_mito$/, "5", $1); gsub(/^anop_X$/, "6", $1); print}' | bgzip > admixture_modified.vcf.gz

tabix -p vcf admixture_modified.vcf.gz
zgrep -v "^#" admixture_modified.vcf.gz | cut -f1 | sort | uniq
bcftools query -l admixture_modified.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture_modified.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out melas_global_gambiaealigned


# STEP 4 RUN ADMIXTURE

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 14062 melas_global_gambiaealigned.bed {} | tee log{}.cv10.seed14062.out'

# Inspected inflection point with grep -h CV *out, and K = 3

# Plot admixture
conda create -n radmix r-essentials r-base
install.packages(c("unikn", "countrycode", "optparse"))
setwd("/mnt/storage11/sophie/darlingi/darlingi_database")

Rscript /mnt/storage11/sophie/gitrepos/anopheles_darlingi/admixture/generate_admix_barplot_colours_v2.R \
-d /mnt/storage11/sophie/darlingi/darlingi_database \
--prefix melas_global_gambiaealigned \
--kval 3 \
-m /mnt/storage11/sophie/darlingi/darlingi_metadata_admixture_v2.csv \
--filter_N 1 \
--label_id sample \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE