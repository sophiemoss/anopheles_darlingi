## Admixture

# Identify which samples are in your vcf:
#bcftools query -l yourfile.vcf.gz

# STEP 1 EXTRA FILTERING
# MAF > 0.01 filter has already been applied

# STEP 2 CONVERT CHR NAMES IN VCF TO INTEGERS

zgrep -v "^#" F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz | cut -f1 | sort | uniq

zcat F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz | awk 'BEGIN{OFS=FS="\t"} /^#/ {print; next} {gsub(/^anop_mito$/, "4", $1); gsub(/^anop_X$/, "5", $1); print}' | bgzip > admixture_modified.vcf.gz

tabix -p vcf admixture_modified.vcf.gz
zgrep -v "^#" admixture_modified.vcf.gz | cut -f1 | sort | uniq
bcftools query -l admixture_modified.vcf.gz

# STEP 3 MAKE BED AND BIM FILES

plink --vcf admixture_modified.vcf.gz --set-missing-var-ids @:# --keep-allele-order --const-fid --allow-extra-chr --make-bed --out darlingi_phenotyped


# STEP 4 RUN ADMIXTURE
# Make a K_runs.txt file with the numbers of the K you wish to test, eg. 1 to 10

cat K_runs.txt | xargs -I {} sh -c 'admixture --cv=10 -j20 -s 82645 darlingi_phenotyped.bed {} | tee log{}.cv10.seed14062.out'

# Inspected inflection point with 
# grep -h CV *out, and K = 1

# CV error (K=10): 0.84149
# CV error (K=1): 0.48912
# CV error (K=2): 0.50858
# CV error (K=3): 0.54359
# CV error (K=4): 0.56919
# CV error (K=5): 0.59600
# CV error (K=6): 0.63277
# CV error (K=7): 0.70947
# CV error (K=8): 0.74807
# CV error (K=9): 0.79845

# Plot admixture
conda create -n radmix r-essentials r-base
install.packages(c("unikn", "countrycode", "optparse"))
setwd("/mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/phenotyped_colony/admixture")

Rscript /mnt/storage11/sophie/gitrepos/anopheles_darlingi/admixture/generate_admix_barplot_colours_v2_hollywgspaper.R \
-d /mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/phenotyped_colony/admixture \
--prefix darlingi_phenotyped \
--kval 1 \
-m /mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/darlingi_resistance_metadata.csv \
--filter_N 1 \
--label_id sample \
--label_region region \
--label_country country \
--label_site site \
--country_code TRUE