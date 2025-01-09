## Create a distance matrix for PCA
plink --vcf rondonia_only_darlingi_filtered.vcf.gz --distance square --double-id --out wholegenome_rondonia_darlingi --threads 6 --allow-extra-chr

## now create for every chromosome
## using just the mitochondria
## create mitochondria only vcf using bcftools, then use plink to make distance matrix
bcftools view -r anop_mito rondonia_only_darlingi_filtered.vcf.gz -Oz -o mito_only_rondonia_only_darlingi_filtered.vcf.gz
tabix -f -p vcf mito_only_rondonia_only_darlingi_filtered.vcf.gz
plink --vcf mito_only_rondonia_only_darlingi_filtered.vcf.gz --distance square --double-id --out mito_only_rondonia_darlingi --threads 6 --allow-extra-chr

## just chromosome X
bcftools view -r anop_X rondonia_only_darlingi_filtered.vcf.gz -Oz -o X_only_rondonia_only_darlingi_filtered.vcf.gz
tabix -f -p vcf X_only_rondonia_only_darlingi_filtered.vcf.gz
plink --vcf X_only_rondonia_only_darlingi_filtered.vcf.gz --distance square --double-id --out X_only_rondonia_darlingi --threads 6 --allow-extra-chr

## just chromosome 2
bcftools view -r 2 rondonia_only_darlingi_filtered.vcf.gz -Oz -o 2_only_rondonia_only_darlingi_filtered.vcf.gz
tabix -f -p vcf 2_only_rondonia_only_darlingi_filtered.vcf.gz
plink --vcf 2_only_rondonia_only_darlingi_filtered.vcf.gz --distance square --double-id --out 2_only_rondonia_darlingi --threads 6 --allow-extra-chr

## just chromosome 3
bcftools view -r 3 rondonia_only_darlingi_filtered.vcf.gz -Oz -o 3_only_rondonia_only_darlingi_filtered.vcf.gz
tabix -f -p vcf 3_only_rondonia_only_darlingi_filtered.vcf.gz
plink --vcf 3_only_rondonia_only_darlingi_filtered.vcf.gz --distance square --double-id --out 3_only_rondonia_darlingi --threads 6 --allow-extra-chr
