# FILTERING OF COMBINED GENOTYPED VCF
##### FILTERING OF HOLLY WGS PAPER SAMPLES
# same filters used for both gambiae pipeline and melas pipeline

########### FILTER 1: remove INDELS with bcftools, here -M2 indicates that bcftools should split multi-allelic sites into multiple biallelic sites, 
# keeping this information
# -m2 is used in conjunction with -M2 to apply the minor-allele-based decomposition. 
# This means that bcftools will decompose multi-allelic sites using the minor allele as the reference allele in the biallelic split.
# here -v tells bcftools to only view SNPS, so indels are excluded

bcftools view -M2 -m2 -v snps darlingi.genotyped.vcf.gz -Oz -o bi_snps_darlingi.genotyped.vcf.gz

# you must bgzip the file before indexing if you did not use -Oz to make the output bgzip
# bgzip bi_snps_gambiae_nov2022.2023_07_05.genotyped.vcf

# tabix index the compressed VCF file, creates .vcf.gz.tbi
tabix -f -p vcf bi_snps_darlingi.genotyped.vcf.gz

# rename chromosomes
# create chr_map.txt with two columns, original chr name in first column and new chr name in second column
bcftools annotate --rename-chrs chr_map.txt bi_snps_darlingi.genotyped.vcf.gz -Oz -o renamedchr_bi_snps_darlingi.genotyped.vcf.gz

tabix -p vcf renamedchr_bi_snps_darlingi.genotyped.vcf.gz

# filter out the contigs from the VCF file, note that to produce a properly bgzipped vcf file you need the -Oz flag
bcftools view renamedchr_bi_snps_darlingi.genotyped.vcf.gz --regions anop_X,2,3,anop_mito | bcftools sort -Oz -o filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

# view the unique chromosome names
bcftools query -f '%CHROM\n' filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | sort | uniq > unique_chromosomes_filtered.txt

tabix -p vcf filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

########### FILTER 2: Filter samples to keep those with 40% of genome with > 10x coverage, and min-ac=1 so that all variants that remain are still variants after sample removal
# Holly WGS paper samples, keeping those with average coverage over 5.

# create file with samples to keep: sample_to_keep.txt

bcftools view -S samples_to_keep.txt filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz --min-ac=1 -Oz -o minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

tabix -p vcf minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## holly_wgs_paper_darlingi_variants = 14186231
# bcftools view -H minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | wc -l


########### FILTER 3: for GATK standard filtering
# GATK filter recommendations: https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants
# excellent explanation of this hard filtering and additional soft filtering
# https://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/#ex3.4

gatk VariantFiltration \
-R AnoDarl_H01.genomic.fasta \
-V minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## note there is a repeated warning here JEXL engine undefined variable, not a problem
## because this is coming from where there are positions with no coverage
## https://gatk.broadinstitute.org/hc/en-us/community/posts/4408733963803-GATK-Variant-Filtration-undefined-variable

## remove the SNPs that have been tagged with the above filters
bcftools view -f 'PASS' gatk_tagged_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -Oz -o gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## holly_wgs_paper_darlingi_variants = 13730595
# bcftools view -H gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | wc -l

## ADDITIONAL SOFT FILTERING

## After hard quality filtering, we have a data set containing only variant sites that we trust with high confidence. 
## However, so far we have only removed variants that had skewed values across all samples – 
## we haven’t yet looked at the individual genotypes at each variant site. 
## For example, we may have kept a variant that has a high quality across all of the individual samples, but there 
## may be some individuals where the read depth (FMT/DP) or the genotype quality (GQ) for the individual genotype is very low, 
## so we don’t trust the calls for these individuals

########### FILTER 4: Filter reads that have a read depth below 5 OR a genotype quality below 20
## The difference is, that instead of the whole variant site, we are now considering single genotypes 
## (the information in the “FORMAT” fields”) using -e 'FMT/DP<3 | FMT/GQ<20' and we will not remove the whole variant site, 
## but only set the respective genotypes to missing (./.) by using the bcftools filter -S . command.
## Here we use the logical operator “|” instead of “||”, since using “||” would mean that every genotype at a variant site 
## is set to missing, even if only one genotype doesn’t fulfill the threshold
bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## check this, the genotype here should have been set to missing ./. as we are filtering for depth below 5
#bcftools query -i 'FMT/DP<5' -f '[GT=%GT\tDP=%DP\tGQ=%GQ\t]\n' DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | less -S

## happy? index vcf
tabix -p vcf DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## holly_wgs_paper_darlingi_variants: 13730595. Same number as all our filters above just set to missing, so genotypes removed but site still in VCF.
# bcftools view -H DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | wc -l

########### FILTER 5: exclude -e all sites at which no alternative alleles are called for any of the samples
bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

tabix -p vcf AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

## holly_wgs_paper_darlingi_variants: 12510546
# bcftools view -H AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | wc -l

########### FILTER 6 ##########

## Remove variants with a high amount of missing genotypes and filter on minor allele frequency
## The data set overall will now have lots of missing data, because we have replaced calls with depth <5 or quality below 20 with ./.  
## Therefore we will now remove all variants that have more than 20% missing genotypes or MAF < 0.01
## Filtering for MAF < 0.01 means that remaining sites in the VCF need to have a minimum minor allele frequency of 1%

bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.01' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

# check this has worked. The minor allele count should be 1 or above in this dataset.
# bcftools query F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -f'%AC\n' | sort -g | head

bcftools query F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz -f'%AC\n' | sort -g | head

# final SNPs: 


# Filtering complete!

########## PHASE VCF FILE ###########

# Phase the filtered vcf using beagle, https://faculty.washington.edu/browning/beagle/beagle_5.2_13Oct21.pdf
# necessary for some selection statistical analyses
# could also use phasing pipeline from malariagen. Have a look at this, and understand gatk parameters above.

# use beagle conda environment. mamba install beagle.
beagle -Xmx500g gt=F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz out=darlingi_filtered_phased

tabix -p vcf darlingi_filtered_phased.vcf.gz

## Check the number of SNPs in the phased file, it should be the same.

# unphased darlingi, chromosome 3: 31996
bcftools query -f '%CHROM\t%POS\n' F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz | awk '$1=="3"' | wc -l
# phased darlingi, chromosome 3: 31996
bcftools query -f '%CHROM\t%POS\n' darlingi_filtered_phased.vcf.gz | awk '$1=="3"' | wc -l


## SnpEff annotation of filtered VCF (this one was unphased, could do either)
snpEff Anopheles_darlingi_2 renamedchr_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz > renamedchr_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf
bgzip F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf
tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.ann.vcf.gz


###########################################################################################################
############################### Phenotyped Only Samples VCF and Filtering #################################
###########################################################################################################

## Subset vcfs to create vcf containing only the phenotyped samples

bcftools view -S phenotyped_samples.txt F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz --min-ac=1 -Oz -o minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz
tabix -p vcf minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz

gatk VariantFiltration \
-R AnoDarl_H01.genomic.fasta \
-V minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz \
-filter "QD < 5.0" --filter-name "QD5" \
-filter "QUAL < 30.0" --filter-name "QUAL30" \
-filter "SOR > 3.0" --filter-name "SOR3" \
-filter "FS > 60.0" --filter-name "FS60" \
-filter "MQ < 40.0" --filter-name "MQ40" \
-filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
-filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
-O gatk_tagged_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz

## remove the SNPs that have been tagged with the above filters
bcftools view -f 'PASS' gatk_tagged_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz -Oz -o gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz
# count variants: 3999253

bcftools filter -S . -e 'FMT/DP<5 | FMT/GQ<20' -O z -o DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz
tabix -p vcf DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz

########### FILTER 5: exclude -e all sites at which no alternative alleles are called for any of the samples
bcftools filter -e 'AC==0' DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz -O z -o AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz
tabix -p vcf AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz

########### FILTER 6 ##########
bcftools filter -e 'F_MISSING > 0.2 || MAF <= 0.02' -O z -o F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz
tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz

# Filtering complete!
# Phase VCF file
beagle -Xmx500g gt=F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz out=phenotyped_darlingi_filtered_phased
tabix -p vcf phenotyped_darlingi_filtered_phased.vcf.gz

## Check the number of SNPs in the phased file, it should be the same.

# unphased darlingi, chromosome 3: 1656594
bcftools query -f '%CHROM\t%POS\n' F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz | awk '$1=="3"' | wc -l
# phased darlingi, chromosome 3: 1656594
bcftools query -f '%CHROM\t%POS\n' phenotyped_darlingi_filtered_phased.vcf.gz | awk '$1=="3"' | wc -l

## SnpEff annotation of filtered VCF (this one was unphased, could do either)
snpEff Anopheles_darlingi F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.vcf.gz > F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.ann.vcf
bgzip F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.ann.vcf
tabix -p vcf F_MISSING_MAF_AC0_DP5_GQ20_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi_phenotyped.genotyped.ann.vcf.gz

## Check the number of SNPs in the Rondonia Samples and the number of SNPs in the Colony Samples
# Colony samples vcf

colony_old_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

bcftools filter -e 'AC==0' -Oz -o colony_samples_only_SNPS.vcf.gz colony_old_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

# Rondonia samples vcf

rondonia_samples_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz

bcftools filter -e 'AC==0' -Oz -o rondonia_samples_only_SNPS.vcf.gz rondonia_samples_F_MISSING_MAF_AC0_DP5_GQ15_gatk_filtered_minac_filtrenamedchr_bi_snps_darlingi.genotyped.vcf.gz
