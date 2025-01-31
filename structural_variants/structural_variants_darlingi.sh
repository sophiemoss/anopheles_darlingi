
## Analysis done here:
## /mnt/storage11/sophie/bijagos_mosq_wgs/2019_melas_fq2vcf_gambiae_aligned/structural_variants

## Using delly call ALL

# read the sample identifiers from all_samples.txt, process the .all.bcf files, in parallel using -P 10,
# then convert this to a vcf and then convert it to a text file without the header lines for easy viewing.
# Produces individual text files for each sample

ls *.mkdup.bam | sed 's/.bam//' > delly_all_samples.txt

for f in *.mkdup.bam; do delly call -t ALL -g AnoDarl_H01.genomic.fasta "$f" -o "${f%.*}.all.bcf"; done

delly merge -o structural_variants.bcf *.all.bcf

bcftools view structural_variants.bcf > structural_variants.vcf

bgzip structural_variants.vcf

## Have a look at the structural_variants.vcf.gz
zgrep -v ^"##" structural_variants.vcf.gz | less

## genotype the structural variants 
cat delly_all_samples.txt | parallel -j 15 --bar "delly call -g AnoDarl_H01.genomic.fasta -v structural_variants.bcf -o {}.all.genotyped.bcf  {}.bam"

# %% merge all of the genotyped bcf files together
bcftools merge -m id -O b -o merged.bcf *.all.genotyped.bcf

# list samples in merged.bcf
bcftools query -l merged.bcf

# convert to vcf to see all structural variants in all samples (unfiltered)
bcftools convert -Oz -o merged.vcf.gz merged.bcf

# check sample names
bcftools view -h merged.vcf.gz | grep '^#CHROM' | cut -f10- > merged_sample_names.txt

# rename

mv merged.vcf.gz merged_genotyped_structural_variants.vcf.gz
tabix -p vcf merged_genotyped_structural_variants.vcf.gz

## Filter to remove SVs with >20% missing data
bcftools view merged_genotyped_structural_variants.vcf.gz | bcftools view -i 'F_PASS(GT!="mis")>0.8' -Oz -o merged_genotyped_structural_variants_miss20.vcf.gz
tabix -p vcf merged_genotyped_structural_variants_miss20.vcf.gz

# Query number of variants
bcftools view merged_genotyped_structural_variants_miss20.vcf.gz | grep -v "^#" | wc -l
# 597055 structural variants in total across entire genome after filtered for missingness.

# Filter the VCF based on genes of interest
bcftools view -R AnDar_allIR_structuralvariants1000.bed merged_genotyped_structural_variants_miss20.vcf.gz -Oz -o genes_merged_genotyped_structural_variants_miss20.vcf.gz
tabix -p vcf genes_merged_genotyped_structural_variants_miss20.vcf.gz
bcftools view genes_merged_genotyped_structural_variants_miss20.vcf.gz | grep -v "^#" | wc -l
# 3310 variants in genes of interest

# Annotate with snpeff

snpEff Anopheles_darlingi_2 genes_merged_genotyped_structural_variants_miss20.vcf.gz > genes_merged_genotyped_structural_variants_miss20.ann.vcf
bgzip genes_merged_genotyped_structural_variants_miss20.ann.vcf
tabix -p vcf genes_merged_genotyped_structural_variants_miss20.ann.vcf.gz

# Extract information from vcf with snpeff
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t%QUAL\t%FILTER\t%INFO/SVTYPE\t%INFO/ANN\t[%GT\t]\n' genes_merged_genotyped_structural_variants_miss20.ann.vcf.gz > structural_variants.ann.tsv

# Use closed excel to open the tsv when no other excels are open, so that it opens the wizard, then make sure in step 1 it says delimted, then step 2 only tsv is ticked.