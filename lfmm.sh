# working out how to run the lfmm model

# creating the .env file for a specific position
# the env file is a vector file which contains information on the number of copies of a resistance allele in each sample, 0, 1 or 2.
# first, this command extracts the genotypes at position 2422652 on chromosome 2L from the VCF and saves them into a file called genotypes.txt.
bcftools query -r 2L:2422652-2422652 -f '[%GT\t]\n' 2022gambiaevcfphased.vcf.gz > 2422652_genotypes.txt
# then convert the genotype codes into numeric allele counts
sed -i 's/0|0/0/g; s/0|1/1/g; s/1|0/1/g; s/1|1/2/g' 2422652_genotypes.txt
# Replace tabs (or spaces) with newlines to ensure each value is on its own line
tr '\t' '\n' < 2422652_genotypes.txt > 2422652_L995F.env
# Verify the resulting .env file
head 2422652_L995F.env
# this file  is then a vector where each row corresponds to a sample, and the values represent the number of copies of the resistance allele for that SNP (position 2422652)