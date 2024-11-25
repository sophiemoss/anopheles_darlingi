## download additional samples from ncbi using fasterq-dump
conda activate fasterq-dump
fasterq-dump SRR000001

## to make a genomics database of sample VCFs, use the following
ls *.g.vcf.gz | sed 's/.g.vcf.gz//' > fastq2vcfsamples.txt

/mnt/storage12/nbillows/Projects/Malaria/fastq2matrix/scripts/merge_vcfs.py import --sample-file fastq2vcfsamples.txt --ref AnoDarl_H01.genomic.fasta --prefix darlingi --vcf-dir .
/mnt/storage11/sophie/fastq2matrix/scripts/merge_vcfs.py import  --sample-file fastq2vcfsamples.txt --ref AnoDarl_H01.genomic.fasta --prefix darlingi --vcf-dir .

## 28/10/24 got to here! error no attribute bed
## now merge VCF files

/mnt/storage12/nbillows/Projects/Malaria/fastq2matrix/scripts/merge_vcfs.py genotype --ref AnoDarl_H01.genomic.fasta --prefix darlingi > mergevcf_log.txt 2>&1

# resulting vcf is called gambiae_bijagos_2022.2023_07_25.genotyped.vcf.gz