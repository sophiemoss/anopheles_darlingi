ls *_1.fastq.gz | sed 's/_1.fastq.gz//' > fastq2vcfsamples.txt

cat fastq2vcfsamples.txt | parallel -j 20 \
"/mnt/storage11/sophie/fastq2matrix/scripts/fastq2vcf.py all --read1 {}_1.fastq.gz --read2 {}_2.fastq.gz \
--ref AnoDarl_H01.genomic.fasta \
--prefix {}" > my_env_fastq2vcf_log.txt 2>&1

cat fastq2vcfsamples.txt | parallel -j 20 \
"fastq2vcf.py all --read1 {}_1.fastq.gz --read2 {}_2.fastq.gz \
--ref AnoDarl_H01.genomic.fasta \
--prefix {}" > nina_env_fastq2vcf_log.txt 2>&1