
cat fastq2vcfsamples.txt | parallel -j 20 \
"/mnt/storage11/sophie/fastq2matrix/scripts/fastq2vcf.py all --read1 {}_1.fastq.gz --read2 {}_2.fastq.gz \
--ref VectorBase-62_AmelasCM1001059_A_Genome.fasta \
--prefix {}" > fastq2vcf_log.txt 2>&1
