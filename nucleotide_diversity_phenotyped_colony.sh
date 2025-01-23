## Create vcf of samples to check nucleotide diversity
## subset vcf to have one of susceptible samples and one of resistant samples. Should be from phased vcf!
## generate nucleotide diversity metric

## original phased vcf is phenotyped_darlingi_filtered_phased.vcf.gz

bcftools view -S resistant_sample_ids.txt phenotyped_darlingi_filtered_phased.vcf.gz -o resistant_phenotyped_darlingi_filtered_phased.vcf.gz

# for multiple vcfs
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

vcftools --gzvcf resistant_phenotyped_darlingi_filtered_phased.vcf.gz --window-pi 10000 --out resistant_phenotyped_darlingi_nuc_div_window_10kb

# calculate nucleotide diversity per site
#vcftools --gzvcf 3L_bijagos_only_melas_phased.vcf.gz --site-pi --out bijagos_3L_only_nuc_div

# calculate mean from the .pi file that was created

awk 'NR>1 { total += $5; count++ } END { print total/count }' resistant_phenotyped_darlingi_nuc_div_window_10kb.windowed.pi

# 0.00091761 (0.09% of nucleotides are different between individuals in this population on average - this makes sense as they are colony samples)

## Tajima's D

vcftools --gzvcf resistant_phenotyped_darlingi_filtered_phased.vcf.gz --TajimaD 20000 --out resistant_phenotyped_darlingi_filtered_phased

vcftools --gzvcf susceptible_phenotyped_darlingi_filtered_phased.vcf.gz --TajimaD 20000 --out susceptible_phenotyped_darlingi_filtered_phased

# plot using plot_tajima_d.py

# calculate average Tajima D for each population, used pandas to calculate the average of the TajimaD column. 
df = pd.read_csv("resistant_phenotyped_darlingi_filtered_phased.Tajima.D", sep = ('\t'))
average_td = df["TajimaD"].mean()

# average Tajima D for resistant samples is 0.5736855705987781 
# average Tajima D for susceptible samples is 0.9723671261031623 