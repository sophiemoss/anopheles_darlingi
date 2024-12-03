## Create vcf of samples to check nucleotide diversity

## subset vcf to have one of Rondonia samples and one of Colony_old samples. Should be from phased vcf!

## generate nucleotide diversity metric

# for multiple vcfs
#for i in $(ls *SRR.vcf.gz);do vcftools --gzvcf $i --window-pi 10000 --out ${i}_nuc_div_window_10kb;done

vcftools --gzvcf colony_old_samples_darlingi_filtered_phased.vcf.gz --window-pi 10000 --out colony_old_nuc_div_window_10kb

# calculate nucleotide diversity per site
#vcftools --gzvcf 3L_bijagos_only_melas_phased.vcf.gz --site-pi --out bijagos_3L_only_nuc_div

# calculate mean from the .pi file that was created

awk 'NR>1 { total += $5; count++ } END { print total/count }' colony_old_nuc_div_window_10kb.windowed.pi

# 0.00091761 (0.09% of nucleotides are different between individuals in this population on average - this makes sense as they are colony samples)

# Nucleotide diversity for Rondonia Samples
vcftools --gzvcf rondonia_samples_darlingi_filtered_phased.vcf.gz --window-pi 10000 --out rondonia_nuc_div_window_10kb

# calculate mean from the .pi file that was created
awk 'NR>1 { total += $5; count++ } END { print total/count }' rondonia_nuc_div_window_10kb.windowed.pi
# 0.00155733 (0.1% of nucleotides are different between individuals in this population on average - these are wild caught so makes sense to be greater % than colony)


## Tajima's D

vcftools --gzvcf colony_old_samples_darlingi_filtered_phased.vcf.gz --TajimaD 20000 --out colony_old_samples_darlingi_filtered_phased

vcftools --gzvcf rondonia_samples_darlingi_filtered_phased.vcf.gz --TajimaD 20000 --out rondonia_samples_darlingi_filtered_phased

# plot using plot_tajima_d.py

# calculate average Tajima D for each population, used pandas to calcualte the average of the TajimaD column. average_td = df["TajimaD"].mean()
# average Tajima D for colony_old samples is 1.0151272463863017 
# average Tajima D for Rondonia samples is -0.8265081070855008 



