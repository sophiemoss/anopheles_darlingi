########################################
######## R CODE: PCA of admixture
########################################
library("SNPRelate")

# Aedes aegypti
snpgdsBED2GDS("gyp_imputed.bed", "gyp_imputed.fam", "gyp_imputed.bim", "gyp_imputed.gds", option = snpgdsOption(autosome.end = 60000))
#note the .bim file must be edited to conform with the standard inputs for chromosome names, etc

# open dataset
genofile <- snpgdsOpen("gyp_imputed.gds")
sample_newcal=read.delim("newcaledonia.txt",header = F)
sample_bravan=read.delim("brazil_vanuatu.txt",header = F)

#all snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=8, sample.id = sample_bravan$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_nc <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_newcal$V1)
SL_bravan <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_bravan$V1)

write.csv(SL_nc$eigenvect,"pca_newcal.txt")
write.csv(SL_bravan$eigenvect,"pca_bravan.txt")


# Aedes albopictus
snpgdsBED2GDS("alb_imputed.bed", "alb_imputed.fam", "alb_imputed.bim", "alb_imputed.gds", option = snpgdsOption(autosome.end = 60000))
#note the .bim file must be edited to conform with the standard inputs for chromosome names, etc

# open dataset
genofile <- snpgdsOpen("alb_imputed.gds")
sample_norfly=read.delim("northfly.txt",header = F)
sample_laetsi=read.delim("lae_torresstrait.txt",header = F)

#all snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=8, sample.id = sample_laetsi$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_nc <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_norfly$V1)
SL_bravan <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_laetsi$V1)

write.csv(SL_nc$eigenvect,"pca_norfly.txt")
write.csv(SL_bravan$eigenvect,"pca_laetsi.txt")



########################################
######## R CODE: latent factor mixed models
########################################
install.packages("LEA")
library("LEA")

# Anopheles gambiae
# First, convert VCF to .geno format.
# vcf2geno converts a vcf to a .geno format which is compatible with the LEA package.
# the force = TRUE parameter overwrites any existing file if necessary.
vcf2geno(input.file = "2022gambiaevcfphased.vcf", output.file = "admix.geno", force = TRUE)
vcf2geno(input.file = "emma_PR.vcf", output.file = "admix.geno", force = TRUE)
# admix_pop=read.delim("vssc395ids.txt",row.names = 1,header = TRUE) # this line appears to be redundant

# Run sNMF to find appropriate K.
# sNMF is run to estimate the ancestry components (K). Here, it tried K values from 1:25, with 10 repetitions each
# and uses 25 CPUs for parallel processing. The entropy = TRUE option calculates cross-entropy, which is then used later 
# to calculate the best k. The project is saved as xxx.snmfProject.

project = NULL
project = snmf(input.file = "mini.admix.geno", K = 1:25, entropy = TRUE,
               repetitions = 10, CPU = 25, project = "new")

project = load.snmfProject("mini.admix.snmfProject")

# Plot cross-entropy to select the best K. This helps to visualise the cross-entropy score, indicating the 
# best K for population structure.

png("mini.admixture_plot.png", width = 800, height = 600)
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
# best = which.min(cross.entropy(project, K = 3)) # this line appears to be redundant, not used again later

## Now run the lfmm using this selected K value for environmental association testing. 
# The *.env file is a vector of values 0-2 indicating number of copies of the resistance allele
# Read in the .env file
L995F.mini.env <- read.table("2422652_L995F_mini.env", quote="\"", comment.char="")
# Use vcf2lfmm to convert the VCF file to the .lfmm format needed for LFMM.
vgsc.mini.lfmm = vcf2lfmm("minivcf_2022gambiaephased.vcf")

L995F.mini.lfmm <- lfmm("minivcf_2022gambiaephased.lfmm", "2422652_L995F_mini.env", K = 3, CPU = 16, rep = 10, project="force", 
                  iterations = 100, burnin = 5, all = FALSE, missing.data = TRUE)

project = load.lfmmProject("minivcf_2022gambiaephased_2422652_L995F_mini.lfmmProject")

# Extract and analyse Z-scores for association testing
# Collect z-scores from the 10 replicates into a matrix (zs_GTC)
zs_GTC = z.scores(L995F.mini.lfmm, K = 3, d = 1)

# Combine z-scores using the median
# The zs.median vector stores the median Z-score per locus
zs.median = apply(zs_GTC, MARGIN = 1, median)

# Lamda is calculated as a genomic inflation factor (GIF) to correct for populations structure
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

# Compute adjusted p-values from the combined z-scores, and save a histogram of these adjusted p-values.
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
png("histogram_plot_fewerits.png")
hist(adj.p.values, col = "red")
dev.off()


## FDR control: Implements the Benjamini-Hochberg FDR method at level q = 0.01, to identify
## significant loci after multiple testing. It sorts the adjusted p-values and finds loci where significance
## meets the FDR threshold. Finally, the adjusted p-values are saved to a CSV file.
## L = number of loci
L = 16452859
#fdr level q
q = 0.01
w = which(sort(adj.p.values) < q * (1:L)/L)
candidates.bh = order(adj.p.values)[w]
candidates.bh
saveRDS(candidates.bh, file = "candidates_bh.rds")
# load candidates.bh with:
# candidates.bh <- readRDS("candidates_bh.rds")

write.csv(adj.p.values, "adj_p_values_k3_vcf_fewer_iterations")
# this is a csv which contains the adjusted p-values for each SNP in the VCF file
# these p-values are based on associations between loci and the environmental variable specified using K populations
# you can then use the candidates.bh array to identify the indices of SNPs in the VCF file which are significant

## Use candidates.bh and the adj_p_values csv to identify which SNPs in the VCF are significantly associated with the env variable.
library(VariantAnnotation)
vcf <- readVcf("2022gambiaevcfphased.vcf")

# Extract the positions from candidates.bh
# Get the chromosome and position information for significant loci from candidates.bh
significant_positions <- data.frame(CHROM = seqnames(vcf)[candidates.bh],POS = start(vcf)[candidates.bh])
# Save these positions to a text file (required by bcftools)
# The resulting significant_positions.txt file will contain two columns (CHROM and POS) for each significant locus
write.table(significant_positions, "significant_positions.txt", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)
# Subet the VCF using bcftools to create a new VCF with only significant loci
system("bcftools view -R significant_positions.txt 2022gambiaevcfphased.vcf.gz -o significant_L995F_2022gambiaevcfphased_fewerits.vcf")

# In this situation there were 244408 candidate positions out of 16452859 original positions (0.01485505 = 1.49% of positions)

# Plot the significant positions, make one with the p-values in and then work out how to plot it 
significant_positions

########################################
######## R CODE: PCA of VSSC outliers
########################################
library("SNPRelate")

snpgdsBED2GDS("gypVSSC.bed", "gypVSSC.fam", "gypVSSC.bim", "gypVSSC.gds", option = snpgdsOption(autosome.end = 60000))

# open dataset
genofile <- snpgdsOpen("gypVSSC.gds")
samples=read.delim("gvpVSSC_homozygotes.txt",header = F)

##all vssc snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=2, sample.id = sample_all$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_hom <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samples$V1)

write.csv(cbind(SL_hom$sample.id, SL_hom$eigenvect),"pca_gypVSSC.txt")


########################################
######## R CODE: Tajima's D, pi, LD, Fst
########################################
library(tidyverse)
library(tidyquant)

# V1016G -  an example. 

#####  Tajima's D
V1016G_D = read.delim("V1016G.Tajima.D",header = TRUE)
V1016G_D_wt = read.delim("totWT.Tajima.D",header = TRUE)

V1016G_D <- V1016G_D  %>% mutate(BIN_STARTmb = (V1016G_D$BIN_START/1000000))
V1016G_D = V1016G_D  %>% mutate(rollD = (zoo::rollapply(V1016G_D$TajimaD, 55, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_D_wt = V1016G_D_wt %>% mutate(rollD = (zoo::rollapply(V1016G_D_wt$TajimaD, 55, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_D <- V1016G_D  %>% mutate(deltaD = (V1016G_D$rollD-V1016G_D_wt$rollD))
V1016G_D <- V1016G_D  %>% mutate(deltaDraw = (V1016G_D$TajimaD-V1016G_D_wt$TajimaD))


#####   PI
# V1016G 
V1016G_pi = read.delim("V1016G.sites.pi",header = TRUE)
V1016G_pi_wt = read.delim("totWT.sites.pi",header = TRUE)

V1016G_pi <- V1016G_pi  %>% mutate(POSmb = (V1016G_pi$POS/1000000))
V1016G_pi <- V1016G_pi  %>% mutate(deltapi = (V1016G_pi$PI-V1016G_pi_wt$PI))
V1016G_pi = V1016G_pi  %>% mutate(rolldeltapi = (zoo::rollapply(V1016G_pi$deltapi, 10, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_pi$rolldeltapi<-ifelse(V1016G_pi$rolldeltapi==0,NA,V1016G_pi$rolldeltapi)


#####   LD
V1016G_LD = read.delim("V1016G.geno.ld",header = TRUE)
V1016G_LD_wt = read.delim("totWT.geno.ld",header = TRUE)

V1016G_LD <- V1016G_LD  %>% mutate(POSmb = ((V1016G_LD$POS1+V1016G_LD$POS2)/2000000))
V1016G_LD = V1016G_LD  %>% mutate(rollD = (zoo::rollapply(V1016G_LD$R.2, 500, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_LD_wt = V1016G_LD_wt %>% mutate(rollD = (zoo::rollapply(V1016G_LD_wt$R.2, 500, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_LD <- V1016G_LD  %>% mutate(deltaLD = (V1016G_LD$R.2-V1016G_LD_wt$R.2))
V1016G_LD <- V1016G_LD  %>% mutate(rolldeltaLD = (V1016G_LD$rollD-V1016G_LD_wt$rollD))
########################################
######## R CODE: PCA of admixture
########################################
library("SNPRelate")

# Aedes aegypti
snpgdsBED2GDS("gyp_imputed.bed", "gyp_imputed.fam", "gyp_imputed.bim", "gyp_imputed.gds", option = snpgdsOption(autosome.end = 60000))
#note the .bim file must be edited to conform with the standard inputs for chromosome names, etc

# open dataset
genofile <- snpgdsOpen("gyp_imputed.gds")
sample_newcal=read.delim("newcaledonia.txt",header = F)
sample_bravan=read.delim("brazil_vanuatu.txt",header = F)

#all snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=8, sample.id = sample_bravan$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_nc <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_newcal$V1)
SL_bravan <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_bravan$V1)

write.csv(SL_nc$eigenvect,"pca_newcal.txt")
write.csv(SL_bravan$eigenvect,"pca_bravan.txt")


# Aedes albopictus
snpgdsBED2GDS("alb_imputed.bed", "alb_imputed.fam", "alb_imputed.bim", "alb_imputed.gds", option = snpgdsOption(autosome.end = 60000))
#note the .bim file must be edited to conform with the standard inputs for chromosome names, etc

# open dataset
genofile <- snpgdsOpen("alb_imputed.gds")
sample_norfly=read.delim("northfly.txt",header = F)
sample_laetsi=read.delim("lae_torresstrait.txt",header = F)

#all snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=8, sample.id = sample_laetsi$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_nc <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_norfly$V1)
SL_bravan <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=sample_laetsi$V1)

write.csv(SL_nc$eigenvect,"pca_norfly.txt")
write.csv(SL_bravan$eigenvect,"pca_laetsi.txt")



########################################
######## R CODE: latent factor mixed models
########################################
install.packages("LEA")
library("LEA")

# Anopheles gambiae
# First, convert VCF to .geno format.
# vcf2geno converts a vcf to a .geno format which is compatible with the LEA package.
# the force = TRUE parameter overwrites any existing file if necessary.
vcf2geno(input.file = "2022gambiaevcfphased.vcf", output.file = "admix.geno", force = TRUE)
# admix_pop=read.delim("vssc395ids.txt",row.names = 1,header = TRUE) # this line appears to be redundant

# Run sNMF to find appropriate K.
# sNMF is run to estimate the ancestry components (K). Here, it tried K values from 1:25, with 10 repetitions each
# and uses 25 CPUs for parallel processing. The entropy = TRUE option calculates cross-entropy, which is then used later 
# to calculate the best k. The project is saved as xxx.snmfProject.

project = NULL
project = snmf(input.file = "admix.geno", K = 1:25, entropy = TRUE,
               repetitions = 10, CPU = 25, project = "new")

project = load.snmfProject("admix.snmfProject")

# Plot cross-entropy to select the best K. This helps to visualise the cross-entropy score, indicating the 
# best K for population structure.

png("admixture_plot.png", width = 800, height = 600)
plot(project, col = "blue", pch = 19, cex = 1.2)
dev.off()
# best = which.min(cross.entropy(project, K = 3)) # this line appears to be redundant, not used again later

## Now run the lfmm using this selected K value for environmental association testing. 
# The *.env file is a vector of values 0-2 indicating number of copies of the resistance allele
# Read in the .env file
L995F.env <- read.table("2422652_L995F_newer.env", quote="\"", comment.char="")
# Use vcf2lfmm to convert the VCF file to the .lfmm format needed for LFMM.
vgsc.lfmm = vcf2lfmm("2022gambiaephased.vcf")

L995F.lfmm <- lfmm("2022gambiaephased.lfmm", "2422652_L995F_newer.env", K = 3, CPU = 16, rep = 10, project="force", 
                  iterations = 100, burnin = 5, all = FALSE, missing.data = TRUE)

project = load.lfmmProject("2022gambiaevcfphased_2422652_L995F_newer.lfmmProject")

# Extract and analyse Z-scores for association testing
# Collect z-scores from the 10 replicates into a matrix (zs_GTC)
zs_GTC = z.scores( , K = 5, d = 1)

# Combine z-scores using the median
# The zs.median vector stores the median Z-score per locus
zs.median = apply(zs_GTC, MARGIN = 1, median)

# Lamda is calculated as a genomic inflation factor (GIF) to correct for populations structure
lambda = median(zs.median^2)/qchisq(0.5, df = 1)
lambda

# Compute adjusted p-values from the combined z-scores, and save a histogram of these adjusted p-values.
adj.p.values = pchisq(zs.median^2/lambda, df = 1, lower = FALSE)
png("histogram_plot.png")
hist(adj.p.values, col = "red")
dev.off()

## FDR control: Implements the Benjamini-Hochberg FDR method at level q = 0.01, to identify
## significant loci after multiple testing. It sorts the adjusted p-values and finds loci where significance
## meets the FDR threshold. Finally, the adjusted p-values are saved to a CSV file.
## L = number of loci
L = 16452859
#fdr level q
q = 0.01
w = which(sort(adj.p.values) < q * (1:L)/L)
big_candidates.bh = order(adj.p.values)[w]
big_candidates.bh
saveRDS(big_candidates.bh, file = "big_candidates_bh.rds")
# load candidates.bh with:
big_candidates.bh <- readRDS("big_candidates_bh.rds")

write.csv(adj.p.values, "adj_p_values_k3_vcf_big")
# this is a csv which contains the adjusted p-values for each SNP in the VCF file
# these p-values are based on associations between loci and the environmental variable specified using K populations
# you can then use the candidates.bh array to identify the indices of SNPs in the VCF file which are significant

## Use candidates.bh and the adj_p_values csv to identify which SNPs in the VCF are significantly associated with the env variable.
library(VariantAnnotation)
vcf <- readVcf("2022gambiaevcfphased.vcf")

# Extract the positions from candidates.bh
# Get the chromosome and position information for significant loci from candidates.bh
significant_positions_bigvcf <- data.frame(CHROM = seqnames(vcf)[big_candidates.bh],POS = start(vcf)[big_candidates.bh])
# Save these positions to a text file (required by bcftools)
# The resulting significant_positions.txt file will contain two columns (CHROM and POS) for each significant locus
write.table(significant_positions_bigvcf, "significant_positions_bigvcf.txt", sep = "\t", 
            col.names = FALSE, row.names = FALSE, quote = FALSE)

# Subet the VCF using bcftools to create a new VCF with only significant loci
system("bcftools view -R significant_positions_bigvcf.txt 2022gambiaevcfphased.vcf.gz -o significant_L995F_2022gambiaevcfphased_bigvcf.vcf")

# In the fewer iterations vcf there were 244408 candidate positions out of 16452859 original positions (0.01485505 = 1.49% of positions)
# In the big vcf there were 226243 candidate positions out of 16452859 original positions (1.38 % of positions) 
significant_positions_bigvcf

# Add p-values to these significant positions
significant_pvalues_bigvcf <-adj.p.values[big_candidates.bh]
# Add them as a new column
significant_positions_bigvcf$p_value <- significant_pvalues_bigvcf
# View updated dataframe
head(significant_positions_bigvcf)
# Save the updated dataframe to a file if needed
write.table(significant_positions_bigvcf, "significant_positions_with_pvalues_bigvcf.txt", sep = "\t", 
            col.names = TRUE, row.names = FALSE, quote = FALSE)

# Now plot the significant p-values over their positions for each chromosome

# Load necessary libraries
library(ggplot2)
library(dplyr)

# Load the data from significant_positions_with_pvalues.txt
significant_positions_bigvcf <- read.table("significant_positions_with_pvalues_bigvcf.txt", header = TRUE, sep = "\t")

# Define chromosome lengths based on GFF information
chromosome_lengths <- data.frame(
  CHROM = c("2L", "2R", "3L", "3R", "Mt", "X"),
  LENGTH = c(49364325, 61545105, 41963435, 53200684, 15363, 24393108)
)

# Join chromosome lengths with significant positions data
significant_positions_bigvcf <- significant_positions_bigvcf %>%
  inner_join(chromosome_lengths, by = "CHROM")

# Create a plot for each chromosome
# We use -log10(p_value) to make small p-values (significant) more prominent on the plot
# Create a list to store plots and save each plot as an individual file
plots <- list()
for (chrom in unique(significant_positions_bigvcf$CHROM)) {
  # Filter data for the current chromosome
  chromosome_data <- significant_positions_bigvcf %>% filter(CHROM == chrom)
  
  # Create the plot for the current chromosome
  p <- ggplot(chromosome_data, aes(x = POS, y = -log10(p_value))) +
    geom_point(color = "blue", alpha = 0.5) +
    scale_x_continuous(limits = c(0, chromosome_data$LENGTH[1])) +
    labs(title = paste("Chromosome", chrom), x = "Position", y = "-log10(p-value)") +
    theme_minimal()
  
  # Store the plot in the list
  plots[[chrom]] <- p
  
  # Save the plot to a file (e.g., PNG format)
  ggsave(filename = paste0("Chromosome_", chrom, "_pvalue_plot_big_vcf.png"), plot = p, width = 10, height = 6)
}













########################################
######## R CODE: PCA of VSSC outliers
########################################
library("SNPRelate")

snpgdsBED2GDS("gypVSSC.bed", "gypVSSC.fam", "gypVSSC.bim", "gypVSSC.gds", option = snpgdsOption(autosome.end = 60000))

# open dataset
genofile <- snpgdsOpen("gypVSSC.gds")
samples=read.delim("gvpVSSC_homozygotes.txt",header = F)

##all vssc snps
PCARV <- snpgdsPCA(genofile, eigen.cnt=2, sample.id = sample_all$V1)
SnpLoad <- snpgdsPCASNPLoading(PCARV, genofile)

# calculate sample eigenvectors from SNP loadings
SL_hom <- snpgdsPCASampLoading(SnpLoad, genofile, sample.id=samples$V1)

write.csv(cbind(SL_hom$sample.id, SL_hom$eigenvect),"pca_gypVSSC.txt")


########################################
######## R CODE: Tajima's D, pi, LD, Fst
########################################
library(tidyverse)
library(tidyquant)

# V1016G -  an example. 

#####  Tajima's D
V1016G_D = read.delim("V1016G.Tajima.D",header = TRUE)
V1016G_D_wt = read.delim("totWT.Tajima.D",header = TRUE)

V1016G_D <- V1016G_D  %>% mutate(BIN_STARTmb = (V1016G_D$BIN_START/1000000))
V1016G_D = V1016G_D  %>% mutate(rollD = (zoo::rollapply(V1016G_D$TajimaD, 55, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_D_wt = V1016G_D_wt %>% mutate(rollD = (zoo::rollapply(V1016G_D_wt$TajimaD, 55, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_D <- V1016G_D  %>% mutate(deltaD = (V1016G_D$rollD-V1016G_D_wt$rollD))
V1016G_D <- V1016G_D  %>% mutate(deltaDraw = (V1016G_D$TajimaD-V1016G_D_wt$TajimaD))


#####   PI
# V1016G 
V1016G_pi = read.delim("V1016G.sites.pi",header = TRUE)
V1016G_pi_wt = read.delim("totWT.sites.pi",header = TRUE)

V1016G_pi <- V1016G_pi  %>% mutate(POSmb = (V1016G_pi$POS/1000000))
V1016G_pi <- V1016G_pi  %>% mutate(deltapi = (V1016G_pi$PI-V1016G_pi_wt$PI))
V1016G_pi = V1016G_pi  %>% mutate(rolldeltapi = (zoo::rollapply(V1016G_pi$deltapi, 10, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_pi$rolldeltapi<-ifelse(V1016G_pi$rolldeltapi==0,NA,V1016G_pi$rolldeltapi)


#####   LD
V1016G_LD = read.delim("V1016G.geno.ld",header = TRUE)
V1016G_LD_wt = read.delim("totWT.geno.ld",header = TRUE)

V1016G_LD <- V1016G_LD  %>% mutate(POSmb = ((V1016G_LD$POS1+V1016G_LD$POS2)/2000000))
V1016G_LD = V1016G_LD  %>% mutate(rollD = (zoo::rollapply(V1016G_LD$R.2, 500, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_LD_wt = V1016G_LD_wt %>% mutate(rollD = (zoo::rollapply(V1016G_LD_wt$R.2, 500, mean, partial=T, fill=NA, na.rm = TRUE)))
V1016G_LD <- V1016G_LD  %>% mutate(deltaLD = (V1016G_LD$R.2-V1016G_LD_wt$R.2))
V1016G_LD <- V1016G_LD  %>% mutate(rolldeltaLD = (V1016G_LD$rollD-V1016G_LD_wt$rollD))

