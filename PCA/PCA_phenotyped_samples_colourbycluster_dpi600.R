# Run this as an R script

library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)

workdir <- "/mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/phenotyped_colony/pca" # Working directory with plink files
prefix <- "3_only_darlingi_phenotyped" # Prefix for plink files
metadata <- "/mnt/storage11/sophie/darlingi/phenotype_darlingi_paper/darlingi_resistance_metadata.csv" # File path to metadata

calc_variance_explained <- function(pc_points) {
    vars <- round(pc_points$eig / sum(pc_points$eig) * 100, 1)
    names(vars) <- paste0("PC", seq_len(length(vars)))
    vars
}

# METADATA
met <- read.table(metadata, sep = ",", stringsAsFactors = FALSE, header = TRUE)

#### DIST#
dist <- read.table(file.path(workdir, paste0(prefix, ".dist")), header = FALSE)
id <- read.table(file.path(workdir, paste0(prefix, ".dist.id")))

desc <- id %>% left_join(met, by = c("V1" = "sample"))

dist_m <- as.matrix(dist)
colnames(dist_m) <- desc$V1
rownames(dist_m) <- desc$V1

# PCA #
cmd <- cmdscale(dist_m, k = 10, eig = TRUE, x.ret = TRUE) # Multidimensional Scaling - might take a while
# saveRDS(cmd, paste0(prefix, ".dist.rds") # save to RDS format
#cmd <- readRDS(file.path(workdir, paste0(prefix, ".dist.rds"))
vars <- calc_variance_explained(cmd) # Calculations of variance explained

# Overlay region, country info
# Convert cmdscale output to a dataframe and set row names as a 'sample' column
df <- as.data.frame(cmd$points, stringsAsFactors = FALSE)
df$sample <- rownames(df)  # Set rownames as the 'sample' column

# Add metadata columns
df$pyrethroid_resistance_status <- gsub("_", " ", desc$pyrethroid_resistance_status)

# Rename PCA columns to start with 'PC' if they start with 'V'
colnames(df) <- gsub("^V", "PC", colnames(df))

color_by <- "pyrethroid_resistance_status" # specify if colored by region or country

## changing colour scheme to be with viridis, with color_by being a discrete variable

my_colours <- c("susceptible" = "#6f9c3d", "resistant" = "#fe6b40")

# plot
png("3_only_phenotyped_PCA_coloured_by_population_dpi_600.png", width = 7000, height = 7000, res = 600) 
ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by))) +
    geom_point(size = 5) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), y = paste0("PC2", " (", vars["PC2"], "%)"), title = "Chromosome 3", color = "Population") +
    scale_color_manual(values = my_colours, labels = c("susceptible" = "Pyrethroid susceptible", "resistant" = "Pyrethroid resistant")) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(labels = label_number()) +
    theme_classic(base_size = 120) +
    theme(legend.position = "bottom", 
    legend.direction = "horizontal", 
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 40, b = 30, l = 10, unit = "pt"),
    legend.text = element_text(size = 120),
    axis.line = element_line(linewidth = 0.5),
    axis.ticks = element_line(linewidth = 0.5))
dev.off()
#
