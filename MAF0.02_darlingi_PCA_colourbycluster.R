# Run this as an R script

library(showtext)
library(dplyr)
library(ggplot2)
library(ape)
showtext_auto()
library(viridis)
library(scales)

workdir <- "/mnt/storage11/sophie/darlingi/darlingi_database" # Working directory with plink files
prefix <- "mito_only_darlingi" # Prefix for plink files
metadata <- "/mnt/storage11/sophie/darlingi/darlingi_metadata.csv" # File path to metadata

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
#df$sample <- gsub("_", " ", desc$sample)
df$population <- gsub("_", " ", desc$population)

# Correctly reorder the dataframe to make 'sample' the first column
#df <- df[, c("sample", setdiff(names(df), "sample"))]

# Rename PCA columns to start with 'PC', if they start with 'V' 
# (assuming all PCA columns currently start with 'V', adjust if different)
names(df) <- gsub("^V", "PC", names(df))

color_by <- "population" # specify if colored by region or country

## changing colour scheme to be with viridis, with color_by being a discrete variable

my_colours <- c("Rondonia State" = "#FFD23F", "Colony old" = "#EE4266", 
                "colony new" = "#3BCEAC")

# plot
png("darlingi_MAF_0.02_mito_only_PCA_coloured_by_cluster.png") 
ggplot(data = df, aes(x = PC1, y = PC2, color = !!sym(color_by))) +
    geom_point(size = 3) +
    labs(x = paste0("PC1", " (", vars["PC1"], "%)"), y = paste0("PC2", " (", vars["PC2"], "%)"), title = "Mitochondrial genome") +
    scale_color_manual(values = my_colours) +
    scale_x_continuous(labels = label_number()) +
    scale_y_continuous(labels = label_number()) +
    theme_classic() +
    theme(legend.position = "bottom", legend.direction = "vertical", plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 10, r = 40, b = 30, l = 10, unit = "pt"))
dev.off()
#
