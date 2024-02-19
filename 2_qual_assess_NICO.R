### STEP 2: Quality assessment and visualization ###

## This script will prepare normalized gene counts.
## PCA plot will also be generated for quality assessment.

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(vsn)
  library(RColorBrewer)
})

# Load .RData file from STEP 1
load("data/NICO_DESeqObject.RData")


# Count data transformations ----------------------------------------------

vsd <- vst(dds_regions, blind = FALSE) # variance stabilizing transformation
ntd <- normTransform(dds_regions, f = log2, pc = 1) # normalized counts transformation
dds <- estimateSizeFactors(dds_regions) # account for sequencing depth

# Plot transformations from 'core' samples in Patient 1 and 2
df <- bind_rows(
  as_tibble(assay(ntd)[, c(3,6)]) %>%
    mutate(transformation = "log2(x + 1)"),
  as_tibble(assay(vsd)[, c(3,6)]) %>%
    mutate(transformation = "vst")
)

colnames(df)[1:2] <- c("x", "y")
ggplot(df, aes(x = x, y = y)) + geom_hex(bins = 80) +
  coord_fixed() + facet_grid(. ~ transformation)

# Visualize effect of transformation on variance
meanSdPlot(assay(ntd))
meanSdPlot(assay(vsd))

# Data quality assessment by PCA -----------------

# Principal component plot of the samples
pcaData <- plotPCA(vsd, intgroup = c("description"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
pcaPlot <- ggplot(pcaData, aes(x = PC1, y = PC2, color = description)) +
  geom_point(size = 3) +
  stat_ellipse() +
  scale_color_brewer(palette = "Set1", labels = c("Core", "CE", "NCE")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) + 
  coord_fixed() +
  theme_classic() + theme(legend.title = element_blank(),
                          legend.position = c(0.9,0.15))
pcaPlot

# Save transformed counts into .RData file (replace dds_regions with dds)
save(vsd, ntd, dds,
     file = "data/NICO_DESeqObject.RData")
