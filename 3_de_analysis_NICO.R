### STEP 3: Running analysis and exporting DESeq2 results ###

# Load libraries
suppressPackageStartupMessages({
  library(DESeq2)
  library(DEGreport)
})

library(AnnotationDbi)
library(org.Hs.eg.db)
library(sva)
library(DESeq2)
library(ggplot2)
library(readr)
library(tibble)
library(dplyr)

# Load .RData file from STEP 2
load("data/NICO_DESeqObject.RData")

# Run differential expression analysis
dds <- DESeq(dds)

# Build results tables
## Apply fold change shrinkage to remove the noise associated with LFC from
## low-count genes - useful for ranking and visualization.
resultsNames(dds)
resLFC_regCvNE <- lfcShrink(dds,
                            coef = "description_core_vs_non_enhancing",
                            type = "apeglm",
                            saveCols = c(1,2,3))
resLFC_regEvNE <- lfcShrink(dds,
                            coef = "description_enhancing_vs_non_enhancing",
                            type = "apeglm",
                            saveCols = c(1,2,3))
resLFC_regEvC <- lfcShrink(dds,
                           contrast = c("description","enhancing","core"),
                           type = "ashr",
                           saveCols = c(1,2,3)) # utilizing a different shrinkage estimator

# View information on each column in results tables
mcols(resLFC_regCvNE, use.names = TRUE)

# Identify significantly different genes across all contrasts

log2cutoff <- 1.5
qvaluecutoff <- 0.05

sigGenes <- unique(c(
  subset(resLFC_regCvNE, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$EnsemblID,
  subset(resLFC_regEvNE, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$EnsemblID,
  subset(resLFC_regEvC, padj<=qvaluecutoff & abs(log2FoldChange)>=log2cutoff)$EnsemblID
))

length(sigGenes)

# Perform Likelihood ratio test to test all levels of the 'description' factor at once
## Determine if the variation in the data is explained significantly using the
## full model compared to the reduced model.
dds_lrt <- DESeq(dds, test = "LRT", reduced = ~patient)
resultsNames(dds_lrt)
res_lrt <- results(dds_lrt, saveCols = c(1,2,3)) # extract results

## P-values are determined by the difference in deviance between the full and
## reduced model formula (not log2 fold changes). Fold changes are not directly
## associated with the actual hypothesis test.

# Subset the LRT results to return genes with padj < 0.05

sig_res_lrt <- res_lrt %>%
  data.frame() %>%
  as_tibble() %>% # row names removed
  filter(padj < qvaluecutoff)

# Get sig gene lists
sigGenes_lrt <- sig_res_lrt %>% 
  pull(EnsemblID)

length(sigGenes_lrt)

# Subset results for faster cluster finding
clustering_sig_genes <- sig_res_lrt %>%
  arrange(padj) %>%
  head(n = 2000)

# Obtain vsd values for those significant genes
cluster_vsd <- vsd[clustering_sig_genes$EnsemblID, ]

# Use the `degPatterns` function from the 'DEGreport' package to show gene clusters across sample groups
clusters <- degPatterns(ma = assay(cluster_vsd),
                        metadata = colData(cluster_vsd),
                        time = "description",
                        col = NULL)

# Extract the groups of genes that share a pattern of expression change
# across the factor levels
cluster_groups <- clusters$df
group1 <- clusters$df %>%
  filter(cluster == 1)
group2 <- clusters$df %>%
  filter(cluster == 2)

# Exporting results to CSV files
resLFC_regCvNE %>%
  as_tibble() %>%
  write_csv(file = "results/description_core_versus_NCE.csv")

resLFC_regEvNE %>%
  as_tibble() %>%
  write_csv(file = "results/description_CE_versus_NCE.csv")

resLFC_regEvC %>%
  as_tibble() %>%
  write_csv(file = "results/description_CE_versus_core.csv")

res_lrt %>%
  as_tibble() %>%
  write_csv(file = "results/LRT_results.csv")

# Also save results into .RData file
save(vsd, ntd, dds, dds_lrt,
     resLFC_regCvNE, resLFC_regEvNE, resLFC_regEvC, res_lrt,
     file = "data/NICO_DESeqObject.RData")
