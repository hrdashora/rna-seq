## This script will use the RNA-seq normalized counts table generated
## by the Xie Lab to visualize expression levels of SWI/SNF complex genes 

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
})

# Create list of chromatin remodeler complex genes
crc_genes <- c("SMARCC1", "SMARCC2", "SMARCD1", # core subunits of mammalian SWI/SNF complex families
          "SMARCA4", "SMARCA2", # ATPases of mammalian SWI/SNF complexes
          "SMARCA5", "SMARCA1", # ATPases of ISWI complexes
          # non-catalytic subunits of IWSI complexes, associatitng with SMARCA5
          "BAZ1A", "BAZ1B", "BAZ2A", "BAZ2B", "RSF1", "BPTF", "CECR2", "POLE3", "CHRAC1",
          "CHD4") # core component of NuRD complex

# Read data from .xlsx
normCounts_crc <- readxl::read_xlsx("data/GSC_RNAseq_normalized_count_table.xlsx") %>%
  filter(Gene %in% crc_genes) %>%
  tidyr::pivot_longer(cols = where(is.numeric), names_to = "Sample", values_to = "Normalized_Counts") %>%
  filter(stringr::str_detect(Gene, "SMARCA"))

# Plot normalized counts
ggpubr::ggbarplot(
  data = normCounts_crc,
  x = "Sample",
  y = "Normalized_Counts",
  facet.by = "Gene",
  fill = "Gene",
  color = "Gene",
  palette = "npg",
  xlab = "Cell Line",
  ylab = "Normalized Count",
  title = "RNA-seq gene expression levels in GSCs",
  subtitle = "SWI/SNF & ISWI complex ATPases",
  caption = "Xie Lab",
  position = position_dodge(width = 0.8)
) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  ggpubr::rotate_x_text(angle = 90) +
  ggpubr::rremove("legend") +
  ggpubr::rremove("xlab")
