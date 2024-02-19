### STEP 4: Plot significant DE genes ###

## This script imports the results CSV files generated in STEP 3 and creates
## volcano plots.

# Load libraries
suppressPackageStartupMessages({
  library(EnhancedVolcano)
})

# Load CSV and remove genes by automatic independent filtering for having
# a low mean normalized count, or an extreme count outlier
resdf_regCEvC <- readr::read_csv(file = "results/description_CE_versus_core.csv") %>%
  dplyr::filter(!is.na(padj))
resdf_regCEvNCE <- readr::read_csv(file = "results/description_CE_versus_NCE.csv") %>%
  dplyr::filter(!is.na(padj))
resdf_regCvNCE <- readr::read_csv(file = "results/description_core_versus_NCE.csv") %>%
  dplyr::filter(!is.na(padj))
resdf_lrt <- readr::read_csv(file = "results/LRT_results.csv") %>%
  dplyr::filter(!is.na(padj))

# Create volcano plots
## NOTE: 'CE' versus 'core' plot is uninformative
volc_regCvNCE <- EnhancedVolcano(toptable = resdf_regCvNCE,
                                 lab = resdf_regCvNCE$GeneSymbol,
                                 x = 'log2FoldChange',
                                 y = 'padj',
                                 title = "Volcano plot of Resected GBM Tissue",
                                 subtitle = bquote(italic("DEGs in Core Region (+) versus NCE Region (-)")),
                                 pCutoff = 0.05,
                                 FCcutoff = 1
                                 )
volc_regCvNCE

volc_regCEvNCE <- EnhancedVolcano(toptable = resdf_regCEvNCE,
                                  lab = resdf_regCEvNCE$GeneSymbol,
                                  x = 'log2FoldChange',
                                  y = 'padj',
                                  title = "Volcano plot of Resected GBM Tissue",
                                  subtitle = bquote(italic("DEGs in CE Region (+) versus NCE Region (-)")),
                                  pCutoff = 0.05,
                                  FCcutoff = 1
                                  )
volc_regCEvNCE

volc_lrt <- EnhancedVolcano(toptable = resdf_lrt,
                            lab = resdf_lrt$GeneSymbol,
                            x = 'log2FoldChange',
                            y = 'padj')
volc_lrt


