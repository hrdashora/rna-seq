## Load bioconductor packages
# library(AllelicImbalance)
# library(bumphunter) 
# library(csaw)
# library(DESeq)
library(edgeR)
# library(IRanges)
# library(Rsamtools)
# library(rtracklayer)
# library(sva)
# library(SummarizedExperiment)
# library(VariantAnnotation)

## Load R packages
# library(dplyr)
# library(extRemes)
# library(forcats)
library(magrittr)
# library(remotes) # powsimR
library(readr)

## Load count data
count_dataframe <- readr::read_tsv(file.path(getwd(), "data", "cookbook", "ch1", "modencodefly_count_table.txt" ))
genes <- count_dataframe[['gene']]
count_dataframe[['gene']] <- NULL
count_matrix <- as.matrix(count_dataframe)
rownames(count_matrix) <- genes
pheno_data <- readr::read_table2(file.path(getwd(), "data", "cookbook", "ch1", "modencodefly_phenodata.txt"))

## Experiments of interest
experiments_of_interest <- c("L1Larvae", "L2Larvae")
columns_of_interest <- which( pheno_data[['stage']] %in% experiments_of_interest ) 

## Form the grouping factor
grouping <- pheno_data[['stage']][columns_of_interest] %>% 
  forcats::as_factor()

## Form subset of count data
counts_of_interest <-  count_matrix[,columns_of_interest]

## Create DGE object
count_dge <- edgeR::DGEList(counts = counts_of_interest, group = grouping)

## Perform differential expression analysis
design <- model.matrix(~ grouping)
eset_dge <- edgeR::estimateDisp(count_dge, design)
fit <- edgeR::glmQLFit(eset_dge, design)
result <- edgeR::glmQLFTest(fit, coef=2)
topTags(result)
