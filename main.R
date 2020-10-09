# Import packages
library(limma)
library(edgeR)
library(magrittr)
library(readr)
library(SummarizedExperiment)
library(DESeq2)
library(org.Mm.eg.db)
library(biomaRt)
library(Glimma)
library(GO.db)
library(EGSEA)
library(EGSEAdata)
library(dplyr)
library(RColorBrewer)
library(gplots)
library(ggplot2)

# Wrangle RNASeq-count datafiles into single dataframe
heart_count_dataframe <- readr::read_tsv(file.path(getwd(), "data", "jain",
                                                   "HeartKLF2-KLF4_DKO_K2K4_DKO_vs_Cre_rowcount_data.tsv"))
lung_count_dataframe <- readr::read_tsv(file.path(getwd(), "data", "jain",
                                                  "KLF2-KLF4_DKO_Lungs_DKO_vs_CRE_rowcount_data.tsv" ))
merge_count_dataframe <- dplyr::inner_join(heart_count_dataframe,
                                           lung_count_dataframe,
                                           by = "ensgid")
merge_count_dataframe <- dplyr::rename(merge_count_dataframe,
                                       ensembl_id = ensgid, 
                                       h.cre.ctrl.s1 = H_Cre_Ctrl_1,
                                       h.cre.ctrl.s2 = H_Cre_Ctrl_2,
                                       h.cre.ctrl.s3 = H_Cre_Ctrl_3,
                                       h.cre.ctrl.s4 = H_Cre_Ctrl_4,
                                       h.k2k4.dko.s1 = H_K2K4_DKO_1,
                                       h.k2k4.dko.s2 = H_K2K4_DKO_2,
                                       h.k2k4.dko.s3 = H_K2K4_DKO_3,
                                       h.k2k4.dko.s4 = H_K2K4_DKO_4,
                                       l.cre.ctrl.s1 = CRE_1_S1,
                                       l.cre.ctrl.s2 = CRE_2_S2,
                                       l.cre.ctrl.s3 = CRE_3_S3,
                                       l.cre.ctrl.s4 = CRE_4_S4,
                                       l.k2k4.dko.s1 = DKO_1_S5,
                                       l.k2k4.dko.s2 = DKO_2_S6,
                                       l.k2k4.dko.s3 = DKO_3_S7,
                                       l.k2k4.dko.s4 = DKO_4_S8)

# Create dataframe of sample information
sample_name <- colnames(merge_count_dataframe)[-1]
tissue_origin <- substring(sample_names,1,1) %>%
  forcats::as_factor()
cell_genotype <- rep(c("ctrl", "dko"), each = 4, times = 2) %>%
  forcats::as_factor()
sample_info <- dplyr::tibble(sample_names, tissue_origin, cell_genotype)
group <- paste(sample_info$tissue_origin,
               sample_info$cell_genotype,
               sep=".") %>%
  forcats::as_factor()

# Convert count dataframe to matrix
ensembl_id <- merge_count_dataframe$ensembl_id
merge_count_dataframe$ensembl_id <- NULL
count_matrix <- as.matrix(merge_count_dataframe)
rownames(count_matrix) <- ensembl_id

# Filter to remove lowly expressed genes
myCPM <- cpm(count_matrix)
## Ensure CPM threshold corresponds to a count of 10-15 for the library sizes
plot(myCPM[,1],count_matrix[,1],ylim=c(0,50),xlim=c(0,3)) # looking at first sample
abline(v=0.5)
thresh <- myCPM > 0.5
keep <- rowSums(thresh) >= 2 # keep genes that have at least 2 TRUES in each row of thresh
count_keep <- count_matrix[keep,]

# Create DGEList object
count_dge <- edgeR::DGEList(counts = count_keep)
## Plot the library sizes as a barplot
barplot(count_dge$samples$lib.size,names=colnames(count_dge),las=2) # not normally distributed
title("Barplot of library sizes")
## Examine distribution on log2 scale
logcounts <- cpm(count_dge,log=TRUE)
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2)
abline(h=median(logcounts),col="blue")
title("Boxplots of logCPMs (unnormalised)")

# Principle components analysis with MDS
par(mfrow=c(1,2))
col_tissue <- c("purple","orange")[sample_info$tissue_origin]
col_gene <- c("blue","dark green")[sample_info$cell_genotype]
plotMDS(count_dge, col = col_tissue, pch=16,cex=2)
legend("top",fill=c("purple","orange"),legend=levels(sample_info$tissue_origin),cex=0.8)
title("Tissue Origin")
plotMDS(count_dge, col = col_gene, pch=16,cex=2)
legend("top",fill=c("blue","dark green"),legend=levels(sample_info$cell_genotype),cex=0.8)
title("Cell Genotype")
par(mfrow=c(1,1))
char_gene <- c(1,4)[sample_info$cell_genotype]
cols <- c("purple","orange")
plotMDS(count_dge,dim=c(1,2),col=col_tissue,pch=char_gene,cex=2)
legend("top",legend=c("heart","lung"),col=cols,pch=16)
legend("bottom",legend=c("control","knock-out"),pch=c(1,4))
title("PCA Bi-plot")

# Heirarchical clustering with heatmaps
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing=TRUE))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
mypalette <- RColorBrewer::brewer.pal(11,"RdYlBu")
morecols <- colorRampPalette(mypalette)
png(file="High_var_genes.heatmap.png")
gplots::heatmap.2(highly_variable_lcpm,
          col=rev(morecols(50)),
          trace="none",
          main="Top 500 most variable genes across samples",
          ColSideColors=col_tissue,
          scale="row")
dev.off()

# Apply normalisation to DGEList object
count_dge <- calcNormFactors(count_dge)
## Mean difference plots, comparing log counts to normalised counts
par(mfrow=c(1,2))
plotMD(logcounts,column = 7)
abline(h=0,col="grey")
plotMD(logcounts,column = 10)
abline(h=0,col="grey")
par(mfrow=c(1,2))
plotMD(count_dge,column = 7)
abline(h=0,col="grey")
plotMD(count_dge,column = 10)
abline(h=0,col="grey")

#save(group,count_dge,count_keep,logcounts,sample_info,file="rnaseq_objects.Rdata")
