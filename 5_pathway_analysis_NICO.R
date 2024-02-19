### STEP 5: Gene Set Enrichment and Over-representation Analysis with MSigDB and GO Terms ###

## Significant genes are those with adjusted p-value < .1
## and the universe is the set of genes with adjusted p-value not NA

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(goseq)
  library(clusterProfiler)
  library(ReactomePA)
})

# Load .RData file from STEP 3
load("NICO_DESeqObject.RData")

# goseq -------------------------------------------------------------------

# goseq requires a simple named vector, which contains two pieces of information.
# 1. Measured genes: all genes for which RNA-seq data was gathered.
# Each element of your vector should be named by a unique gene identifier.
# 2. Differentially expressed genes: each element of your vector should be
# either a 1 or a 0, where 1 indicates that the gene is differentially expressed
# and 0 that it is not.

rowmean.threshold <- 1
fdr.threshold <- 0.1
rs <- rowMeans(counts(ddssva))
dds <- ddssva[ rs > rowmean.threshold, ]
dds <- DESeq(dds)
res_CvNCE <- results(dds, contrast = c("description","core","non_enhancing"),
                     independentFiltering = FALSE) # using mean counts threshold instead of IF
res_CEvNCE <- results(dds, contrast = c("description","enhancing","non_enhancing"),
                     independentFiltering = FALSE)
res_CEvC
assayed.genes_CvNCE <- rownames(res_CvNCE)
de.genes_CvNCE <- rownames(res_CvNCE)[ which(res_CvNCE$padj < fdr.threshold) ]
assayed.genes_CEvNCE <- rownames(res_CEvNCE)
de.genes_CEvNCE <- rownames(res_CEvNCE)[ which(res_CEvNCE$padj < fdr.threshold) ]

# Construct named vectors suitable for goseq
gene.vector_CvNCE <- as.integer(assayed.genes_CvNCE %in% de.genes_CvNCE)
names(gene.vector_CvNCE) <- assayed.genes_CvNCE
head(gene.vector_CvNCE)

gene.vector_CEvNCE <- as.integer(assayed.genes_CEvNCE %in% de.genes_CEvNCE)
names(gene.vector_CEvNCE) <- assayed.genes_CEvNCE
head(gene.vector_CEvNCE)

# Check which genome/gene identifier combinations are in the local database
supportedOrganisms()[supportedOrganisms()$Genome == "hg38",]
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txsByGene <- transcriptsBy(txdb, "gene")
lengthData <- median(width(txsByGene))

# Fitting the PWF
pwf_CvNCE <- nullp(gene.vector_CvNCE, genome = "hg19", id = "ensGene")
head(pwf_CvNCE)
pwf_CEvNCE <- nullp(gene.vector_CEvNCE, genome = "hg19", id = "ensGene")
head(pwf_CEvNCE)

# Using the Wallenius approximation
GO.wall_CvNCE <- goseq(pwf_CvNCE, "hg19", "ensGene")
head(GO.wall_CvNCE)
enriched.GO_CvNCE <- GO.wall_CvNCE$category[p.adjust(GO.wall_CvNCE$over_represented_pvalue,
                                                     method = "BH") < .05]
head(enriched.GO_CvNCE)

GO.wall_CEvNCE <- goseq(pwf_CEvNCE, "hg19", "ensGene")
head(GO.wall_CEvNCE)
enriched.GO_CEvNCE <- GO.wall_CEvNCE$category[p.adjust(GO.wall_CEvNCE$over_represented_pvalue,
                                                     method = "BH") < .05]
head(enriched.GO_CEvNCE)

# Access information top 10 GO terms
for(go in enriched.GO_CvNCE[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

for(go in enriched.GO_CEvNCE[1:10]){
  print(GOTERM[[go]])
  cat("--------------------------------------\n")
}

# Plot Top 10 GO Terms
df_CvNCE <- GO.wall_CvNCE %>%
  dplyr::filter(ontology == "BP") %>%
  dplyr::arrange(over_represented_pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(enrichment_score = -log10(over_represented_pvalue))
  
p_CvNCE <- ggbarplot(df_CvNCE, x = "term", y = "enrichment_score",
               color = "white",
               fill = "#E64B35FF",
               sort.val = "asc",
               subtitle = "Core Region versus NCE Region",
               ylab = "Enrichment Score",
               xlab = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  rotate() +
  rremove("legend")

p_CvNCE

df_CEvNCE <- GO.wall_CEvNCE %>%
  dplyr::filter(ontology == "BP") %>%
  dplyr::arrange(over_represented_pvalue) %>%
  dplyr::slice_head(n = 10) %>%
  dplyr::mutate(enrichment_score = -log10(over_represented_pvalue))

p_CEvNCE <- ggbarplot(df_CEvNCE, x = "term", y = "enrichment_score",
                     color = "white",
                     fill = "#4DBBD5FF",
                     sort.val = "asc",
                     subtitle = "CE Region versus NCE Region",
                     ylab = "Enrichment Score",
                     xlab = "") +
  scale_y_continuous(expand = expansion(mult = c(0, .1))) +
  rotate() +
  rremove("legend")

p_CEvNCE


# clusterProfiler ---------------------------------------------------------

# Prepare geneLists for clusterProfiler
# Load CSVs of DESeq results and remove genes by automatic independent
# filtering for having a low mean normalized count, or an extreme count outlier
resdf_regCEvNCE <- readr::read_csv(file = "results/description_CE_versus_NCE.csv") %>%
  dplyr::filter(!is.na(padj))
resdf_regCvNCE <- readr::read_csv(file = "results/description_core_versus_NCE.csv") %>%
  dplyr::filter(!is.na(padj))

## feature 1: numeric vector
geneList_CEvNCE <- resdf_regCEvNCE %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(log2FoldChange)
geneList_CvNCE <- resdf_regCvNCE %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(log2FoldChange)
geneList <- reslfc %>% as.data.frame() %>% dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(log2FoldChange)

## feature 2: named vector
names(geneList_CEvNCE) <- resdf_regCEvNCE %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(EntrezID)
names(geneList_CvNCE) <- resdf_regCvNCE %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(EntrezID)
names(geneList) <- reslfc %>% as.data.frame() %>% dplyr::filter(!is.na(padj)) %>%
  dplyr::filter(!is.na(EntrezID)) %>%
  dplyr::pull(EntrezID)

# Remove duplicate gene names
keep <- which(!duplicated(names(geneList_CEvNCE)))
geneList_CEvNCE <- geneList_CEvNCE[keep]

keep <- which(!duplicated(names(geneList_CvNCE)))
geneList_CvNCE <- geneList_CvNCE[keep]

keep <- which(!duplicated(names(geneList)))
geneList <- geneList[keep]

## feature 3: decreasing order
geneList_CEvNCE <- sort(geneList_CEvNCE, decreasing = TRUE)
geneList_CvNCE <- sort(geneList_CvNCE, decreasing = TRUE)
geneList <- sort(geneList, decreasing = TRUE)

# GO Gene Set Enrichment Analysis
ego_CEvNCE <- gseGO(geneList     = geneList_CEvNCE,
                    ont          = "ALL",
                    OrgDb        = org.Hs.eg.db,
                    keyType      = "ENTREZID",
                    minGSSize    = 10,
                    maxGSSize    = 500,
                    pvalueCutoff = 0.1,
                    pAdjustMethod = "BH",
                    verbose      = TRUE)

ego_CvNCE <- gseGO(geneList     = geneList_CvNCE,
                   ont          = "ALL",
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   verbose      = TRUE)

ego <- gseGO(geneList     = geneList,
                   ont          = "ALL",
                   OrgDb        = org.Hs.eg.db,
                   keyType      = "ENTREZID",
                   minGSSize    = 10,
                   maxGSSize    = 500,
                   pvalueCutoff = 0.1,
                   pAdjustMethod = "BH",
                   verbose      = TRUE)

# Isolate BP ontology and reduce redundancy of GO terms
bpego_CEvNCE <- ego_CEvNCE %>% dplyr::filter(ONTOLOGY == "BP")
bpego_CvNCE <- ego_CvNCE %>% dplyr::filter(ONTOLOGY == "BP")
bpego <- ego %>% dplyr::filter(ONTOLOGY == "BP")

bpego_CEvNCE <- simplify(bpego_CEvNCE, cutoff = 0.7, by = "p.adjust", select_fun = min)
bpego_CvNCE <- simplify(bpego_CvNCE, cutoff = 0.7, by = "p.adjust", select_fun = min)
bpego <- simplify(bpego, cutoff = 0.7, by = "p.adjust", select_fun = min)

bpego_CEvNCE <- enrichplot::pairwise_termsim(bpego_CEvNCE)
bpego_CvNCE <- enrichplot::pairwise_termsim(bpego_CvNCE)
bpego <- enrichplot::pairwise_termsim(bpego)

bpegox_CEvNCE <- setReadable(bpego_CEvNCE, 'org.Hs.eg.db', 'ENTREZID')
bpegox_CvNCE <- setReadable(bpego_CvNCE, 'org.Hs.eg.db', 'ENTREZID')
bpegox <- setReadable(bpego, 'org.Hs.eg.db', 'ENTREZID')

geneInCategory(bpegox_CEvNCE)
geneInCategory(bpegox_CvNCE)
geneInCategory(bpegox)

# Visualize GO enrichment results
## Dot plot
label_names <- c(activated = "Activated in CE", suppressed = "Suppressed in CE")
dotplot(bpego_CEvNCE, showCategory = 5, label_format = 60, font.size = 12, split = ".sign",
        title = "BP GO Terms for DEGs between CE and NCE Regions") +
  facet_grid(.~.sign, labeller = as_labeller(label_names))
label_names <- c(activated = "Activated in Core", suppressed = "Suppressed in Core")
dotplot(bpego_CvNCE, showCategory = 5, label_format = 60, font.size = 12, split = ".sign",
        title = "BP GO Terms for DEGs between Core and NCE Regions") +
  facet_grid(.~.sign, labeller = as_labeller(label_names))

label_names <- c(activated = "NER", suppressed = "Core/CE")
dotplot(bpego, showCategory = 5, label_format = 60, font.size = 12, split = ".sign",
        title = "BP GO Terms for DEGs") +
  facet_grid(.~.sign, labeller = as_labeller(label_names))
library(enrichplot)
library(clusterProfiler)

bpegox %>% filter(NES >= 1) %>%
  mutate(qscore = -log(p.adjust, base = 10)) %>%
  dotplot(x = "qscore")

bpegox %>% filter(NES <= -1) %>%
  mutate(qscore = -log(p.adjust, base = 10)) %>%
  dotplot(x = "qscore")

## Gene-concept network
cnetplot(bpegox_CEvNCE, categorySize = "p.adjust", foldChange = geneList_CEvNCE)
cnetplot(bpegox_CvNCE, categorySize = "p.adjust", foldChange = geneList_CvNCE)

## Heatmap-like functional classification
heatplot(bpegox_CEvNCE, foldChange = geneList_CEvNCE, showCategory = 5) +
  ggtitle("heatplot for GSEA")
heatplot(bpegox_CvNCE, foldChange = geneList_CvNCE, showCategory = 5) +
  ggtitle("heatplot for GSEA")
## Tree plot
treeplot(bpego_CEvNCE)
treeplot(bpego_CvNCE)
## Enrichment map
emapplot(bpego_CEvNCE)
emapplot(bpego_CvNCE)
## Upset plot
upsetplot(bpego_CEvNCE)
upsetplot(bpego_CvNCE)
upsetplot(bpego)
## Ridgeline plot
ridgeplot(bpego_CEvNCE, label_format = 60) + ggtitle("ridgeplot for GSEA")
ridgeplot(bpego_CvNCE, label_format = 60) + ggtitle("ridgeplot for GSEA")
ridgeplot(bpego)
# Save plots
# ggsave(filename = "", plot = NULL, device = png,
#        width = 30, height = 30, units = "cm", dpi = "retina")
# ggsave(filename = "", plot = NULL, device = cairo_ps,
#        width = 30, height = 30, units = "cm",
#        dpi = "retina", fallback_resolution = 320)

# Reactome pathway gene set enrichment analysis
y_CEvNCE <- gsePathway(geneList_CEvNCE,
                       pvalueCutoff = 0.2,
                       pAdjustMethod = "BH",
                       verbose = TRUE)
head(y_CEvNCE)

label_names <- c(activated = "Activated in CE", suppressed = "Suppressed in CE")
dotplot(y_CEvNCE, showCategory = 5, label_format = 60, font.size = 12, split = ".sign",
        title = "Reactome Terms for DEGs between CE and NCE Regions") +
  facet_grid(.~.sign, labeller = as_labeller(label_names))

y_CvNCE <- gsePathway(geneList_CvNCE,
                       pvalueCutoff = 0.2,
                       pAdjustMethod = "BH",
                       verbose = TRUE)
head(y_CvNCE)

label_names <- c(activated = "Activated in Core", suppressed = "Suppressed in Core")
dotplot(y_CvNCE, showCategory = 5, label_format = 60, font.size = 12, split = ".sign",
        title = "Reactome Terms for DEGs between Core and NCE Regions") +
  facet_grid(.~.sign, labeller = as_labeller(label_names))

# Pathway visualization
viewPathway("E2F mediated regulation of DNA replication", 
            readable = TRUE, 
            foldChange = geneList_CEvNCE)



# MSigDB ------------------------------------------------------------------
library(msigdb)
library(ExperimentHub)
library(GSEABase)

# Download data from the msigdb R package
eh = ExperimentHub()
query(ExperimentHub(), 'msigdb')
msigdb.hs <- getMsigdb(org = "hs", id = "EZID", version = "7.5")
msigdb.hs <- appendKEGG(msigdb.hs, version = "7.5")
length(msigdb.hs)
#calculate the number of signatures in each category
table(sapply(lapply(msigdb.hs, collectionType), bcCategory))
#calculate the number of signatures in each subcategory
table(sapply(lapply(msigdb.hs, collectionType), bcSubCategory))

# Subset collections from the MSigDB
c <- subsetCollection(msigdb.hs, collection = c('c2','c5'))
c_ids <- geneIds(c)
c2 <- subsetCollection(msigdb.hs, 'c2')
c2_ids <- geneIds(c2)

# Create expression data using vst data from DESeq2
entrez.genes <- rowData(dds_intext, use.names = FALSE)$EntrezID
exprs_mat <- counts(dds_intext, normalized = TRUE)

rownames(exprs_mat) <- entrez.genes
colnames(exprs_mat) <- rownames(colData(dds_intext))
exprs_mat <- exprs_mat[!is.na(rownames(exprs_mat)), ]
head(exprs_mat)

# library(biomaRt)
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# genes <- getBM(
#   filters = "ensembl_gene_id",
#   attributes = c("ensembl_gene_id", "entrezgene_id"),
#   values = ensembl.genes,
#   mart = mart)
# convert.genes <- genes %>% na.omit() %>% distinct(ensembl_gene_id, .keep_all = TRUE)
# keep <- intersect(ensembl.genes, convert.genes$ensembl_gene_id)

# Convert gene sets into a list of gene indices
c_indices <- limma::ids2indices(c2_ids, rownames(exprs_mat))

# Create a design matrix and contrasts
description <- colData(vsd_reg)$description
subject <- colData(vsd_reg)$subject
design <- model.matrix(~0 + description)
design
contr.matrix <- limma::makeContrasts(
  CvNCE = descriptioncore - descriptionnon_enhancing,
  CEvNCE = descriptionenhancing - descriptionnon_enhancing,
  CvCE = descriptioncore - descriptionenhancing,
  levels = colnames(design))
contr.matrix

# Apply the romer method to the gene sets, COMPUTATIONALLY INTENSIVE
c_cam.CvNCE <- limma::romer(y = exprs_mat,
                            index = c_indices,
                            design = design,
                            contrast = contr.matrix[,1])
head(c_cam.CvNCE, 5)

c_cam.CEvNCE <- limma::romer(y = exprs_mat,
                            index = c_indices,
                            design = design,
                            contrast = contr.matrix[,2])
head(c_cam.CEvNCE, 5)

# Make a barcodeplot for particular gene signatures
vfit <- limma::lmFit(exprs_mat, design)

vfit_CvNE <- limma::contrasts.fit(vfit, contrasts = contr.matrix[,"CvNCE"])
efit_CvNE <- limma::eBayes(vfit_CvNE)
rslt_CvNE <- limma::decideTests(vfit_CvNE)

vfit_CEvNCE <- limma::contrasts.fit(vfit, contrasts = contr.matrix[,"CEvNCE"])
efit_CEvNCE <- limma::eBayes(vfit_CEvNCE)
rslt_CEvNCE <- limma::decideTests(vfit_CEvNCE)

limma::barcodeplot(efit_CvNE$t[,1],
                   index = c_indices$BUFFA_HYPOXIA_METAGENE)

# Use GSVA for single sample gene set enrichment, starting with gene expression data matrix
exprs_mat[1:5, 1:5] # vst gene expression matrix

# Calculate GSVA enrichment scores
library(GSVA)
gsva.es <- gsva(expr = exprs_mat, # provide matrix of expression values
                gset.idx.list = c2, # provide GeneSetCollection object
                method = "gsva",
                kcdf = "Gaussian", # input expression values are log transformed
                verbose = T)

# Calculate enrichment scores for MSigDB lists
# Select relevant sets
library(msigdbr)
library(enrichplot)
gene_sets <- msigdbr(species = "human")
head(gene_sets)
cgp_gene_sets <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
head(cgp_gene_sets)
hall_gene_sets <- msigdbr(species = "human", category = "H")
unique(hall_gene_sets$gs_name)

msigdbr_t2g <- cgp_gene_sets %>%
  dplyr::distinct(gs_name, entrez_gene)
x_cgp <- GSEA(geneList, TERM2GENE = msigdbr_t2g)

msigdbr_t2g <- hall_gene_sets %>%
  dplyr::distinct(gs_name, entrez_gene)
x_hall <- GSEA(geneList, TERM2GENE = msigdbr_t2g)

y_cgp <- setReadable(x_cgp, 'org.Hs.eg.db', 'ENTREZID')
y_hall <- setReadable(x_hall, 'org.Hs.eg.db', 'ENTREZID')

# gseaplot2(y_cgp, geneSetID = "LEIN_ASTROCYTE_MARKERS",
#           base_size = 11,
#           pvalue_table = F)
# ggsave(filename = "gsea_NEvC.eps", plot = last_plot(), device = cairo_ps,
#        width = 25, height = 25, units = "cm",
#        dpi = "retina", fallback_resolution = 320)
gseaplot2(y, geneSetID = "HALLMARK_HYPOXIA",
          title = "Hypoxia Hallmark Genes in Cancer (MSigDB)",
          base_size = 11,
          pvalue_table = T)
gseaplot2(y, geneSetID = "LEIN_OLIGODENDROCYTE_MARKERS",
          title = "Oligodendrocyte Markers (Lein et al.)",
          base_size = 11,
          pvalue_table = F)
gseaplot2(y_cgp, geneSetID = c("VERHAAK_GLIOBLASTOMA_NEURAL",
                           "VERHAAK_GLIOBLASTOMA_MESENCHYMAL"))
gseaplot2(y, geneSetID = c("LEIN_OLIGODENDROCYTE_MARKERS", "BUFFA_HYPOXIA_METAGENE"))
cnetplot(y, foldChange = geneList_NEvC, showCategory = "LEIN_ASTROCYTE_MARKERS",
         layout = "kk",
         node_label = "all")
ggsave(filename = "net_NEvC.eps", plot = last_plot(), device = cairo_ps,
       width = 25, height = 25, units = "cm",
       dpi = "retina", fallback_resolution = 320)
heatplot(y, foldChange = geneList_CvNCE, showCategory = c("LEIN_ASTROCYTE_MARKERS"))
heatplot(y, foldChange = geneList_CvNCE, showCategory = 10)
dotplot(y_hall, showCategory  = 40, font.size = 8, label_format = 60)
dotplot(y, showCategory= 10, font.size = 8, label_format = 60, split=".sign") + facet_grid(.~.sign)

z <- mutate(egor_NEvC, richFactor = Count / as.numeric(sub("/\\d+", "", BgRatio)))
z

ggplot(z, showCategory = 20, 
       aes(richFactor, fct_reorder(Description, richFactor))) + 
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() + 
  xlab("rich factor") +
  ylab(NULL) + 
  ggtitle("Enriched MSigDB Curated Lists")

# Create heatmap with transformed counts data to visualize geneset
load("NICO_DESeqObjects.RData")
genes <- gene_sets %>%
  dplyr::filter(gs_name == "LEIN_ASTROCYTE_MARKERS") %>%
  pull(ensembl_gene)

genes_neural <- gene_sets %>%
  dplyr::filter(gs_name == c("VERHAAK_GLIOBLASTOMA_NEURAL")) %>%
  pull(ensembl_gene)
genes_mes <- gene_sets %>%
  dplyr::filter(gs_name == c("VERHAAK_GLIOBLASTOMA_MESENCHYMAL")) %>%
  pull(ensembl_gene)
genes <- c(genes_neural, genes_mes)


mat <- subset(counts(dds_intext, normalized = TRUE), rownames(counts(dds_intext, normalized = TRUE)) %in% genes)
mat  <- mat - rowMeans(mat) # amount by which each gene deviates in a specific sample from the gene’s average across all samples
pheatmap(mat)
anno <- as_tibble(colData(dds_intext), rownames = NA) %>% # need to juggle rownames in this manner to maintain pheatmap functionality
  rownames_to_column("sample") %>%
  dplyr::select(sample, description) %>%
  column_to_rownames("sample") %>%
  dplyr::filter(description %in% c("core","non_enhancing")) %>%
  arrange(description)
mat_order <- mat[ ,rownames(anno)] # reorder the matrix based on the annotation
library(RColorBrewer)
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(255)
library(pheatmap)
pheatmap(mat_order,
         cluster_rows = FALSE,
         show_rownames = F,
         cluster_cols = FALSE,
         annotation_col = anno,
         scale = "row",
         color = colors)

# Use GSVA to make a heatmap of pathways, not genes
mat <- gsva.es - rowMeans(gsva.es)
anno <- as_tibble(colData(ntd_reg), rownames = NA) %>% # need to juggle rownames in this manner to maintain pheatmap functionality
  rownames_to_column("sample") %>%
  dplyr::select(sample, description) %>%
  column_to_rownames("sample") %>%
  dplyr::filter(description %in% c("core","non_enhancing")) %>%
  arrange(description)
mat_order <- mat[ ,rownames(anno)] # reorder the matrix based on the annotation
colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(255)
pheatmap(mat_order, cluster_rows = FALSE, show_rownames = F,
         cluster_cols = FALSE, annotation_col = anno, color = colors)

msigdbr_list <- split(x = gene_sets$entrez_gene, f = gene_sets$gs_name)
GSVA::gsva(expr, gset.idx.list = msigdbr_list, method = "ssgsea")

