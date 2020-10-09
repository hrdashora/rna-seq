load("rnaseq_objects.Rdata")

# Create design matrix and contrasts
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)
contr.matrix <- makeContrasts(
  H.DKOVsCtrl = h.dko - h.ctrl,
  L.DKOVsCtrl = l.dko - l.ctrl,
  levels = design)

# Voom transform the data
par(mfrow=c(1,1))
v <- voom(count_dge,design,plot = TRUE)
## Compare boxplots of normalised data
par(mfrow=c(1,2))
boxplot(logcounts, xlab="", ylab="Log2 counts per million",las=2,
        main="Unnormalised logCPM")
abline(h=median(logcounts),col="blue") # blue horizontal line that corresponds to the median logCPM
boxplot(v$E, xlab="", ylab="Log2 counts per million",las=2,
        main="Voom transformed logCPM")
abline(h=median(v$E),col="blue")
vfit <- limma::lmFit(v, design)
vfit <- contrasts.fit(vfit, contr.matrix)
efit <- eBayes(vfit)
res_efit <- decideTests(efit)
#plotSA(efit)
#summary(res_efit)
#topTable(efit,coef="H.CtrlVsDKO",sort.by="p")

# Testing relative to a threshold
## Interested in genes that have a absolute logFC of 0.585 (FC = 1.5)
fit_treat <- treat(efit,lfc=0.585)
res_treat <- decideTests(fit_treat)
#summary(res_treat)
par(mfrow=c(1,1))
plotMD(fit_treat,coef=1,status=res_treat[,"H.DKOVsCtrl"], values=c(-1,1))
abline(h=0,col="grey")
# glMDPlot(fit_treat, coef=1, counts=count_dge$counts, groups=group,
#          status=res_treat, side.main="external_gene_name", main="H.CtrlVsDKO",
#          folder="md")

# Create gene annotation table
ensLookup <- rownames(efit)
#ensLookup <- gsub("\\.[0-9]*$", "", ens)
mart <- useMart("ENSEMBL_MART_ENSEMBL")
mart <- useDataset("mmusculus_gene_ensembl", mart)
annotLookup <- getBM(
  mart=mart,
  attributes=c("ensembl_gene_id",
               "ensembl_gene_id_version",
               "entrezgene_id",
               "mgi_symbol",
               "go_id",
               "entrezgene_description",
               "transcript_length"),
  filter="ensembl_gene_id_version",
  values=ensLookup,
  uniqueRows=TRUE)
annotLookup <- annotLookup[!duplicated(annotLookup$ensembl_gene_id_version),]
annotLookup <- data.frame(ensLookup[match(annotLookup$ensembl_gene_id_version,
                                          ensLookup)],
                          annotLookup)
colnames(annotLookup) <- c("original_id",
                           c("ensembl_gene_id",
                             "ensembl_gene_id_version",
                             "entrezgene_id",
                             "mgi_symbol",
                             "go_id",
                             "entrezgene_description",
                             "transcript_length"))
ord_annotLookup <- annotLookup[match(rownames(efit),annotLookup$original_id),]
efit$genes <- ord_annotLookup[,-1] # drop duplicated version column
#fit_treat$genes <- ord_annotLookup[,-1] 

h_limma_res <- topTable(efit,coef="H.DKOVsCtrl",sort.by="p",n="Inf")
write.csv(h_limma_res,file="Heart_DKOVsCtrl.csv",row.names=FALSE)
l_limma_res <- topTable(efit,coef="L.DKOVsCtrl",sort.by="p",n="Inf")
write.csv(l_limma_res,file="Lung_DKOVsCtrl.csv",row.names=FALSE)

# Create Venn Diagram
res_common <- which(res_treat[,1]!=0 & res_treat[,2]!=0)
#head(fit_treat$genes$mgi_symbol[res_common], n=20)
vennDiagram(res_treat[,1:2], circle.col=c("turquoise", "salmon"))
write.fit(fit_treat, res_treat, file="treatresults.txt")

# Volcano plot, highlight top 100 most DE genes
volcanoplot(efit,coef=1,highlight=100,
            names=efit$genes$mgi_symbol)

# Gene ontology testing
## Load in the mouse c2 gene sets
Mm.c2 <- readRDS("~/R/rna-seq/data/Mm.c2.all.v7.1.entrez.rds")
## Recreate voom objects based on Enterez IDs
cases <- efit$genes[!is.na(efit$genes$entrezgene_id),]
count_keep_go <- count_keep[match(cases$ensembl_gene_id_version,rownames(count_keep)),]
rownames(count_keep_go) <- cases$entrezgene_id
count_dge_go <- edgeR::DGEList(counts = count_keep_go)
count_dge_go <- calcNormFactors(count_dge_go)
v_go <- voom(count_dge_go,design,plot = TRUE)
vfit_go <- limma::lmFit(v_go, design)
vfit_go <- contrasts.fit(vfit_go, contr.matrix)
efit_go <- eBayes(vfit_go)
## Testing with goana
h_go <- goana(efit_go, coef="H.DKOVsCtrl",species = "Mm")
l_go <- goana(efit_go, coef="L.DKOVsCtrl",species = "Mm")
topGO(h_go, n=10)
topGO(l_go, n=10)
## Testing with camera
idx <- ids2indices(Mm.c2,id=rownames(v_go))
cam.H.DKOVsCtrl <- camera(v_go,idx,design,
                          contrast=contr.matrix[,1],
                          inter.gene.cor=0.05) 
cam.L.DKOVsCtrl <- camera(v_go,idx,design,
                          contrast=contr.matrix[,2],
                          inter.gene.cor=0.05)
write.csv(cam.L.DKOVsCtrl,file="cam.H.DKOVsCtrl.csv")
write.csv(cam.L.DKOVsCtrl,file="cam.L.DKOVsCtrl.csv")
head(cam.H.DKOVsCtrl,5)
head(cam.L.DKOVsCtrl,5)

barcodeplot(efit_go$t[,1],
            index=idx$EGUCHI_CELL_CYCLE_RB1_TARGETS,
            main="T-statistic: EGUCHI_CELL_CYCLE_RB1_TARGETS")
# Genes in this gene set tend to be up-regulated between
# DKO and Control in heart endothelial cells

# 4-way comparison plot
heart_res_df <- readr::read_csv(file.path(getwd(),"Heart_DKOVsCtrl.csv"))
lung_res_df <- readr::read_csv(file.path(getwd(),"Lung_DKOVsCtrl.csv"))
fourway <- tibble(ensembl_id = heart_res_df$ensembl_gene_id_version,
                  entrez_id = heart_res_df$entrezgene_id,
                  symbol = heart_res_df$mgi_symbol,
                  h.logFC = heart_res_df$logFC,
                  h.adj.P.val = heart_res_df$adj.P.Val)
lung_merge <- tibble(ensembl_id = lung_res_df$ensembl_gene_id_version,
                     l.logFC = lung_res_df$logFC,
                     l.adj.P.val = lung_res_df$adj.P.Val)
fourway <- inner_join(fourway,lung_merge, by = "ensembl_id")
fourway <- fourway[(fourway$h.adj.P.val < 0.05 & fourway$l.adj.P.val < 0.05),]

fourway <- fourway %>% mutate(status = case_when(
  abs(h.logFC) > 0.585 & abs(l.logFC) > 0.585 ~ "Both",
  abs(h.logFC) > 0.585 & abs(l.logFC) <= 0.585 ~ "Heart",
  abs(h.logFC) <= 0.585 & abs(l.logFC) > 0.585 ~ "Lung",
  TRUE ~ "Neither"
))
ggplot2::ggplot(data = fourway, aes(x = h.logFC, y = l.logFC)) +
  geom_point(aes(color = status), size = 1.25, alpha = 0.25) +
  scale_color_manual(values = c("Both" = "purple",
                                "Heart" = "red",
                                "Lung"="blue",
                                "Neither"="grey"
                                ),
                     labels = c("Both",
                                "Heart Only",
                                "Lung Only",
                                "Neither"))+
  geom_hline(aes(yintercept = -0.585), linetype="dashed") +
  # geom_text(aes(x = 7,y = -0.6,label = "-1.5 FC", vjust = 1),
  #           size = 3.25,
  #           colour = "black") +
  geom_hline(aes(yintercept = 0.585), linetype="dashed") +
  # geom_text(aes(x = 7, y = 0.6,label = "1.5 FC", vjust = 0),
  #           size = 3.25,
  #           colour = "black") +
  geom_vline(aes(xintercept = -0.585), linetype="dashed") +
  geom_vline(aes(xintercept = 0.585), linetype="dashed") +
  labs(x = "Heart (Log2 Fold Change)" , y = "Lung (Log2 Fold Change)",
       title = "DIfferential Expression of Significant Genes in M.m. Endothelial Cells",
       subtitle = "K2K4 DKO - CRE Ctrl, Adj. p-value < 0.05, N = 6035",
       color = "FC > 1.5")
  
# Calculate RPKM
myRPKM <- edgeR::rpkm(y = count_dge,
                      gene.length = efit$genes$transcript_length)
write.csv(myRPKM,file="RPKM.csv")

# Principal Component Analysis
mds <- plotMDS(count_dge)

select <- order(rowMeans(v$E), decreasing = TRUE)[1:500]
highexprgenes_counts <- v$E[select,]
data_for_PCA <- t(highexprgenes_counts)
cmd <- cmdscale(dist(data_for_PCA), k=3, eig=TRUE)
eig_pc <- cmd$eig * 100 / sum(cmd$eig)
barplot(eig_pc,
        las=1,
        xlab="Dimensions", 
        ylab="Proportion of explained variance (%)", y.axis=NULL,
        col="darkgrey")
plot(cmd$points[,1], -cmd$points[,2], type="n", xlab="Dimension 1", ylab="Dimension 2", main="")
text(cmd$points[,1], -cmd$points[,2], rownames(cmd$points), cex=0.8) 

project.pca <- prcomp(data_for_PCA)
#summary(project.pca)
project.pca.proportionvariances <- ((project.pca$sdev^2) / (sum(project.pca$sdev^2)))*100
barplot(project.pca.proportionvariances, cex.names=1, xlab=paste("Principal component (PC), 1-", length(project.pca$sdev)), ylab="Proportion of variation (%)", main="Scree plot", ylim=c(0,100))
plot(project.pca$x, type="n", main="Principal components analysis bi-plot", xlab=paste("PC1, ", round(project.pca.proportionvariances[1], 2), "%"), ylab=paste("PC2, ", round(project.pca.proportionvariances[2], 2), "%"))
points(project.pca$x, col="black", pch=16, cex=1)

