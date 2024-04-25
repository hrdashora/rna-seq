## STEP 6: Generate files for CIBERSORT deconvolution ##

# Visualize data tables from GLIOVIS output as well

# Load libraries
library(tidyverse)
library(ggpubr)
library(pheatmap)
library(reshape2)

# Import file that contains TPM information
dir <- getwd()
tpm_path <- file.path(dir,"tpm_matrix_BrainTumorCenter_2004449.csv")
tpm <- readr::read_csv(file = tpm_path)
# Name matrix by Gene Symbol and select duplicated gene symbol with highest rowMeans
samples <- colnames(tpm)[2:58]
tpm <- tpm %>% mutate(mean_tpm = rowMeans(select(., all_of(samples)))) %>%
  arrange(mean_tpm) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  filter(!is.na(GeneSymbol)) %>%
  arrange(GeneSymbol) %>%
  dplyr::select(-...1, -mean_tpm, -EnsemblID, -EntrezID) %>%
  relocate(GeneSymbol)
  
# Write data frame to tab-delimited file for CIBERSORT
write.table(tpm, file = "tpm_matrix_CIBERSORT.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Write data frame to CSV for GLIOVIS
tpm_m <- as.matrix(data.frame(tpm[-1], row.names = tpm$GeneSymbol))
tpm_t <- t(tpm_m) %>% as_tibble(rownames = "Sample")
write_csv(tpm_t, file = "tpm_matrix_GLIOVIS.csv") # complete file modifications in excel

# Import GLIOVIS data and merge with sampleTable from STEP 1
# Ensure that patients 7 and 18 are removed
metatable <- coldata %>% rownames_to_column(var = "name") %>% select(-group1, -group2)
decon_newman <- read_csv("results/deconvoluteme_gliovis_newman.csv") %>%
  left_join(x = ., y = metatable, by = c("Sample" = "name")) %>%
  group_by(description) %>%
  na.omit()
decon_bindea <- read_csv("results/deconvoluteme_gliovis_bindea.csv") %>%
  left_join(x = ., y = metatable, by = c("Sample" = "name")) %>%
  group_by(description) %>%
  na.omit()
sub <- read_csv("results/subtypeme_gliovis.csv") %>%
  left_join(x = ., y = metatable, by = c("Sample" = "name")) %>%
  group_by(description, majority.call) %>%
  na.omit()
est <- read_csv("results/estimateme_gliovis.csv") %>%
  left_join(x = ., y = metatable, by = c("Sample" = "name")) %>%
  group_by(description) %>%
  na.omit()

sub_count <- sub %>% summarise(type_count = n())
est_count <- pivot_longer(est, cols = c(ImmuneScore, StromalScore, ESTIMATEScore), names_to = "score")

ggp <- ggplot(sub_count,
              aes(x = description,
                  y = type_count,
                  fill = majority.call)) +
  geom_bar(position = "fill", stat = "identity") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "GBM Subtype Classification",
       y = "Proportion of Samples (%)",
       fill = "Subtype") +
  scale_x_discrete(name = NULL,
                   labels = c("Core","Enhancing","Non-enhancing")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggp

ggp <- ggplot(est_count,
              aes(x = score,
                  y = value,
                  fill = description)) + 
  geom_boxplot() + 
  labs(title = "Infiltration of Immune Cells in Tumor Tissue",
       subtitle = "via ESTIMATE Algorithm",
       y = "Score",
       fill = "Location") +
  scale_fill_discrete(labels = c("Core","Enhancing","Non-enhancing")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())
ggp

ggpubr::compare_means(ESTIMATEScore ~ description, data = est)
ggpubr::compare_means(ImmuneScore ~ description, data = est)
ggpubr::compare_means(StromalScore ~ description, data = est)

### YOU ARE HERE ###

decon_num <- as.matrix(decon[,2:28])
rownames(decon_num) <- decon$Sample
decon_num_scale = scale(decon_num)

subj_df <- data.frame("Subject" = decon$subject)
rownames(subj_df) <- rownames(decon_num) # name matching

cat_df = data.frame("Location" = decon$description)
rownames(cat_df) = rownames(decon_num)
cat_order_df <- as.data.frame(cat_df[order(cat_df$Location), ])
rownames(cat_order_df) = rownames(decon_num)[order(cat_df$Location)]
colnames(cat_order_df) <- c("Location")

sub_samp_ordered <- decon_num_scale[rownames(cat_order_df), ]

pheatmap(sub_samp_ordered,
         cluster_cols = T,
         cluster_rows = F,
         annotation_row = cat_order_df,
         labels_row = "",
         main = "pheatmap default")

compare_means(Activated.B.cell ~ description, data = decon)

ggplot(decon,
       aes(x = description,
           y = ImmuneScore,
           fill = description)) + 
  geom_boxplot() + 
  labs(title = "Infiltration of Immune Cells in Tumor Tissue",
       subtitle = "via ESTIMATE Algorithm",
       caption = "Data source: NICO-Myriad",
       y = "Immune Score",
       fill = "Location") +
  scale_x_discrete(name = "Location",
                   labels = c("Core","Enhancing","Non-enhancing")) +
  scale_fill_discrete(labels = c("Core","Enhancing","Non-enhancing")) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

x <- decon %>% dplyr::select(-region, -subject, -number, -group2, -description) %>% melt(.)
facets <- c("Activated.B.cell", "Activated.CD4.T.cell", 
            "Regulatory.T.cell", "CD56bright.natural.killer.cell",
            "Natural.killer.cell")

labs <- c("Activated B Cell", 
          "Activated CD4 T cell", 
          "Regulatory T cell",
          "CD56-bright Natural Killer Cell",
          "Natural Killer Cell")
names(labs) <- facets

dec <- ggboxplot(x[x$variable %in% facets, ],
          x = "group1",
          y = "value",
          color = "group1",
          palette = "npg",
          facet.by = "variable") +
  labs(x = NULL, y = "Enrichment Score") +
  facet_wrap(~variable, labeller = labeller(variable = labs), scales = "free_y", ncol = 5) +
  scale_x_discrete(labels = c("interior" = "Interior", 
                              "exterior" = "Exterior")) +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_compare_means(method = "wilcox.test",
                     ref.group = "interior",
                     label = "p.signif",
                     label.y.npc = "top")

ggsave(filename = "decon.png", plot = dec, device = png,
       width = 30, height = 10, units = "cm", dpi = "retina")

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(ggplot2)
library(ggpubr)
library(pheatmap)
library(RColorBrewer)
clin <- read_delim("results/clinical_data.txt")
ciber <- read_csv("results/CIBERSORTx_Job12_Results.csv") %>%
  select(-`Absolute score (sig.score)`, -Correlation, -`P-value`, -RMSE) %>%
  left_join(clin, by = c("Mixture" = "Tumor_Sample_Barcode"))
ciber_long <- tidyr::pivot_longer(ciber, !c(Mixture,Patient_ID,Group_Description))

ggplot(ciber_long, aes(Mixture, name, fill = value)) +
  geom_tile() +
  theme(axis.text.x = element_text(angle = 90))


df <- clin %>% arrange(Group_Description) %>%
  column_to_rownames(var = "Tumor_Sample_Barcode") %>%
  select(Group_Description)

mat <- as.matrix(ciber[,2:23])
rownames(mat) <- ciber$Mixture
mat_order <- mat[rownames(df),]
pheat <- pheatmap(mat = mat_order,
                  annotation_row = df,
                  border_color = F,
                  cluster_rows = F,
                  cluster_cols = F,
                  scale = "column"
                  )
pheat

coul <- brewer.pal(12,"Set3") 
coul <- colorRampPalette(coul)(22)

ggplot(ciber_long %>% filter(name == c("Macrophages M0")),
       aes(Group_Description, value, fill = name)) +
  geom_boxplot() +
  scale_fill_manual(name = "Immune Cell Types", values = coul) +
  scale_x_discrete(labels = c("Core","CE","NCE"))+
  #scale_y_continuous(trans = "log10") +
  labs(x = "Glioblastoma Regions", y = "Proportion") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90), legend.key.size = unit(5,'mm')) +
  guides(fill = guide_legend(ncol = 1))

ggplot(ciber_long, aes(Patient_ID, value, fill = name)) +
  geom_col(position = "fill") +
  scale_fill_manual(name = "Immune Cell Types", values = coul) +
  labs(x = "Patient", y = "Relative Proportion") +
  theme_classic() +
  theme(axis.text.x = element_blank(), legend.key.size = unit(5,'mm')) +
  guides(fill = guide_legend(ncol = 1))

compare_means(value ~ Group_Description, data = ciber_long, group.by = "name") %>%
  dplyr::arrange(p)

ciber_data <- dplyr::left_join(ciber_long, clin, by = c("Mixture" = "Tumor_Sample_Barcode"))
p <- ggplot(ciber_data, aes(x=Group_Description, y=value)) +
  facet_wrap(~name, scales = "free") +
  geom_boxplot(alpha=0.8) +
  geom_point() +
  stat_compare_means(size = 2, vjust = .5) +
  theme(axis.text.x = element_blank())
p
