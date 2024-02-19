### STEP 1: Generating files for differential expression analysis ###

## This script will prepare a DESeq2 object for downstream analysis using the
## RNA-seq counts table prepared from the CITI Group pipeline and sample
## metadata from obtained from cBioPortal (CCF internal version).

# Load libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  })

# Import files with counts information and clinical data ------------------
dir <- getwd()
data_path <- file.path(dir,"data","raw_count_matrix_clean_BrainTumorCenter_2004449.csv")

## Un-comment and run the following file transfer command while in terminal, ensure
## that server connection is intact on desktop. Use personal file paths as required.
# rsync -r --progress --size-only \
# /Volumes/MicrobialCompositionandAnalyticsCore/CITI_DATA_LATHIA/NICO_STUDY/BrainTumorCenter_2004449/RNASeq/data/Deliverables/CountsAndMatrices/raw_count_matrix_clean_BrainTumorCenter_2004449.csv \
# /Users/dashorh/rna-seq/data

data <- readr::read_csv(file = data_path)
meta_path <- file.path(dir,"meta", "ccf_bci_gbm_lathia_clinical_data.tsv")
meta <- readr::read_tsv(file = meta_path)

metaTable <- meta %>%
  # select releveant columns from the clinical data
  dplyr::select(`Patient ID`,
                `Sample ID`,
                `final diagnosis`,
                `Group Description`,
                `Notes`,
                `record_id`,
                `Sample Name`) %>%
  # filter for samples with histologically confirmed diagnosis of GBM
  dplyr::filter(`final diagnosis` == "GBM") %>%
  # filter for rows with non-recurrent surgical resections (keep rows with NA)
  dplyr::filter((`Notes` != "RECURRENT TUMOR") %>% replace_na(TRUE)) %>%
  # re-factor the description to make labels consistent and remove ambiguity from hybrid labels
  dplyr::mutate(fct_collapse(.f = factor(`Group Description`),
                             core = c("central tumor", "core", "Core", "necrotic core", "enhancing rim, core"),
                             enhancing = c("enhancing rim", "non-enhancing rim, enhancing rim"),
                             non_enhancing = c("non-enhancing rim","non-enhancing tumor", "core, non-enhancing rim"))) %>%
  # rename and rearrange columns for better readability
  dplyr::rename(patient = `Patient ID`,
                sample_id = `Sample ID`,
                sample_name = `Sample Name`,
                record = `record_id`,
                diagnosis = `final diagnosis`,
                description = `fct_collapse(...)`) %>%
  dplyr::select(sample_name, patient, sample_id, record, diagnosis, description) %>%
  # set row names to match the column names in raw counts data
  tibble::column_to_rownames(var = "sample_name") %>%
  # convert charcter columns to factors
  dplyr::mutate(across(where(is.character), as.factor))

# Utilize prepared matrix of read counts to create DESeqDataSet
raw_counts <- data %>% column_to_rownames(var = "...1") %>%
  # remove gene annotation columns
  dplyr::select(-EnsemblID, -GeneSymbol, -EntrezID) %>%
  # remove patients from further analysis according to metadata
  dplyr::select(all_of(rownames(metaTable)))

all(rownames(metaTable) %in% colnames(raw_counts))
all(colnames(raw_counts) == rownames(metaTable)) # columns of count matrix and rows of the column data are in the same order

# Create design contrasting 3 tumor regions
dds_regions <- DESeqDataSetFromMatrix(countData = raw_counts,
                                      colData = metaTable,
                                      design = ~ patient + description) 

# Additional feature data may be added to the metadata columns
featureData <- dplyr::select(data, c(last_col(2):last_col(0)))
mcols(dds_regions) <- DataFrame(mcols(dds_regions), featureData)

# Pre-filter low count genes, require raw count of at least 10 in a minimal number of samples
smallestGroupSize <- min(table(metaTable$description))
keep <- rowSums(counts(dds_regions)) >= smallestGroupSize
dds_regions <- dds_regions[keep, ]

# Specify the reference level
dds_regions$description <- relevel(dds_regions$description, ref = "non_enhancing")

# Save DESeq2 object in .RData file
save(dds_regions, file = "data/NICO_DESeqObject.RData")
