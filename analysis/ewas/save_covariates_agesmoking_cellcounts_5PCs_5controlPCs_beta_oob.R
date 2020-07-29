## Write covariates to use in OSCA

# Libraries ---------------------------------------------------------------
library(tidyverse)

# Samplesheet ---------------------------------------------------------------
samplesheet <- readRDS("../../data/raw_data/samplesheet_NL_cleaned.rds")
pca <- readr::read_delim("../../data/processed/pca_beta_NL_oob_normalized.eigenvec", col_names = FALSE, delim = " ")
samples <- pca$X1
pca$X1 <- pca$X2 <- NULL
colnames(pca) <- paste0("PC", 1:ncol(pca))
pca <- pca %>%
  dplyr::mutate(Sample_Name = samples) %>% dplyr::select(Sample_Name, dplyr::everything())

# Covariates ---------------------------------------------------------------

qcovariates <- c("Predicted_Age", "smokingScore", "NK", "CD4T", "CD8T", "Gran", "Mono",
"PC1_control", "PC2_control", "PC3_control", "PC4_control", "PC5_control")
covariates <- c("Sex", "Batch")
phenotype <- c("Inferred_C9_status")
levels <- c("wt", "long")

# Select ---------------------------------------------------------------

# Case-only
samplesheet <- samplesheet %>% dplyr::filter(Pheno == 2)

# Select covariates/phenotype
samplesheet <- samplesheet[complete.cases(samplesheet[,c(covariates, qcovariates, phenotype)]),]
samplesheet <- samplesheet %>%
  dplyr::filter((!!as.name(phenotype) %in% levels))

covar <- apply(samplesheet[,covariates,drop=FALSE], 2, FUN = function(x) as.numeric(as.factor(x)))
covar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, covar)

qcovar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, samplesheet[,qcovariates])
qcovar <- qcovar %>% dplyr::left_join(pca[,c("Sample_Name", "PC1", "PC2", "PC3", "PC4", "PC5")], by = c("IID" = "Sample_Name"))

phenotype <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name,
                        ifelse(samplesheet[[phenotype]] == levels[2], 1, 0))

## Save
write_tsv(qcovar, path = "../../data/temp/qcovar_predictedage_smoking_cellcounts_5controlPCs_5PCs_beta_oob.txt", col_names=FALSE)
