## Write covariates for EPIC data
# Libraries
library(tidyverse)
library(Matrix)

# Data ---------------------------------------------------------------

samplesheet <- readRDS("../../data/raw_data/samplesheet_UK_2017.06_383samples_Ireland_2017.03_298samples_cleaned.rds")
missingness_matrix <- readRDS("../../data/raw_data/beta_dasen_NAs_UK_2017.06_383samples_Ireland_2017.03_298samples_missingness_matrix.rds")

# Covariates ---------------------------------------------------------------

qcovariates <- c("Predicted_Age", "smokingScore")
covariates <- c("Sex", "Batch")
phenotype <- c("C9ORF72.Status")
levels <- c("normal", "expanded")

samplesheet <- samplesheet %>%
  dplyr::filter((!!as.name(phenotype) %in% levels))

samplesheet <- samplesheet[complete.cases(samplesheet[,c(covariates, qcovariates, phenotype)]),]

covar <- apply(samplesheet[,covariates,drop=FALSE], 2, FUN = function(x) as.numeric(as.factor(x)))
covar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, covar)

qcovar <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name, samplesheet[,qcovariates])

phenotype <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name,
             ifelse(samplesheet[[phenotype]] == levels[2], 1, 0))

# Check missingness ---------------------------------------------------------------
missingness <- missingness_matrix[,samplesheet$Sample_Name]
props <- Matrix::rowMeans(missingness)
props_0.05 <- props[props > 0.05]

## Save
write_tsv(qcovar, path = "../../data/temp/qcovar_predictedage_smoking_EPIC.txt", col_names=FALSE)
write_tsv(covar, path = "../../data/temp/covar_sex_batch_EPIC.txt", col_names=FALSE)
# C9 status
write_tsv(phenotype, path = "../../data/temp/pheno_c9_EPIC.txt", col_names=FALSE)
write_lines(names(props_0.05), path = "../../data/temp/probes_miss0.05_EPIC.txt")
