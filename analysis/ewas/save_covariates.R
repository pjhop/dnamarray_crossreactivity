## Write covariates to use in OSCA

# Libraries ---------------------------------------------------------------
library(tidyverse)
library(Matrix)

# Samplesheet ---------------------------------------------------------------
samplesheet <- readRDS("../../data/raw_data/samplesheet_NL_cleaned.rds")
missingness_matrix <- readRDS("../../data/raw_data/beta_cleaned_dasen_NAs_missingness_matrix.rds")

# Covariates ---------------------------------------------------------------

qcovariates <- c("Predicted_Age", "smokingScore")
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

phenotype <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name,
                        ifelse(samplesheet[[phenotype]] == levels[2], 1, 0))

keep <- data.frame(IID = samplesheet$Sample_Name, FID = samplesheet$Sample_Name)

# Check missingness ---------------------------------------------------------------
missingness <- missingness_matrix[,samplesheet$Sample_Name]
props <- Matrix::rowMeans(missingness)
props_0.05 <- props[props > 0.05]

## Save
write_tsv(qcovar, path = "../../data/temp/qcovar_predictedage_smoking.txt", col_names=FALSE)
write_tsv(covar, path = "../../data/temp/covar_sex_batch.txt", col_names=FALSE)
write_tsv(keep, path = "../../data/temp/keep.txt", col_names=FALSE)

# Only sex
covar <- covar %>% dplyr::select(-Batch)
write_tsv(covar, path = "../../data/temp/covar_sex.txt", col_names=FALSE)#

# C9 status
write_tsv(phenotype, path = "../../data/temp/pheno_c9.txt", col_names=FALSE)

# Probes with >0.05 missingness
write_lines(names(props_0.05), path = "../../data/temp/probes_miss0.05.txt")
