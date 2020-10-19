## ISVA
library(tidyverse)
library(sva)
library(meffil)

## Data
c9_status <- readr::read_tsv("../../data/temp/pheno_c9.txt", col_names=FALSE)
covar <- readr::read_tsv("../../data/temp/qcovar_predictedage_smoking.txt", col_names=FALSE)
qcovar <- readr::read_tsv("../../data/temp/covar_sex_batch.txt", col_names=FALSE)
rmprobes <- readr::read_lines("../../data/temp/probes_miss0.05.txt")
anno <- readr::read_tsv("../../data/extdata/annotated_450k.opi", col_names = FALSE)
anno <- anno %>% filter(X1 %in% c("X", "Y"))
beta <- readRDS("../../data/raw_data/beta_cleaned_NL.rds")

## Subset
beta <- beta[!rownames(beta) %in% c(rmprobes, anno$X2), c9_status$X1]
beta <- impute::impute.knn(beta)
beta <- beta$data

# Perform SVA
covars <- covar %>% dplyr::left_join(qcovar[,c(1,3:ncol(covar))], by = c("X1"))

ewas.ret <- meffil.ewas(beta, variable=c9_status$X3,
  covariates=covars[,3:ncol(covars)], winsorize.pc=NA)

saveRDS(ewas.ret, file = "../../data/output/ewas/ewas_meffil_c9_missing0.05.rds")
