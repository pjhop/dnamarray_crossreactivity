## Libraries  ---------------------------
library(minfi)
library(tidyverse)
library(matrixStats)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(DNAmCrosshyb)

## Data  ---------------------------
rgset <- readRDS("../../data/raw_data/rgset_UK_2017.06_383samples_Ireland_2017.03_298samples.rds")

## Get OOB betas/intensities after NATEN normalization
normalized_oob <- get_OOB(rgset, keep = "both", normalized = TRUE)

## Save  ---------------------------
saveRDS(normalized_oob$beta, file = "../../data/processed/beta_UK_2017.06_383samples_Ireland_2017.03_298samples_oob_normalized.rds")
saveRDS(normalized_oob$total, file = "../../data/processed/total_intensities_oob_normalized_UK_2017.06_383samples_Ireland_2017.03_298samples.rds")

## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
