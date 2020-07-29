## Libraries ---------------------------
library(minfi)
library(wateRmelon)
library(readr)
library(magrittr)
library(optparse)
library(dplyr)

## Arguments ---------------------------
remove_samples <- read_lines("../../data/raw_data/sample_qc_outliers_NL.txt")

## Load data ---------------------------
mset <- readRDS("../../data/raw_data/mset_NL.rds")

## Remove failed samples ---------------------------

ncol1 <- ncol(mset)
mset <- mset[,!(colnames(mset) %in% remove_samples)]
ncol2 <- ncol(mset)
cat(sprintf("\nRemoved %s samples from the mset!", ncol1-ncol2))

## Get beta-values  ---------------------------
cat("\nRunning dasen function..")
mset <- dasen(mset)
M <- mset@assayData$methylated
U <- mset@assayData$unmethylated
rm(mset)
Total <- M + U
cat("\nDone!")

## Save betas  ---------------------------
cat("\nSaving the total intensities!")
saveRDS(M, file = "../../data/processed/M_intensities_NL.rds")
saveRDS(U, file = "../../data/processed/U_intensities_NL.rds")
saveRDS(Total, file = "../../data/processed/total_intensities_NL.rds")

## SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
