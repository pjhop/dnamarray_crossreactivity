## Perform a meta-analysis on the test-statistics of the three batches

## Libraries ---------------------------

library(tidyverse)
library(magrittr)
library(meta)
library(BiocParallel)
library(optparse)

## Data ---------------------------

# Load input: Batches, Nr_cases, Nr_controls, filepath
filepaths <- tibble(Batch = c("NL_2015.10_1000samples", "Netherlands_ReRun_193samples", "Netherlands_ReRun_2016.06_1736samples"),
                    path = c("../../data/output/ewas/mlma_loco_c9_covarsex_qcovaragesmoking_missing0.05_NL_2015.10_1000samples.linear",
                             "../../data/output/ewas/mlma_loco_c9_covarsex_qcovaragesmoking_missing0.05_Netherlands_ReRun_193samples.linear",
                             "../../data/output/ewas/mlma_loco_c9_covarsex_qcovaragesmoking_missing0.05_Netherlands_ReRun_2016.06_1736samples.loco.mlma"))
samplesheet <- readRDS("../../data/raw_data/samplesheet_NL_cleaned.rds")
overview <- samplesheet %>%
    dplyr::filter(Pheno == 2) %>%
    dplyr::group_by(Batch) %>%
    dplyr::count(Inferred_C9_status) %>%
    dplyr::filter(!is.na(Inferred_C9_status)) %>%
    tidyr::spread(key = Inferred_C9_status,
                  value = n)
overview <- overview %>% dplyr::left_join(filepaths, by = "Batch")
colnames(overview) <- c("Name", "Nr_cases", "Nr_controls", "filepath")

## Functions ---------------------------

# Load test statistics
load_test_statistics <- function(file) {
	stats <- readr::read_tsv(file)
	stats
}

get_probes <- function(file) {
	probes <- as.data.frame(file)[,"Probe"]
	probes
}

get_stats <- function(batch, stats) {
	# For each stats
  stats_ <- stats[[batch]]
	stats_ <- tibble::tibble(Probe = stats_$Probe, beta = stats_$b, se = stats_$se, Batch = batch)
	stats_
}

meta_analysis <- function(probe, stats) {
  stats_ <- stats %>% dplyr::filter(Probe == probe)
  b <- stats_$beta
  se <- stats_$se
  n.e <- stats_$Nr_cases
  n.c <- stats_$Nr_controls
  met <- meta::metagen(TE = b, seTE = se, n.e = n.e, n.c = n.c)
  tib <- tibble::tibble(Probe = probe, b = met$TE.fixed, se = met$seTE.fixed, p = met$pval.fixed,
                Q = met$Q, Q.p = pchisq(Q, df = met$df.Q, lower.tail = FALSE))
  tib
}

## Run analyses ---------------------------

# Load test-statistics
all_stats <- purrr::map(overview$filepath, .f = load_test_statistics)
names(all_stats) <- overview$Name
all_probes <- purrr::map(all_stats, get_probes)
probes_intersect <- Reduce(intersect, all_probes)

## Combine all betas and SEs
stats <- purrr::map_df(names(all_stats), .f = get_stats, stats = all_stats)

## Add info
stats <- stats %>%
  dplyr::left_join(overview[,c("Name", "Nr_cases", "Nr_controls")],
                    by = c("Batch" = "Name"))

## Run meta-analyses in parallel
cat("\nRunning meta-analysis...")
meta_pvals <- BiocParallel::bplapply(unique(stats$Probe), FUN = meta_analysis, stats = stats,
                                     BPPARAM = BiocParallel::MulticoreParam(workers = 10))
meta_pvals <- do.call(rbind, meta_pvals)


## Save data
readr::write_tsv(meta_pvals, path = "../../data/output/ewas/meta_analysis.tsv")


## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
