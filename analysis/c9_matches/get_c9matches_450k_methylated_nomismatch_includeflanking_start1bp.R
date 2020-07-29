# Identify 450k probes that map to the C9 hexanucleotide (GGCCCC) repeat using the
# probe_mapper_nonref.R script

# Libraries ---------------------------------------------------------------------
library(tidyverse)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
data(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# Functions ---------------------------------------------------------------------
source("../../R/probe_mapper_nonref.R")

# Function to annotate the data
annotate <- function(matches, anno) {
  matches <- matches %>%
    dplyr::left_join(anno[,c("Name", "Type", "ProbeSeqA", "ProbeSeqB", "Color")],
                     by = c("Probe" = "Name")) %>%
    dplyr::mutate(predicted_color = ifelse(sbe_site %in% c("A", "T"), "Red",
                                           ifelse(sbe_site %in% c("C", "G"), "Grn", NA)))
}

# Analysis ---------------------------------------------------------------------

# In silico repeat
repeat_sequence <- paste(rep("GGCCCC", 10), collapse="")
flanking1 <- "ACCGCAGCCCCGCCCCGGGCCCGCCCCCGGGCCCGCCCCGACCACGCCCC"
flanking2 <- "TAGCGCGCGACTCCTGAGTTCCAGAGCTTGCTACAGGCTGCGGTTGTTTC"
repeat_sequence <- paste0(flanking1, repeat_sequence, flanking2)
next_base = "C"
prev_base = "A"

# Load probe sequences
probe_sequences <- readr::read_tsv("../../data/processed/probe_sequences_450k.txt")

# Annotation
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

### Run probe_mapper, assuming the repeat is methylated, no mismatches
matches_methylated_nomismatch <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                              probe_sequences = probe_sequences,
                                              max_width = 35,
                                              min_width = 1, stepsize = 1,
                                              allow_mismatch = FALSE , allow_indel = FALSE, use_Y = FALSE,
                                              methylation_status = "methylated", nr_cores = 1)

saveRDS(matches_methylated_nomismatch, file = "../../data/processed/c9_matches/450k_matches_methylated_nomismatch_includeflanking_start1bp.rds")


## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
