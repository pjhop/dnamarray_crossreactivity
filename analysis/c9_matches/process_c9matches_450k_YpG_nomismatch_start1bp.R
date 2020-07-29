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
repeat_sequence <- paste(rep("GGCCCCGGCCCCGGCCCC", 10), collapse="")
next_base = "G"
prev_base = "C"

# Load probe sequences
probe_sequences <- readr::read_tsv("../../data/processed/probe_sequences_450k.txt")

# Annotation
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))

### Run probe_mapper,
matches_YpG_nomismatch <- readRDS("../../data/processed/c9_matches/450k_matches_YpG_nomismatch_start1bp.rds")

cat("Cleaning matches\n")
matches_YpG_nomismatch <- clean_matches(matches_YpG_nomismatch, anno = anno, probe_sequences = probe_sequences)
cat("Selecting maximum width for each probe\n")
matches_YpG_nomismatch_max <- select_max_width(matches_YpG_nomismatch)
cat("Selecting unique sequence for each match\n")
matches_YpG_nomismatch_final <- select_unique_sequence(matches_YpG_nomismatch_max)
cat("Annotation\n")
matches_YpG_nomismatch_final <- annotate(matches_YpG_nomismatch_final, anno = anno)
cat("Saving\n")
saveRDS(matches_YpG_nomismatch_final, file = "../../data/processed/c9_matches/450k_matches_YpG_nomismatch_final_start1bp.rds")

## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
