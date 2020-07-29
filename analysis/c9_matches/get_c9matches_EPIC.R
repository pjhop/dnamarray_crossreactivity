# Identify probes that map to the C9 hexanucleotide (GGCCCC) repeat using the
# probe_mapper_nonref.R script

# Libraries ---------------------------------------------------------------------
library(tidyverse)
library(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
data(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)

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
next_base = "G"
prev_base = "C"

# Load probe sequences
probe_sequences <- readr::read_tsv("../../data/processed/probe_sequences_EPIC.txt")

# Annotation
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19))

### Run probe_mapper, assuming the repeat is methylated, no mismatches
matches_methylated_nomismatch <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                              probe_sequences = probe_sequences,
                                              max_width = 35,
                                              min_width = 6, stepsize = 1,
                                              allow_mismatch = FALSE , allow_indel = FALSE, use_Y = FALSE,
                                              methylation_status = "methylated", nr_cores = 1)

matches_methylated_nomismatch <- clean_matches(matches_methylated_nomismatch, anno = anno, probe_sequences = probe_sequences)
matches_methylated_nomismatch_max <- select_max_width(matches_methylated_nomismatch)
matches_methylated_nomismatch_final <- select_unique_sequence(matches_methylated_nomismatch_max)
matches_methylated_nomismatch_final <- annotate(matches_methylated_nomismatch_final, anno = anno)
saveRDS(matches_methylated_nomismatch_final, file = "../../data/processed/c9_matches/EPIC_matches_methylated_nomismatch_final.rds")


### Run probe_mapper, assuming the repeat is unmethylated, no mismatches
matches_unmethylated_nomismatch <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                              probe_sequences = probe_sequences,
                                              max_width = 35,
                                              min_width = 6, stepsize = 1,
                                              allow_mismatch = FALSE , allow_indel = FALSE, use_Y = FALSE,
                                              methylation_status = "unmethylated", nr_cores = 1)

matches_unmethylated_nomismatch <- clean_matches(matches_unmethylated_nomismatch, anno = anno, probe_sequences = probe_sequences)
matches_unmethylated_nomismatch_max <- select_max_width(matches_unmethylated_nomismatch)
matches_unmethylated_nomismatch_final <- select_unique_sequence(matches_unmethylated_nomismatch_max)
matches_unmethylated_nomismatch_final <- annotate(matches_unmethylated_nomismatch_final, anno = anno)
saveRDS(matches_unmethylated_nomismatch_final, file = "../../data/processed/c9_matches/EPIC_matches_unmethylated_nomismatch_final.rds")

### Run probe_mapper, assuming the repeat is methylated, allowing 1 mismatch/INDEL

matches_methylated_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                           probe_sequences = probe_sequences,
                                           max_width = 35,
                                           min_width = 6, stepsize = 1,
                                           min_distance = 6,
                                           allow_mismatch = TRUE , allow_indel = TRUE, use_Y = FALSE,
                                           methylation_status = "methylated", nr_cores = 1)

matches_methylated_mismatch_indel <- clean_matches(matches_methylated_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_methylated_mismatch_indel_max <- select_max_width(matches_methylated_mismatch_indel)
matches_methylated_mismatch_indel_final <- select_unique_sequence(matches_methylated_mismatch_indel_max)
matches_methylated_mismatch_indel_final <- annotate(matches_methylated_mismatch_indel_final, anno = anno)
saveRDS(matches_methylated_mismatch_indel_final, file = "../../data/processed/c9_matches/EPIC_matches_methylated_mismatch_indel_final.rds")

### Run probe_mapper, assuming the repeat is unmethylated, allowing 1 mismatch/INDEL

matches_unmethylated_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                           probe_sequences = probe_sequences,
                                           max_width = 35,
                                           min_width = 6, stepsize = 1,
                                           min_distance = 6,
                                           allow_mismatch = TRUE , allow_indel = TRUE, use_Y = FALSE,
                                           methylation_status = "unmethylated", nr_cores = 1)

matches_unmethylated_mismatch_indel <- clean_matches(matches_unmethylated_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_unmethylated_mismatch_indel_max <- select_max_width(matches_unmethylated_mismatch_indel)
matches_unmethylated_mismatch_indel_final <- select_unique_sequence(matches_unmethylated_mismatch_indel_max)
matches_unmethylated_mismatch_indel_final <- annotate(matches_unmethylated_mismatch_indel_final, anno = anno)
saveRDS(matches_unmethylated_mismatch_indel_final, file = "../../data/processed/c9_matches/EPIC_matches_unmethylated_mismatch_indel_final.rds")

### Run probe_mapper, coding CpG-sites as YpG, no mismatch

matches_YpG_nomismatch <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                       probe_sequences = probe_sequences,
                                       max_width = 35,
                                       min_width = 6, stepsize = 1,
                                       allow_mismatch = FALSE , allow_indel = FALSE, use_Y = TRUE,
                                       nr_cores = 1)

matches_YpG_nomismatch <- clean_matches(matches_YpG_nomismatch, anno = anno, probe_sequences = probe_sequences)
matches_YpG_nomismatch_max <- select_max_width(matches_YpG_nomismatch)
matches_YpG_nomismatch_final <- select_unique_sequence(matches_YpG_nomismatch_max)
matches_YpG_nomismatch_final <- annotate(matches_YpG_nomismatch_final, anno = anno)
saveRDS(matches_YpG_nomismatch_final, file = "../../data/processed/c9_matches/EPIC_matches_YpG_nomismatch_final.rds")

### Run probe_mapper, coding CpG-sites as YpG, allow 1 mismatch/INDEL

matches_YpG_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                           probe_sequences = probe_sequences,
                                           max_width = 35,
                                           min_width = 6, stepsize = 1,
                                           min_distance = 6,
                                           allow_mismatch = TRUE , allow_indel = TRUE, use_Y = TRUE,
                                           nr_cores = 1)

matches_YpG_mismatch_indel <- clean_matches(matches_YpG_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_YpG_mismatch_indel_max <- select_max_width(matches_YpG_mismatch_indel)
matches_YpG_mismatch_indel_final <- select_unique_sequence(matches_YpG_mismatch_indel_max)
matches_YpG_mismatch_indel_final <- annotate(matches_YpG_mismatch_indel_final, anno = anno)
saveRDS(matches_YpG_mismatch_indel_final, file = "../../data/processed/c9_matches/EPIC_matches_YpG_mismatch_indel_final.rds")

## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
