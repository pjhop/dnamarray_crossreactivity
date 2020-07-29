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
next_base = "G"
prev_base = "C"

# Load probe sequences
probe_sequences <- readr::read_tsv("../../data/processed/probe_sequences_450k.txt")

# Annotation
anno <- as.data.frame(getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19))


### Run probe_mapper, assuming the repeat is methylated, allowing 1 mismatch/INDEL

matches_methylated_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                                  probe_sequences = probe_sequences,
                                                  max_width = 35,
                                                  min_width = 6, stepsize = 1,
                                                  min_distance = 0,
                                                  allow_mismatch = TRUE , allow_indel = TRUE, use_Y = FALSE,
                                                  methylation_status = "methylated", nr_cores = 1)

matches_methylated_mismatch_indel <- clean_matches(matches_methylated_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_methylated_mismatch_indel_max <- select_max_width(matches_methylated_mismatch_indel)
matches_methylated_mismatch_indel_final <- select_unique_sequence(matches_methylated_mismatch_indel_max)
matches_methylated_mismatch_indel_final <- annotate(matches_methylated_mismatch_indel_final, anno = anno)
saveRDS(matches_methylated_mismatch_indel_final, file = "../../data/processed/c9_matches/450k_matches_methylated_mismatch_indel_0bp_final.rds")

### Run probe_mapper, assuming the repeat is unmethylated, allowing 1 mismatch/INDEL

matches_unmethylated_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                                    probe_sequences = probe_sequences,
                                                    max_width = 35,
                                                    min_width = 6, stepsize = 1,
                                                    min_distance = 0,
                                                    allow_mismatch = TRUE , allow_indel = TRUE, use_Y = FALSE,
                                                    methylation_status = "unmethylated", nr_cores = 1)

matches_unmethylated_mismatch_indel <- clean_matches(matches_unmethylated_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_unmethylated_mismatch_indel_max <- select_max_width(matches_unmethylated_mismatch_indel)
matches_unmethylated_mismatch_indel_final <- select_unique_sequence(matches_unmethylated_mismatch_indel_max)
matches_unmethylated_mismatch_indel_final <- annotate(matches_unmethylated_mismatch_indel_final, anno = anno)
saveRDS(matches_unmethylated_mismatch_indel_final, file = "../../data/processed/c9_matches/450k_matches_unmethylated_mismatch_indel_0bp_final.rds")

### Run probe_mapper, coding CpG-sites as YpG, allow 1 mismatch/INDEL

matches_YpG_mismatch_indel <- probe_mapper(repeat_sequence, next_base = next_base, prev_base = prev_base,
                                           probe_sequences = probe_sequences,
                                           max_width = 35,
                                           min_width = 6, stepsize = 1,
                                           min_distance = 0,
                                           allow_mismatch = TRUE , allow_indel = TRUE, use_Y = TRUE,
                                           nr_cores = 1)

matches_YpG_mismatch_indel <- clean_matches(matches_YpG_mismatch_indel, anno = anno, probe_sequences = probe_sequences)
matches_YpG_mismatch_indel_max <- select_max_width(matches_YpG_mismatch_indel)
matches_YpG_mismatch_indel_final <- select_unique_sequence(matches_YpG_mismatch_indel_max)
matches_YpG_mismatch_indel_final <- annotate(matches_YpG_mismatch_indel_final, anno = anno)
saveRDS(matches_YpG_mismatch_indel_final, file = "../../data/processed/c9_matches/450k_matches_YpG_mismatch_indel_0bp_final.rds")


## Print SessionInfo
cat("\n\n\nSessionInfo:\n\n")
sessionInfo()
