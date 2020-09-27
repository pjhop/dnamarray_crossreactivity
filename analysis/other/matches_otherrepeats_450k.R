library(DNAmCrosshyb)
## run for several sequences 
sequences <- list(
  ATTCT = list(
    sequence = paste(rep("ATTCT", 15), collapse=""),
    prev_base = "T",
    next_base = "A"
  ),
  CAG = list(
    sequence = paste(rep("CAG", 20), collapse=""),
    prev_base = "G",
    next_base = "C"
  ),
  CTG = list(
    sequence = paste(rep("CTG", 20), collapse=""),
    prev_base = "G",
    next_base = "C"
  ),
  GAA = list(
    sequence = paste(rep("GAA", 20), collapse=""),
    prev_base = "A",
    next_base = "G"
  ),
  GCC = list(
    sequence = paste(rep("GCC", 20), collapse=""),
    prev_base = "C",
    next_base = "G"
  ),
  GCG = list(
    sequence = paste(rep("GCG", 20), collapse=""),
    next_base = "G",
    prev_base = "G"
  ),
  CGG = list(
    sequence = paste(rep("CGG", 20), collapse=""),
    next_base = "C",
    prev_base = "G"
  )
)

run <- function(lst) {
  print(lst[["sequence"]])
  matches <- map_probes_sequence(sequence = lst[["sequence"]], next_base = lst[["next_base"]], prev_base = lst[["prev_base"]],
                                 array = "450k", min_width = 10, max_width = 35, allow_indel = TRUE, 
                                 allow_mismatch = TRUE, step_size = 1, use_Y = TRUE)
  matches
}

matches_all <- purrr::map(sequences, .f = run)
names(matches_all) <- names(sequences)
saveRDS(matches_all, file = "../../data/misc/matches_repeat_sequences_10_35_useY_allowmismatch_allowINDEL_450k.rds")



run <- function(lst) {
  matches <- map_probes_sequence(sequence = lst[["sequence"]], next_base = lst[["next_base"]], prev_base = lst[["prev_base"]],
                                 array = "450k", min_width = 10, max_width = 35, allow_indel = TRUE, 
                                 allow_mismatch = TRUE, step_size = 1, use_Y = FALSE, methylation_status = "methylated")
  matches
}

matches_all_methylated <- purrr::map(sequences, .f = run)
names(matches_all_methylated) <- names(sequences)
saveRDS(matches_all_methylated, file = "../../data/misc/matches_repeat_sequences_10_35_methylated_allowmismatch_allowINDEL_450k.rds")

run <- function(lst) {
  matches <- map_probes_sequence(sequence = lst[["sequence"]], next_base = lst[["next_base"]], prev_base = lst[["prev_base"]],
                                 array = "450k", min_width = 10, max_width = 35, allow_indel = TRUE, 
                                 allow_mismatch = TRUE, step_size = 1, use_Y = FALSE, methylation_status = "unmethylated")
  matches
}

matches_all_unmethylated <- purrr::map(sequences, .f = run)
names(matches_all_unmethylated) <- names(sequences)
saveRDS(matches_all_unmethylated, file = "../../data/misc/matches_repeat_sequences_10_35_unmethylated_allowmismatch_allowINDEL_450k.rds")
