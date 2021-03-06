### Functions to test whether 450k/EPIC probes match to a given bisulfite-converted DNA sequence
## Allows matching of specified 3' subsequence length of the probe (0-50)
## Includes option to allow mismatches/INDELs

probe_mapper <- function(sequence, next_base, prev_base, probe_sequences,max_width = 50, min_width = 15, stepsize = 5,
                    allow_mismatch = FALSE, max_mismatch = 1, allow_indel = FALSE, min_distance = 6,
                    use_Y = TRUE, methylation_status = "methylated", nr_cores = 1) {

  convert_vec <- c("C" = "G", "G" = "C", "A" = "T", "T" = "A", "Y" = "R", "R" = "Y", "?" = "?")

  ## In silico bisulfite conversion
  bs_converted_sequence_forward <- Biostrings::DNAString(bisulfite_convert(sequence, nextbase = next_base, use_Y = use_Y,
                                                                           methylation_status = methylation_status))
  bs_converted_sequence_reverse <- Biostrings::DNAString(bisulfite_convert(complement_strand(sequence),
                                                                           nextbase = unname(convert_vec[prev_base]),
                                                                            use_Y = use_Y, methylation_status = methylation_status))
  bs_converted_sequence_forward_complement <- Biostrings::reverseComplement(bs_converted_sequence_forward)
  bs_converted_sequence_reverse_complement <- Biostrings::reverseComplement(bs_converted_sequence_reverse)

  # Convert into biostrings DNAStringSet
  stringset <- Biostrings::DNAStringSet(probe_sequences$Probe_sequence)

  widths <- seq(from = min_width, to = max_width, by = stepsize)
  if(widths[length(widths)] != max_width) widths <- c(widths, max_width)

  # Add parallel option
  if(nr_cores > 1) {
    match_total <- BiocParallel::bplapply(X = widths, FUN = matches,
                                          stringset = stringset,
                                          forward = bs_converted_sequence_forward,
                                          reverse = bs_converted_sequence_reverse,
                                          forward_complement = bs_converted_sequence_forward_complement,
                                          reverse_complement = bs_converted_sequence_reverse_complement,
                                          allow_mismatch = allow_mismatch, max_mismatch = max_mismatch,
                                          allow_indel = allow_indel, min_distance = min_distance,
                                          BPPARAM = BiocParallel::MulticoreParam(nworkers = nr_cores))
    match_total <- dplyr::bind_rows(match_total)
  } else {
    match_total <- purrr::map_df(widths, .f = matches, stringset = stringset,
                                 forward = bs_converted_sequence_forward,
                                 reverse = bs_converted_sequence_reverse,
                                 forward_complement = bs_converted_sequence_forward_complement,
                                 reverse_complement = bs_converted_sequence_reverse_complement,
                                 allow_mismatch = allow_mismatch, max_mismatch = max_mismatch,
                                 allow_indel = allow_indel, min_distance = min_distance)
  }


  if(nrow(match_total) > 0) {
    if(!allow_mismatch && !allow_indel) {
      match_total$mismatch_pos <- match_total$width_incl_indel <- match_total$indel_pos <- NA
    } else if(!allow_mismatch && allow_indel) {
      match_total$width_incl_indel <- match_total$indel_pos <- NA
    } else if(allow_mismatch && !allow_indel) {
      match_total$mismatch_pos <- NA
    }

    # Add sequence that is matched
    match_total <- subseq_sequences_wrapper(match_total, forward = bs_converted_sequence_forward, reverse = bs_converted_sequence_reverse,
                                            forward_complement = bs_converted_sequence_forward_complement,
                                            reverse_complement = bs_converted_sequence_reverse_complement)

    # Add next base
    match_total <- get_nextbase(match_total, forward = bs_converted_sequence_forward, reverse = bs_converted_sequence_reverse,
                                forward_complement = bs_converted_sequence_forward_complement, reverse_complement = bs_converted_sequence_reverse_complement,
                                next_base = next_base, prev_base = prev_base, use_Y = use_Y, methylation_status = methylation_status)
  }
  # Return
  match_total
}

# width = 20
# allow_mismatch = TRUE
# stringset <- Biostrings::DNAStringSet(probe_sequences$Probe_sequence)

library(magrittr)
matches <- function(width, stringset, forward, reverse, forward_complement, reverse_complement,
                    allow_mismatch, max_mismatch, allow_indel, min_distance) {
   print(width)
   ## Create dictionary

    stringset <- Biostrings::subseq(stringset, start = 50 - width + 1, end = 50)
    dict <- Biostrings::PDict(Biostrings::reverseComplement(stringset))

    # Forward
    matches_forward <- Biostrings::matchPDict(dict, forward, fixed= "pattern")
    if(sum(IRanges::elementNROWS(matches_forward) > 0)) {
      positions_forward <- create_indices(matches_forward, strand = "forward")
      positions_forward$strand <- "forward"
      positions_forward$width <- width
    } else {
      positions_forward <- NULL
    }

    # Reverse
    matches_reverse <- Biostrings::matchPDict(dict, reverse, fixed= "pattern")
    if(sum(IRanges::elementNROWS(matches_reverse) > 0)) {
      positions_reverse <- create_indices(matches_reverse, strand = "reverse")
      positions_reverse$strand <- "reverse"
      positions_reverse$width <- width
    } else {
      positions_reverse <- NULL
    }

    # Forward complement
    matches_forward_complement <- Biostrings::matchPDict(dict, forward_complement, fixed= "pattern")
    if(sum(IRanges::elementNROWS(matches_forward_complement) > 0)) {
      positions_forward_complement <- create_indices(matches_forward_complement, strand = "forward_complement")
      positions_forward_complement$strand <- "forward_complement"
      positions_forward_complement$width <- width
    } else {
      positions_forward_complement <- NULL
    }

    # Reverse complement
    matches_reverse_complement <- Biostrings::matchPDict(dict, reverse_complement, fixed= "pattern")
    if(sum(IRanges::elementNROWS(matches_reverse_complement) > 0)) {
      positions_reverse_complement <- create_indices(matches_reverse_complement, strand = "reverse_complement")
      positions_reverse_complement$strand <- "reverse_complement"
      positions_reverse_complement$width <- width
    } else {
      positions_reverse_complement <- NULL
    }

    matches_total_nomismatch_noindel <- dplyr::bind_rows(positions_forward, positions_forward_complement,
                                                         positions_reverse, positions_reverse_complement)

    if(allow_mismatch) {
      dict <- Biostrings::PDict(Biostrings::reverseComplement(stringset), max.mismatch = 1)

      # Forward
      matches_forward <- Biostrings::matchPDict(dict, forward, fixed= "pattern", max.mismatch = 1)
      if(sum(IRanges::elementNROWS(matches_forward) > 0)) {
        positions_forward <- add_mismatch_position(matches_forward, "forward", forward,
                                                   matches_total_nomismatch_noindel, width = width, stringset = stringset)
      } else {
        positions_forward <- NULL
      }

      # Reverse
      matches_reverse <- Biostrings::matchPDict(dict, reverse, fixed= "pattern", max.mismatch = 1)
      if(sum(IRanges::elementNROWS(matches_reverse) > 0)) {
        positions_reverse <- add_mismatch_position(matches_reverse, "reverse", reverse,
                                                   matches_total_nomismatch_noindel, width = width, stringset = stringset)
      } else {
        positions_reverse <- NULL
      }

      # Forward complement
      matches_forward_complement <- Biostrings::matchPDict(dict, forward_complement, fixed= "pattern", max.mismatch = 1)
      if(sum(IRanges::elementNROWS(matches_forward_complement) > 0)) {
        positions_forward_complement <- add_mismatch_position(matches_forward_complement,
                                                              "forward_complement", forward_complement,
                                                               matches_total_nomismatch_noindel, width = width, stringset = stringset)
      } else {
        positions_forward_complement <- NULL
      }

      # Reverse complement
      matches_reverse_complement <- Biostrings::matchPDict(dict, reverse_complement, fixed= "pattern", max.mismatch = 1)
      if(sum(IRanges::elementNROWS(matches_reverse_complement) > 0)) {
        positions_reverse_complement <- add_mismatch_position(matches_reverse_complement,
                                                              "reverse_complement", reverse_complement,
                                                              matches_total_nomismatch_noindel, width = width, stringset = stringset)
      } else {
        positions_reverse_complement <- NULL
      }

      matches_total_mismatch <- dplyr::bind_rows(positions_forward, positions_reverse,
                                          positions_forward_complement, positions_reverse_complement)

      if(nrow(matches_total_mismatch) > 0) {
        matches_total_mismatch$j <- NULL
        matches_total_mismatch <- matches_total_mismatch %>% dplyr::filter(mismatch_pos >= min_distance)
      }

      matches_total <- dplyr::bind_rows(matches_total_nomismatch_noindel, matches_total_mismatch)
    } else {
      matches_total <- matches_total_nomismatch_noindel
    }

    if(allow_indel) {
      ## Forward
      matches_forward <- Biostrings::countPDict(Biostrings::reverseComplement(stringset), forward,
                                                fixed= "pattern", with.indels = TRUE,
                                                max.mismatch = 1)
      if(sum(matches_forward) > 0) {
        positions_forward <- add_indel_position(matches_forward, stringset = stringset, strand = "forward",
                                                bs_genome = forward, width = width, min_distance = min_distance)
      } else {
        positions_forward <- NULL
      }

      ## Reverse
      matches_reverse <- Biostrings::countPDict(Biostrings::reverseComplement(stringset), reverse,
                                                fixed= "pattern", with.indels = TRUE,
                                                max.mismatch = 1)
      if(sum(matches_reverse) > 0) {
        positions_reverse <- add_indel_position(matches_reverse, stringset = stringset, strand = "reverse",
                                                bs_genome = reverse, width = width, min_distance = min_distance)
      } else {
        positions_reverse <- NULL
      }

      ## Forward complement
      matches_forward_complement <- Biostrings::countPDict(Biostrings::reverseComplement(stringset), forward_complement,
                                                fixed= "pattern", with.indels = TRUE,
                                                max.mismatch = 1)
      if(sum(matches_forward_complement) > 0) {
        positions_forward_complement <- add_indel_position(matches_forward_complement, stringset = stringset, strand = "forward_complement",
                                                bs_genome = forward_complement, width = width, min_distance = min_distance)
      } else {
        positions_forward_complement <- NULL
      }

      ## Reverse complement
      matches_reverse_complement <- Biostrings::countPDict(Biostrings::reverseComplement(stringset), reverse_complement,
                                                           fixed= "pattern", with.indels = TRUE,
                                                           max.mismatch = 1)
      if(sum(matches_reverse_complement) > 0) {
        positions_reverse_complement <- add_indel_position(matches_reverse_complement, stringset = stringset, strand = "reverse_complement",
                                                           bs_genome = reverse_complement, width = width, min_distance = min_distance)
      } else {
        positions_reverse_complement <- NULL
      }

      ## Merge
      matches_total_indel <- dplyr::bind_rows(positions_forward, positions_reverse,
                                                 positions_forward_complement, positions_reverse_complement)
      if(nrow(matches_total_indel) > 0) {
        matches_total_indel$j <- matches_total_indel$id_temp <- NULL
        matches_total_indel <- matches_total_indel %>% dplyr::mutate(id = paste(i, start, end, strand, sep = "_"))
        matches_total <- dplyr::bind_rows(matches_total, matches_total_indel)
      }

    }
    if(!allow_mismatch && !allow_indel) matches_total <- matches_total_nomismatch_noindel

    ## Add probe to dataframe
    if(nrow(matches_total) > 0) {
      matches_total <- matches_total %>%
        dplyr::left_join(probe_sequences[,c("i", "Probe")], by = "i") %>%
        dplyr::select(Probe, i, dplyr::everything())
    }
      matches_total


}




create_indices <- function(match, strand, j = FALSE) {
  index_matches <- IRanges::elementNROWS(match)
  index_matches <- tibble::tibble(i = 1:length(index_matches), nr_matches = index_matches)
  index_matches <- index_matches %>% dplyr::filter(nr_matches > 0)
  new_index <- rep(index_matches$i, times = index_matches$nr_matches)
  starts <- unlist(Biostrings::startIndex(match))
  ends <- unlist(Biostrings::endIndex(match))
  if(j) {
    index_matches2 <- tibble::tibble(i = new_index, j = 1:length(new_index), start = starts, end = ends)
  } else {
    index_matches2 <- tibble::tibble(i = new_index, start = starts, end = ends)
  }
  # Add identifier
  index_matches2$strand = strand
  index_matches2 <- index_matches2 %>% dplyr::mutate(id = paste(i, start, end, strand, sep = "_"))
  index_matches2
}


get_mismatch_position <- function(start_position, positions, bs_genome, stringset) {
  positions_tmp <- positions %>% dplyr::filter(start == start_position)
  bs_genome_tmp <- bs_genome[unique(positions_tmp$start):unique(positions_tmp$end)]
  stringset_tmp <- stringset[positions_tmp$i]
  check <- matrix(!Biostrings::hasLetterAt(Biostrings::reverseComplement(stringset_tmp),
                               bs_genome_tmp,
                               at = 1:length(bs_genome_tmp), fixed = FALSE),
                  nrow = length(stringset_tmp),
                  ncol = stringr::str_length(bs_genome_tmp))
  mismatch_pos <- apply(check, 1, FUN = which)
  positions_tmp$mismatch_pos <- mismatch_pos
  positions_tmp
}

add_mismatch_position <- function(matches, strand, bs_genome, matches_total, width, stringset) {
  # Get positions
  positions <- create_indices(matches, j = TRUE, strand = strand)
  if(nrow(matches_total) > 0) {
    positions <- positions %>% dplyr::filter(!(id %in% matches_total$id)) # Should not be a full match
  }
  positions <- positions %>%
    dplyr::mutate(width = width) %>% # Add width column
    dplyr::filter(!(end > length(bs_genome) | start == 0))  # remove mismatch at last position or at zero position

  if(nrow(positions) > 0) {
    # Get positions
    positions_mismatch <- purrr::map_df(unique(positions$start), .f = get_mismatch_position,
                                        positions = positions, bs_genome = bs_genome, stringset = stringset)
    positions_mismatch <- positions_mismatch %>% dplyr::filter(mismatch_pos != width)
  } else {
    positions_mismatch <- NULL
  }

  positions_mismatch
}

get_indel_matches <- function(i, stringset, bs_genome, strand, width) {
  match <- Biostrings::matchPattern(stringset[[i]], bs_genome,
                                    fixed = FALSE, with.indels = TRUE, max.mismatch = 1)
  matches <- tibble::tibble(i = i, j = 1:length(match), start = Biostrings::start(match),
                             end = Biostrings::end(match), strand = strand, width = width,
                             width_incl_indel = end - start + 1)
  matches
}

add_indel_position <- function(matches, stringset, strand, bs_genome, width, min_distance) {
  i_match <- which(matches > 0)
  ## Get all matches
  matches_check <- purrr::map_df(i_match, .f = get_indel_matches,
                                         stringset = Biostrings::reverseComplement(stringset),
                                         bs_genome = bs_genome, strand = strand, width = width)
  matches_check <- matches_check %>%
    dplyr::filter(width_incl_indel != width) %>%
    dplyr::mutate(id_temp = paste(start, end, sep = "_"))

  positions_indel <- purrr::map_df(unique(matches_check$id_temp), .f = get_indel_position,
                                   positions = matches_check, bs_genome = bs_genome, width = width, stringset = stringset)
  positions_indel <- positions_indel %>% dplyr::filter(indel_pos >= min_distance)
  positions_indel
}

get_indel_position <- function(id, positions, bs_genome, width, stringset) {
  positions_tmp <- positions %>% dplyr::filter(id_temp == id)
  bs_genome_tmp <- bs_genome[unique(positions_tmp$start):unique(positions_tmp$end)]
  stringset_tmp <- stringset[positions_tmp$i]

  # Note doesn't check last nucleotide if width < specified width!
  ## If no mismatches - means it at the last, when width < specified width
  check <- matrix(!Biostrings::hasLetterAt(Biostrings::reverseComplement(stringset_tmp),
                               bs_genome_tmp,
                               at = 1:length(bs_genome_tmp), fixed = FALSE),
                  nrow = length(stringset_tmp),
                  ncol = stringr::str_length(bs_genome_tmp))
  indel_pos <- suppressWarnings(apply(check, 1, FUN = function(x) min(which(x))))
  #last_indel <- suppressWarnings(apply(check, 1, FUN = function(x) max(which(x))))
  #total_mismatches <- apply(check, 1, FUN = function(x) sum(x))
  positions_tmp$indel_pos <- indel_pos
  if(unique(positions_tmp$width_incl_indel) < width) {
    positions_tmp <- positions_tmp %>% dplyr::filter(indel_pos != Inf)
  }
  positions_tmp
}

subseq_sequences_wrapper <- function(matches, forward, reverse, forward_complement, reverse_complement) {

  # Get unique start_end positions
  matches <- matches %>% dplyr::mutate(start_end_strand = paste(start, end, strand, sep = "_"))

  test <- matches %>% dplyr::filter(!duplicated(start_end_strand)) %>% dplyr::select(start, end, strand, start_end_strand)
  test <- test %>%
    dplyr::mutate(sequence_bs = dplyr::case_when(
      strand == "forward" ~ purrr::map2_chr(.x = start, .y = end, .f = subseq_sequences, sequence = forward),
      strand == "reverse" ~ purrr::map2_chr(.x = start, .y = end, .f = subseq_sequences, sequence = reverse),
      strand == "forward_complement" ~ purrr::map2_chr(.x = start, .y = end, .f = subseq_sequences, sequence = forward_complement),
      strand == "reverse_complement" ~ purrr::map2_chr(.x = start, .y = end, .f = subseq_sequences, sequence = reverse_complement)
    ))

  matches <- matches %>% dplyr::left_join(test[,c("start_end_strand", "sequence_bs")], by = "start_end_strand") %>%
                         dplyr::select(-start_end_strand)
  matches
}

get_nextbase <- function(matches, forward, reverse, forward_complement, reverse_complement,
  next_base, prev_base, use_Y, methylation_status) {

  # Bisulfite-convert
  convert_vec_bisulfite <- c("C" = "T", "G" = "G", "A" = "A", "T" = "T")
  convert_vec_bisulfite_Y <- c("C" = "Y", "G" = "G", "A" = "A", "T" = "T")
  convert_vec_bisulfite_methylated <- c("C" = "C", "G" = "G", "A" = "A", "T" = "T")
  convert_vec_bisulfite_unmethylated <- c("C" = "T", "G" = "G", "A" = "A", "T" = "T")
  convert_vec_complement <- c("C" = "G", "G" = "C", "A" = "T", "T" = "A",
                              "Y" = "R", "R" = "Y", "?" = "?")

  # Get unique start_end positions
  matches <- matches %>% dplyr::mutate(start_end_strand = paste(start, end, strand, sep = "_"))

  test <- matches %>% dplyr::filter(!duplicated(start_end_strand)) %>%
    dplyr::select(start, end, strand, start_end_strand)
  test_1 <- test %>% dplyr::filter(start == 1)
  test_2 <- test %>% dplyr::filter(start != 1)

  if(use_Y) {
    test_1 <- test_1 %>%
      dplyr::mutate(sbe_site = dplyr::case_when(
        start == 1 & strand == "forward" & as.character(forward[1]) == "G" ~ convert_vec_bisulfite_Y[prev_base],
        start == 1 & strand == "reverse" & as.character(reverse[1]) == "G" ~ unname(convert_vec_bisulfite_Y[convert_vec_complement[next_base]]),
        start == 1 & strand == "forward" ~ unname(convert_vec_bisulfite[prev_base]),
        start == 1 & strand == "reverse" ~ unname(convert_vec_bisulfite[convert_vec_complement[next_base]]),
        start == 1 & strand == "forward_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[next_base]]),
        start == 1 & strand == "reverse_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[prev_base]])
      ))
  } else if (methylation_status == "methylated") {
    test_1 <- test_1 %>%
      dplyr::mutate(sbe_site = dplyr::case_when(
        start == 1 & strand == "forward" & as.character(forward[1]) == "G" ~ prev_base,
        start == 1 & strand == "reverse" & as.character(reverse[1]) == "G" ~ unname(convert_vec_complement[next_base]),
        start == 1 & strand == "forward" ~ unname(convert_vec_bisulfite[prev_base]),
        start == 1 & strand == "reverse" ~ unname(convert_vec_bisulfite[convert_vec_complement[next_base]]),
        start == 1 & strand == "forward_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[next_base]]),
        start == 1 & strand == "reverse_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[prev_base]])
      ))
  } else if (methylation_status == "unmethylated") {
    test_1 <- test_1 %>%
      dplyr::mutate(sbe_site = dplyr::case_when(
        start == 1 & strand == "forward" & as.character(forward[1]) == "G" ~ convert_vec_bisulfite_unmethylated[prev_base],
        start == 1 & strand == "reverse" & as.character(reverse[1]) == "G" ~ unname(convert_vec_bisulfite_unmethylated[convert_vec_complement[next_base]]),
        start == 1 & strand == "forward" ~ unname(convert_vec_bisulfite[prev_base]),
        start == 1 & strand == "reverse" ~ unname(convert_vec_bisulfite[convert_vec_complement[next_base]]),
        start == 1 & strand == "forward_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[next_base]]),
        start == 1 & strand == "reverse_complement" ~ unname(convert_vec_complement[convert_vec_bisulfite_Y[prev_base]])
      ))
  }


  test_2 <- test_2 %>%
    dplyr::mutate(sbe_site = dplyr::case_when(
      (start != 1 & strand == "forward") ~ purrr::map2_chr(.x = start - 1, .y = start - 1, .f = subseq_sequences, sequence = forward),
      (start != 1 & strand == "reverse")  ~ purrr::map2_chr(.x = start - 1, .y = start - 1, .f = subseq_sequences, sequence = reverse),
      (start != 1 & strand == "forward_complement")  ~ purrr::map2_chr(.x = start - 1, .y = start - 1, .f = subseq_sequences, sequence = forward_complement),
      (start != 1 & strand == "reverse_complement") ~ purrr::map2_chr(.x = start - 1, .y = start - 1, .f = subseq_sequences, sequence = reverse_complement)
    ))

  test <- dplyr::bind_rows(test_1, test_2)

  matches <- matches %>% dplyr::left_join(test[,c("start_end_strand", "sbe_site")], by = "start_end_strand") %>%
    dplyr::select(-start_end_strand)
  matches
}

subseq_sequences <- function(x, y, sequence) {
  sub <- as.character(Biostrings::subseq(sequence, start = x, end = y))
  sub
}

## Get 'cleaned' matches (i.e. remove multiple i's per probe))
clean_matches <- function(matches, anno, probe_sequences) {

  matches <- matches %>%
    dplyr::left_join(anno[,c("Name", "Type")], by = c("Probe" = "Name")) %>%
    dplyr::left_join(probe_sequences[,c("i", "ProbeSeq")], by = "i")

  matches <- matches %>% dplyr::mutate(
    Type2 = dplyr::case_when(
      Type == "I" & ProbeSeq == "ProbeSeqA" ~ "I_Unmethylated",
      Type == "I" & ProbeSeq == "ProbeSeqB" ~ "I_Methylated",
      TRUE ~ "II"
    )
  )
  # Add ID
  matches <- matches %>%
                dplyr::mutate(probe_start_end_strand = paste(Probe, Type2, start, end, strand, sep = "_"))

  # Check
  matches <- matches %>%
    dplyr::select(-c(i, id)) %>%
    dplyr::distinct()

  matches
}


## Select match with highest width for each probe
select_max_width <- function(matches) {

  # Add ID
  matches <- matches %>%
    dplyr::group_by(Probe, Type2) %>%
    dplyr::filter(width == max(width)) %>%
    dplyr::ungroup()

  matches

}

select_unique_sequence <- function(matches) {

  # Add ID
  matches <- matches %>%
    dplyr::group_by(Probe, Type2, width, strand,sbe_site, sequence_bs) %>%
    dplyr::arrange(start) %>%
    dplyr::summarize(
                     start = start[1],
                     end = end[1],
                     indel_pos = indel_pos[1],
                     width_incl_indel = width_incl_indel[1],
                     mismatch_pos = mismatch_pos[1]) %>%
    dplyr::ungroup()

  matches
}

# TODO: add checks (mostly same as the complement_strand function)
complement_strand <- function(seq, format = c("string", "vector"), reverse = TRUE) {
  stopifnot(is.character(seq))
  stopifnot(is.logical(reverse))
  format <- match.arg(format)

  if(all(!stringr::str_detect(seq, "^[ATCGYR]+$"))) {
    stop("Basepairs should be A, C, T, Y or R")
  }
  convert_vec <- c("C" = "G", "G" = "C", "A" = "T", "T" = "A", "Y" = "R", "R" = "Y", "?" = "?")
  if(format == "string") {
    string <- stringr::str_split(seq, pattern = "")[[1]]
    complement <- stringr::str_c(convert_vec[string], collapse = "")
  } else if (format == "vector") {
    complement <- convert_vec[string]
    names(complement) = NULL
  }
  if(reverse & format == "string") {
    stringi::stri_reverse(complement)
  } else if(reverse & format == "vector") {
    rev(complement)
  } else {
    complement
  }
}

bisulfite_convert <- function(sequence, nextbase, methylation_status = c("methylated", "unmethylated"),
                              use_Y = FALSE) {
  stopifnot(is.character(sequence))
  methylation_status <- match.arg(methylation_status)

  convert_vec <- c("C" = "T", "G" = "G", "A" = "A", "T" = "T")
  sequence <- stringr::str_split(sequence, pattern = "")[[1]]

  if(use_Y) {
    sequence_bs <- vector("character", length = length(sequence))
    for(i in seq_along(sequence)) {
      if(i == length(sequence) && sequence[i] == 'C') {
        if(nextbase == "G") {
          sequence_bs[i] <- "Y"
        } else {
          sequence_bs[i] <- "T"
        }
      } else if(sequence[[i]] %in% c("T", "A", "G")) {
        sequence_bs[[i]] <- unname(convert_vec[sequence[[i]]])
      } else if(sequence[[i]] == "C" && dplyr::lead(sequence)[[i]] != "G") {
        sequence_bs[[i]] <- "T"
      } else {
        sequence_bs[[i]] <- "Y"
      }
    }
    sequence_bs <- stringr::str_c(sequence_bs, collapse = "")
    sequence_bs
  } else if(methylation_status == "methylated") {
    sequence_bs <- vector("character", length = length(sequence))
    for(i in seq_along(sequence)) {
      if(i == length(sequence) && sequence[i] == 'C') {
        if(nextbase == "G") {
          sequence_bs[i] <- "C"
        } else {
          sequence_bs[i] <- "T"
        }
      } else if(sequence[[i]] %in% c("T", "A", "G")) {
        sequence_bs[[i]] <- unname(convert_vec[sequence[[i]]])
      } else if(sequence[[i]] == "C" && dplyr::lead(sequence)[[i]] != "G") {
        sequence_bs[[i]] <- "T"
      } else {
        sequence_bs[[i]] <- "C"
      }
    }
    sequence_bs <- stringr::str_c(sequence_bs, collapse = "")
    sequence_bs
  } else if(methylation_status == 'unmethylated') {
    sequence_bs <- convert_vec[sequence]
    sequence_bs <- stringr::str_c(sequence_bs, collapse = "")
    sequence_bs
  }
}
