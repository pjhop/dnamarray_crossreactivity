#!/bin/bash
module load R/3.5.1
Rscript ../../R/create_binary_files.R \
--beta ../../data/raw_data/beta_cleaned_NL.rds \
--type beta \
--sex sep \
--array 450k \
--opi ../../data/extdata/annotated_450k.opi \
--removesamples ../../data/raw_data/NL_related_drop.txt \
--out ../../data/processed/beta_cleaned_dasen_NL
