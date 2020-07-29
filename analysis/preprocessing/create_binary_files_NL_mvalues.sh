#!/bin/bash
module load R/3.5.1
Rscript ../../R/create_binary_files.R \
--beta ../../data/raw_data/mvalues_cleaned_NL.rds \
--type mvalues \
--sex sep \
--array 450k \
--opi ../../data/extdata/annotated_450k.opi \
--removesamples ../../data/raw_data/NL_related_drop.txt \
--out ../../data/processed/mvalues_cleaned_dasen_NL
