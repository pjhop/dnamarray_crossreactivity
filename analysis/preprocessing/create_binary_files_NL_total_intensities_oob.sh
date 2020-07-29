#!/bin/bash
module load R/3.5.1
Rscript ../../R/create_binary_files.R \
--beta ../../data/processed/total_intensities_oob_normalized_NL.rds \
--type other \
--sex sep \
--array 450k \
--opi ../../data/extdata/annotated_450k.opi \
--removesamples ../../data/raw_data/NL_drop_total.txt \
--out ../../data/processed/total_intensities_oob_normalized_NL
