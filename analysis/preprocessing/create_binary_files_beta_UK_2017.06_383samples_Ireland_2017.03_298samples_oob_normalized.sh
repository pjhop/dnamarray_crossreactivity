#!/bin/bash
module load R/3.5.1
Rscript ../../R/create_binary_files.R \
--beta ../../data/processed/beta_UK_2017.06_383samples_Ireland_2017.03_298samples_oob_normalized.rds \
--type beta \
--sex sep \
--array EPIC \
--opi ../../data/extdata/annotated_epic.opi \
--out ../../data/processed/beta_UK_2017.06_383samples_Ireland_2017.03_298samples_oob_normalized
