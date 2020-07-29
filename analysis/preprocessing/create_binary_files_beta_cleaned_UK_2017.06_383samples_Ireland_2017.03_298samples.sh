#!/bin/bash
module load R/3.5.1
Rscript ../../R/create_binary_files.R \
--beta ../../data/raw_data/beta_cleaned_UK_2017.06_383samples_Ireland_2017.03_298samples.rds \
--type beta \
--sex sep \
--array EPIC \
--opi ../../data/extdata/annotated_epic.opi \
--out ../../data/processed/beta_cleaned_UK_2017.06_383samples_Ireland_2017.03_298samples
