#!/bin/bash

~/osca_Linux --befile ../../data/processed/beta_cleaned_UK_2017.06_383samples_Ireland_2017.03_298samples \
--pheno ../../data/temp/pheno_c9_EPIC.txt \
--covar ../../data/temp/covar_sex_batch_EPIC.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking_EPIC.txt \
--exclude-probe ../../data/temp/probes_miss0.05_EPIC.txt \
--mlma-loco --out ../../data/output/ewas/mlma_EPIC_loco_c9_covarsexbatch_qcovaragesmoking_missing0.05 > ../logs/mlma_EPIC_loco_c9_covarsexbatch_qcovaragesmoking_missing0.05.log
