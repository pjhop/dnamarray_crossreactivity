#!/bin/bash

~/osca_Linux --befile ../../data/processed/beta_cleaned_dasen_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking.txt \
--exclude-probe ../../data/temp/probes_miss0.05.txt \
--keep ../../data/temp/NL_2015.10_1000samples_samples.txt \
--mlma-loco --out ../../data/output/ewas/mlma_loco_c9_covarsex_qcovaragesmoking_missing0.05_NL_2015.10_1000samples > ../logs/mlma_loco_c9_covarsex_qcovaragesmoking_missing0.05_NL_2015.10_1000samples.log

