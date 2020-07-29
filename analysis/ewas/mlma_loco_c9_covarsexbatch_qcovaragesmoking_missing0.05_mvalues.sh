#!/bin/bash

~/osca_Linux --befile ../../data/processed/mvalues_cleaned_dasen_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking.txt \
--exclude-probe ../../data/temp/probes_miss0.05.txt \
--mlma-loco --out ../../data/output/ewas/mlma_loco_c9_covarsexbatch_qcovaragesmoking_missing0.05_mvalues > ../logs/mlma_loco_c9_covarsexbatch_qcovaragesmoking_missing0.05_mvalues.log

