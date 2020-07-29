#!/bin/bash

~/osca_Linux --befile ../../data/processed/beta_cleaned_dasen_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking.txt \
--exclude-probe ../../data/temp/probes_miss0.05.txt \
--moa --out ../../data/output/ewas/mlma_moa_c9_covarsexbatch_qcovaragesmoking_missing0.05 > ../logs/mlma_moa_c9_covarsexbatch_qcovaragesmoking_missing0.05.log
