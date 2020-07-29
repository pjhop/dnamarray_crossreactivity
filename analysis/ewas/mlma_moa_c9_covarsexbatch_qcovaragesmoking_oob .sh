#!/bin/bash

~/osca_Linux --befile ../../data/processed/beta_NL_oob_normalized \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking.txt \
--orm ../../data/processed/orm_beta_cleaned_dasen_NL \
--moa --out ../../data/output/ewas/mlma_moa_c9_covarsexbatch_qcovaragesmoking_oob > ../logs/mlma_moa_c9_covarsexbatch_qcovaragesmoking_oob.log
