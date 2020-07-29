#!/bin/bash

~/osca_Linux --befile ../../data/processed/total_intensities_oob_normalized_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking.txt \
--mlma-loco --out ../../data/output/ewas/mlma_loco_c9_covarsexbatch_qcovaragesmoking_total_intensities_oob > ../logs/mlma_loco_c9_covarsexbatch_qcovaragesmoking_total_intensities_oob.log
