#!/bin/bash

~/osca_Linux --befile ../../data/processed/total_intensities_oob_normalized_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking_cellcounts_5controlPCs_5PCs_total_intensities_oob.txt \
--linear --out ../../data/output/ewas/mlma_linear_c9_covarsexbatch_qcovaragesmokingcellcounts5controlPCs5PCs_total_intensities_oob > ../logs/mlma_linear_c9_covarsexbatch_qcovaragesmokingcellcounts5controlPCs5PCs_total_intensities_oob.log
