#!/bin/bash

~/osca_Linux --befile ../../data/processed/beta_cleaned_dasen_NL \
--pheno ../../data/temp/pheno_c9.txt \
--covar ../../data/temp/covar_sex_batch.txt \
--qcovar ../../data/temp/qcovar_predictedage_smoking_cellcounts_5controlPCs_5PCs.txt \
--exclude-probe ../../data/temp/probes_miss0.05.txt \
--linear --out ../../data/output/ewas/mlma_linear_c9_covarsexbatch_qcovaragesmokingcellcounts5controlPCs5PCs_missing0.05 > ../logs/mlma_linear_c9_covarsexbatch_qcovaragesmokingcellcounts5controlPCs5PCs_missing0.05.log

