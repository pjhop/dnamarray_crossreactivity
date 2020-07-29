#!/bin/bash

#$ -cwd
#$ -M p.j.hop-2@umcutrecht.nl
#$ -m beas
#$ -o ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_oob.o
#$ -e ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_oob.e

~/osca_Linux --befile ../../data/processed/beta_UK_2017.06_383samples_Ireland_2017.03_298samples_oob_normalized \
--pheno ../../data/temp/pheno_c9_EPIC.txt \
--covar ../../data/temp/covar_sex_batch_EPIC.txt \
--qcovar ../../data/temp/qcovar_predicted_age_cellcounts_5PCs_5controlPCs_EPIC.txt \
--linear --out ../../data/output/ewas/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_oob > ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_oob.log
