#!/bin/bash

#$ -cwd
#$ -M p.j.hop-2@umcutrecht.nl
#$ -m beas
#$ -o ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_missing0.05.o
#$ -e ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_missing0.05.e

~/osca_Linux --befile ../../data/processed/beta_cleaned_UK_2017.06_383samples_Ireland_2017.03_298samples \
--pheno ../../data/temp/pheno_c9_EPIC.txt \
--covar ../../data/temp/covar_sex_batch_EPIC.txt \
--qcovar ../../data/temp/qcovar_predicted_age_cellcounts_5PCs_5controlPCs_EPIC.txt \
--exclude-probe ../../data/temp/probes_miss0.05_EPIC.txt \
--linear --out ../../data/output/ewas/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_missing0.05 > ../logs/mlma_EPIC_linear_c9_covarsexbatch_qcovaragecellcounts5PCs5controlPCs_missing0.05.log
