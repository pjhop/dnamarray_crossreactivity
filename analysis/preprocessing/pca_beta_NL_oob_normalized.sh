#!/bin/bash

# Create ORM
~/osca_Linux --befile ../../data/processed/beta_NL_oob_normalized --make-orm \
--keep ../../data/temp/keep.txt --out ../../data/processed/orm_beta_NL_oob_normalized > ../logs/orm_beta_NL_oob_normalized.log

# Perform PCA
~/osca_Linux --pca  --orm ../../data/processed/orm_beta_NL_oob_normalized \
--out ../../data/processed/pca_beta_NL_oob_normalized >> ../logs/pca_beta_NL_oob_normalized.log
