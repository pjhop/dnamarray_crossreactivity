#!/bin/bash

# Create ORM
~/osca_Linux --befile ../../data/processed/total_intensities_oob_normalized_NL --make-orm \
--keep ../../data/temp/keep.txt --out ../../data/processed/orm_total_intensities_oob_NL > ../logs/orm_total_intensities_oob_NL.log

# Perform PCA
~/osca_Linux --pca  --orm ../../data/processed/orm_total_intensities_oob_NL \
--out ../../data/processed/pca_total_intensities_oob_NL >> ../logs/pca_total_intensities_oob_NL.log
