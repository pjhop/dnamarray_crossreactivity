# Crossreactive probes on Illumina DNA methylation arrays: a large study on ALS shows that a cautionary approach is warranted in interpreting epigenome-wide association studies

Scripts accompanying the manuscript 'Crossreactive probes on Illumina DNA methylation arrays: a large study on ALS shows that a cautionary approach is warranted in interpreting epigenome-wide association studies'. This repository contains the scripts used for analyses. For convenience, we wrapped
several functions into an R package, which is available at: https://github.com/pjhop/DNAmCrosshyb.

## Data availability
The data used in downstream analyses can be downloaded from: https://sandbox.zenodo.org/record/684090. 
All downstream analyses in c9_analysis.Rmd (https://github.com/pjhop/dnamarray_crossreactivity/blob/master/analysis/c9_analysis.Rmd) and in supplementary_note.Rmd (https://github.com/pjhop/dnamarray_crossreactivity/blob/master/analysis/supplementary_note.Rmd) can be reproduced using the deposited data as follows:
- Clone this repository: <git clone  https://github.com/pjhop/dnamarray_crossreactivity.git>
- Download the data ('data.zip') and place it in the top directory of the dnamarray_crossreactivity folder. 
- Unzip the data.zip folder

## Overview
**R**: Contains some reusable scripts for plotting and mapping probes to sequences.  
**analysis**: Contains scripts used for the analyses presented in the paper  
* **preprocessing**: Folder that contains several scripts used for preprocessing the data. This includes
transforming the data to the format used in the *OSCA* software, PCA, extracting signal intensities
and normalizing out-of-band (OOB) signal intensities.  
* **ewas**: Folder that contains scripts to perform EWAS analyses, most of which is done using the *OSCA* software.  
Also contains scripts to perform meta-analyses.
* **c9_matches**: Folder that contains scripts to map 450k and EPIC probes to the *C9orf72* repeat expansion.
Includes scripts to do so in different ways, for example: assuming that the repeat is
methylated/unmethylated, representing CpGs by YpG, including flanking regions etc.
* **other**: Folder that contains script to match 450k/EPIC probes to various repeat sequences.  
* **c9_analysis.Rmd**: Main Rmarkdown file in which all downstream analyses are performed and figures/tables are generated.
* **supplementary_note.Rmd**: Rmarkdown that generates the supplementary note.

### R
**ewas_plot_functions.R**: Contains several plotting functions (e.g. qqplots, manhattan plots,
  comparing tstats). See the markdown-file for usage.  
**probe_mapper_nonref.R**: Map probes to a non-reference sequence.
This function is used in the scripts in the analysis/c9_matches folder.
Note that an updated version is included as a function in the `DNAmCrosshyb` package,
called `map_probes_sequence`.  
**create_binary_files.R**: R-wrapper to create binary files using the OSCA software.  

### analysis/c9_analysis.Rmd

### analysis/preprocessing

**create_binary_files_\***: Convert matrices containing beta-values/m-values/signal intensities
to binary format used in the *OSCA* software.  
**pca_\***: Perform PCA using *OSCA* software.  
**naten_\***: Extract and normalize out-of-band (OOB) signal intensities. The `get_OOB()` function from the `DNAmCrosshyb` package is used for this.  
**save_total_intensities_NL.R**: Extract total signal intensities and normalize.

### analysis/ewas

**mlma_loco_\***: Perform a mixed linear model (leave-one-chromosome out) as implemented in the *OSCA* software. Different sets of samples/covariates are used for sensitivity analyes (e.g. removing PCA outliers, using m-values instead of beta-values etc.).  
**mlma_linear_\***: Performed a linear model using the *OSCA* software.  
**save_covariates_\***: Save sets of covariates used in the different OSCA analyses.  
**meta_analysis.R**: Perform a meta-analysis (fixed effects inverse variance weighted) on
the three NL 450k batches.  

### analysis/c9_matches

**generate_probe_sequences.R**: Generate all combinations of probe sequences.
Since type II probes can contain 'R' bases that represent either an A or a G (IUPAC code),
there are several possible probe sequences for these probes.  
**get_c9matches_\***: Map 450k/EPIC probes to the C9 repeat expansion using various assumptions (e.g. allowing mismatches/INDELS, assume the repeat is methylated/unmethylated, include
  flanking regions etc.).  
**process_c9matches_\***: Some cleaning of the matches, includes removing duplicate matches (type I probes consist of two sequences, and type II probes of up to 8 sequences), selecting the longest match for each probe and annotating them with predicted color channels.  
