# Meta-analysis identifies shared microbial responses in various diseases and specific inflammatory bowel disease signal

This repo contains the jupyter notebooks, scripts and data for the creation of all figures not created using qiime2.

# How to use
Since some datasets used for the meta-analysis are not publicly available, datasets for the analysis are provided from the per disease cohort case/control effect sizes table (located [here](https://github.com/amnona/paper-metaanalysis/blob/main/ratios/ratios_no_bloom.biom) ). The python notebook for creating this table (without the per-study datasets) is located [here](https://github.com/amnona/paper-metaanalysis/blob/main/scripts/ratios-create-table.ipynb).

## Available notebooks description:
All jupyter notebooks for the various figures are located in the [scripts/](https://github.com/amnona/paper-metaanalysis/tree/main/scripts) directory

* ratios-create-table.ipynb: Create the normalized effect-size table from the per-study samples. Data for this script is not available since it includes per-sample datasets. The output of this script is available and is used in the ratios.ipynb analysis (see below).
* all-samples.ipynb: Merge and subsample all the per-sample data from all studies for the PCoA and PERMANOVA analysis. Data for this script is not available since it includes per-sample datasets.
* distances2.ipynb: Draw the inter-study distance matrix and distance distributions.
* dysbiosis.ipynb: Calculate the dysbiosis index and compare performance. Data for this script is not available since it includes per-sample datasets (for testing the performance on each sample).
* parse-classifier-results.ipynb : Process the classifier analysis results (which are run in parallel using [this bash script](https://github.com/amnona/paper-metaanalysis/blob/main/classifier/run_classifier_batch.sh) ). The notebook creates the AUC figures.
* picrust2.ipynb: Parse the results of the picrust2 prediction table.
* ratios.ipynb: Identify the non-specific and IBD-specific ASVs and create appropriate heatmaps.
* taxonomy.ipynb: Create the taxonomy bar-plots for the non-specific bacteria using the qiime2 assigned taxonomies.

## Setup
Using the scripts requires the latest release of the [Calour](https://github.com/biocore/calour) python module (see [installation instructions](https://github.com/biocore/calour/blob/master/INSTALL.md#install-the-latest-manually-from-github-repository) ) and the [calour_utils](https://github.com/amnona/calour_utils) python module (install from git).

# Questions/Comments?
please open an [Issue](https://github.com/amnona/paper-metaanalysis/issues) in this github repository

