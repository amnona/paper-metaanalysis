# Meta-analysis identifies shared microbial responses in various diseases and specific inflammatory bowel disease signal

This repo contains the jupyter notebooks, scripts and data for the creation of all figures not created using qiime2.

# How to use
Since some datasets used for the meta-analysis are not publicly available, datasets for the analysis are provided from the per disease cohort case/control effect sizes table (located [here](https://github.com/amnona/paper-metaanalysis/blob/main/ratios/ratios_no_bloom.biom) ). The python notebook for creating this table (without the per-study datasets) is located [here](https://github.com/amnona/paper-metaanalysis/blob/main/scripts/ratios-create-table.ipynb).

All jupyter notebooks for the various figures are located in the [scripts/](https://github.com/amnona/paper-metaanalysis/tree/main/scripts) directory

## Setup
Using the scripts requires the latest release of the [Calour](https://github.com/biocore/calour) python module (see [installation instructions](https://github.com/biocore/calour/blob/master/INSTALL.md#install-the-latest-manually-from-github-repository) ) and the [calour_utils](https://github.com/amnona/calour_utils) python module (install from git).

# Questions/Comments?
please open an [Issue](https://github.com/amnona/paper-metaanalysis/issues) in this github repository

