# EREfinder pipeline

This repository contains an R pipeline to:

1. retrieve one promoter/region per gene from Ensembl
2. run EREfinder
3. summarize ERE-like motif hits
4. compare genes of interest against random gene sets
5. visualize permutation results

## External dependency: EREfinder

This pipeline requires the EREfinder command-line executable.
EREfinder can be downloaded here: https://github.com/JonesLabIdaho/EREfinder

Install EREfinder separately, then provide the path to the executable:

erefinder_path <- "/path/to/EREfinder"

## Main script

The full pipeline is in:

`R/erefinder_pipeline_full.R`

## Example

See:

`examples/example_run.R`
