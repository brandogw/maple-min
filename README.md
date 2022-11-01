# Maple-Min: Mutation Analysis for Parallel Evolution

Bare bones version of Gordon Rix's [MAPLE pipeline](https://github.com/gordonrix/maple). Only works with Illumina reads.

Maple is a [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline for analysis of
mutation-rich next generation sequencing data from highly parallelized targeted evolution experiments.

## Setup

Maple requires conda, which can be installed by following [these instructions](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).
Miniconda is lighter weight and provides all that is needed by Maple. Create a conda environment that enables usage of mamba, which is better capable of creating environments with snakemake:

    conda create mamba -n mambaEnv -c conda-forge


Next, use mamba to create the environment that we need. If using mamba within its own environment, the path to the environment should
be specified with the `--prefix In this example, the default path prefix ~/miniconda3/envs is used, but please ensure this is the correct location:

    conda activate mambaEnv
    mamba env create --prefix ~/miniconda3/envs/m2h-seq --file envs/environment.yml
    conda deactivate
    conda activate ~/miniconda3/envs/m2h-seq


Install additional packages

    sudo apt install build-essential libz-dev

Finally, run the following snakemake command to download all the binary files:

    snakemake --snakefile rules/install.smk --directory ~/miniconda3/envs/m2h-seq -j 4 all

Install basemount to grab files straight from the Illumina server. Follow the instrucions found on the website. You will need to log in and authenticate via the basespace website. This will then mount a virtual drive on your computer, which then you can grab files directly from basespace.

    sudo bash -c "$(curl -L https://basemount.basespace.illumina.com/install)"
    cd
    mkdir Basespace
    basemount Basespace/


## Best operating proceedure

All the files you want to manipulate are in the config folder. Mainly update the following:
- config/samples.csv: mainly this can contain more meta data later, right now list the samples to download and analyze in the csv.
- config.yaml: Make sure all field are updated
- ref/ : Make sure references are up to date

    snakemake --cores 8

