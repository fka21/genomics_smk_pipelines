# Suite of Snakemeake pipelines used during _de novo_ genome assemblies
Various tools for evaluating PacBio HiFi reads, quality checking of contigs generated by different assemblers and so on.

## Prerequisites
[Snakemake](https://snakemake.readthedocs.io/en/stable/) should be installed together with [Anaconda/Miniconda3](https://docs.anaconda.com/miniconda/).

Generally all required conda environments are implemented in each pipeline and are installed during deploying the pipeline.

The directoru structure is fixed and parts of it should be defined _a priori_(details in each pipeline can be found).

## Usage

Move to the working directory (i.e. where the Snakefile is located) and run the following command:

```
snakemake --cores [user-defined] --use-conda all
```