# Assembly evaluation pipeline

A pipeline used to scaffold contigs using HiC data. It uses [bwa](https://github.com/lh3/bwa) for aligning HiC reads to the contigs, performs some post-alignment steps (see image below) and using [YAHS](https://github.com/c-zhou/yahs) it performs the scaffolding step. It also does some QC steps post-scaffolding.

![](https://github.com/fka21/genome_assembly_smk_pipelines/blob/main/scaffolding/dag.svg)

Metrics evaluated:
* Using [PretextView](https://github.com/sanger-tol/PretextView) a snapshot of HiC interactions is done.
* Using [hic-qc](https://github.com/phasegenomics/hic_qc) the HiC libraries are assessed. 

> **NOTE**
>   
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequisites

[Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/) should also be present and added to the `$PATH` variable.


## Usage

With the complete `bin/` and all other prerequisites simply run:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The runs from different evaluation methods are grouped in directories. Scaffolded assemblies are found in `02_scaffolding/` directory. To save space `.bam` files which are not required for scaffolding are removed after the pipeline is finished successfully.

> **NOTE**
> 
> [PretextView](https://github.com/sanger-tol/PretextView) is used to generate `.pretext` files in `05_pretext/` which in turn can be used to curate the scaffolding manually.
