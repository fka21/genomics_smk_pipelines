# Single-end RNAseq alignment to transcriptome using pseudoaligner

A pipeline used to align short paired-end Illumina reads to a genome.

![](https://github.com/fka21/genomics_smk_pipelines/blob/main/rnaseq_pipelines/rnaseq_single2transcriptome/dag.svg)

Steps included in the pipeline:
* Initialization steps (symlinking databases, installing tools, indexing genome, conversion of GFF3 to GTF)
* Trimming of reads using [fastp]()
* Filtering reads for rRNA contamination using [sortmerna]() and for contaminants using [kraken2]().
* Alignment to the genome is performed using [STAR]()
* Multiple quality checks are also performed and gathered with [MultiQC]()
  

> **NOTE**  
> 
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequisites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), should also be present and added to the `$PATH` variable.

Please upload to `workflow/resources` directory where the different assemblies and reads will be located:
* Assemblies with a `.fa` suffix and GFF3 formatted annotation files should be placed in `workflow/resources/genome/` directory
* Paired-end reads should be stored in `workflow/resources/reads/`

> **NOTE**
> Please edit the `config/config.yml` with the desired genome and annotation file name. Also, check if the suffix in the config file is the same as the reads. Same holds true to the desire databases.

## Usage

With the complete `bin/` and all other prerequisites simply run:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The runs from different evaluation methods are grouped in `[method]_report/` directories. Within each directory multiple directories are found for each input assembly.

> **NOTE**
> 
> The `inspector_report/` outputs can be utilized to further correct  for local misassemblies. Please consult [inspector](https://github.com/Maggi-Chen/Inspector) documentation for further information.




