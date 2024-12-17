# Single-end RNAseq alignment to genome

A pipeline used to align short single-end Illumina reads to a transcriptome.

![](https://github.com/fka21/genomics_smk_pipelines/blob/main/rnaseq_pipelines/rnaseq_single2genome/dag.png)

Steps included in the pipeline:
* Initialization steps (symlinking databases, installing tools, indexing genome, and GFF3 conversion)
* Trimming of reads using [fastp](https://github.com/OpenGene/fastp)
* Filtering reads for rRNA contamination using [sortmerna]() and for contaminants using [kraken2](https://github.com/DerrickWood/kraken2).
* Alignment to the transcriptome is performed using [STAR](https://github.com/alexdobin/STAR)
* Multiple quality checks are also performed and gathered with [MultiQC](https://seqera.io/multiqc/)
  

> **NOTE 1**  
> 
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.

> **NOTE 2**  
> 
> The example files are empty, please use your own files for testing and don't forget to setup the config.yml and the Snakefile for it.

## Prerequisites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), should also be present and added to the `$PATH` variable.

Please upload to `workflow/resources` directory where the different assemblies and reads will be located:
* Assemblies with a `.fa` suffix and GFF3 formatted annotation files should be placed in `workflow/resources/genome/` directory and the `config/config.yml` should be modified to include the new assembly.
* Paired-end reads should be stored in `workflow/resources/reads/`

> **NOTE**
> Please edit the `config/config.yml` with the desired genome and annotation file name. Also, check if the suffix in the config file is the same as the reads. Same holds true to the desire databases.

## Usage

With the complete `bin/` and all other prerequisites simply run:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The quantification files are located in the `05_featurecounts/` directory. MultiQC report are found in `00_qc/multiqc.html`.


