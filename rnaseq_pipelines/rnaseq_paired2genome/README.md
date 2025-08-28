# Paired RNAseq alignment to genome

A pipeline used to align short paired-end Illumina reads to a transcriptome.

![](https://github.com/fka21/genomics_smk_pipelines/blob/main/rnaseq_pipelines/rnaseq_paired2genome/dag.svg)

Steps included in the pipeline:
* Initialization steps (symlinking databases, installing tools, indexing genome, and GFF3 conversion)
* Trimming of reads using [fastp](https://github.com/OpenGene/fastp)
* Filtering reads for rRNA contamination using [sortmerna]() and for contaminants using [kraken2](https://github.com/DerrickWood/kraken2).
* Alignment to the transcriptome is performed using [STAR](https://github.com/alexdobin/STAR)
* Post alignment BAM file filtering (removing multimappers, low quality mappings using [samtools](https://samtools.github.io))
* Correction for GC biases of the alignment using [deeptools](https://deeptools.readthedocs.io/en/latest/)
* Multiple quality checks are also performed and gathered with [MultiQC](https://seqera.io/multiqc/)
* Count quantification on the filtered BAM files using [featureCounts](https://subread.sourceforge.net/featureCounts.html)
  

> **NOTE**  
> 
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequisites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), should also be present and added to the `$PATH` variable.

Please upload to `workflow/resources` directory where the different assemblies and reads will be located:
* Assemblies with a `.fa` suffix and GFF3 formatted annotation files should be placed in `workflow/resources/genome/` directory and the `config/config.yml` should be modified to include the new assembly.
* Paired-end reads should be stored in `workflow/resources/reads/`

> **NOTE**
> Please edit the `config/config.yml` with the desired genome and annotation file name. Also, check if the suffix in the config file is the same as the reads. Same holds true to the desire databases. This is also where you can set the read decontamination filter flags and the estimated genome size.

## Usage

Prior to running it is usefull to check if all prerequisites are in place:

```
snakemake -np
```

If everything seems to be in place, then use the following command:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The quantification files are located in the `05_featurecounts/` directory. MultiQC report are found in `00_qc/multiqc.html`.


