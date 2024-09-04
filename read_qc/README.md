# Read evaluation pipeline

A pipeline is used to evaluate PacBio HiFi reads before feeding them into an assembler. It amalgamates reports from several independent evaluation methods.

Metrics evaluated:
* Basic read stats of the reads using [NanoPlot](https://github.com/wdecoster/NanoPlot)
* Contamination evaluation with [kraken2](https://github.com/DerrickWood/kraken2) and visualization with [krona](https://github.com/marbl/Krona), [kat](https://github.com/TGAC/KAT)
* K-mer counting with [KMC](https://github.com/refresh-bio/KMC)
* K-mer profile analysis with [genomescope2.0](https://github.com/tbenavi1/genomescope2.0), and [smudgeplot](https://github.com/KamilSJaron/smudgeplot)
* 

> [!NOTE]  
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inquire the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequisites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), [genomescope2.0](https://github.com/tbenavi1/genomescope2.0) should also be present and added to the `$PATH` variable.

Please create a `data/` directory where the different reads will be located:
* Concatenated PacBio HiFi reads in a **fasta** format should be placed in `data/reads`

For the contaminant screening the custom database should be present in the `bin/db/` directory. 

## Usage

With the complete `bin/` and all other prerequisites simply run:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The runs from different evaluation methods are grouped in `[method]_report/` directories. Within each directory multiple directories are found for each input assembly.
