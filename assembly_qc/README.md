# Assembly evaluation pipeline

A pipeline used to evaluate _de novo_ genome assemblies. Theoretically can be used on both contigs and scaffolds. It amalgamates reports from several independent evaluation methods. It is tailored towards the evaluation of assemblies generated using solely PacBio HiFi data.

Metrics evaluated:
* Genome gene complement completeness by utilizing [BUSCO](https://busco.ezlab.org/)
* Genome k-mer completeness by utilizing [merqury](https://github.com/marbl/merqury)
* Genome correctness based on read evidence by utilizing [inspector](https://github.com/Maggi-Chen/Inspector)
* Genome contamination levels by utilizing NCBI's Foreign Contamination Screening ([fcs](https://github.com/ncbi/fcs))

> [!NOTE]  
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequisites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), [inspector](https://github.com/Maggi-Chen/Inspector) should also be present and added to the `$PATH` variable.

Please create a `data/` directory where the different assemblies and reads will be located:
* Assemblies with a `.fa` suffix should be placed in `data/assemblies/` directory
* Concatenated PacBio HiFi reads in a **fasta** format should be placed in `data/reads`

For the contaminant screening the database should be present in the `bin/db/` directory. This can be done by using one of the scripts (following the [fcs-gx](https://github.com/ncbi/fcs/wiki/FCS-GX-input#fcs-gx-database-location) guide):

```
SOURCE_DB_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest"
LOCAL_DB="/path/to/db/folder"
python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/test-only" 
```

The BUSCO lineage should be *a priori* downloaded and be present in the working directory of the Snakefile with `busco_downloads/lineages/[user_determined_lineage]`.

> [!CAUTION]
> The pipeline will only choose the first directory if multiple BUSCO lineage directories exist within the path.

One can access lineages with the following commands.
```
# To view available lineages
busco --list-datasets

# To download the desired lineage
busco --download [user_determined_lineage]
```

## Usage

With the complete `bin/` and all other prerequisites simply run:

```
snakemake --cores [user-defined] --use-conda all
```

## Output

The runs from different evaluation methods are grouped in `[method]_report/` directories. Within each directory multiple directories are found for each input assembly.
