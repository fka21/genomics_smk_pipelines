# Assembly evaluation pipeline

A pipeline used to evaluate _de novo_ genome assemblies. It amalgamates reports from several independent evaluation methods.

Metrics evaluated:
* Genome gene complement completeness by utilsing [BUSCO](https://busco.ezlab.org/)
* Genome k-mer completeness by utilising [merqury](https://github.com/marbl/merqury)
* Genome correctness based on read evidence by utilising [inspector](https://github.com/Maggi-Chen/Inspector)
* Genome contamination levels by utilising NCBI's Foreign Contamination Screening ([fcs](https://github.com/ncbi/fcs))

> [!NOTE]  
> For the different tools I used the parameters I deemed appropriate for my own analyses. To customize the pipeline towards specific needs, please inqure the documentation of the desired tool and adjust the `Snakefile` accordingly.
## Prerequsites

Apart from [Snakemake](https://snakemake.readthedocs.io/en/stable/) and [Anaconda](https://docs.anaconda.com/miniconda/), [inspector](https://github.com/Maggi-Chen/Inspector) should also be present and added to the `$PATH` variable.

For the contaminant screening the database should be present in the `bin/db/` directory. This can be done by using one of the scripts (following the [fcs-gx](https://github.com/ncbi/fcs/wiki/FCS-GX-input#fcs-gx-database-location) guide):

```
SOURCE_DB_MANIFEST="https://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/FCS/database/test-only/test-only.manifest"
LOCAL_DB="/path/to/db/folder"
python3 fcs.py db get --mft "$SOURCE_DB_MANIFEST" --dir "$LOCAL_DB/test-only" 
```

## Usage

With the complete `bin/` and all other prerequsites simply run:

```
snakemake --cores [user-defined] --use-conda all
```
