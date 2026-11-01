# Reference databases

| Config key | Used by | Size (approx.) | How to obtain |
|---|---|---|---|
| `databases.checkm2` | Binette, CheckM2 | 3 GB | `checkm2 database --download` |
| `databases.gtdbtk` | GTDB-Tk (release matching the GTDB-Tk version, R226 for 2.5/2.6) | 110 GB | GTDB-Tk reference package |
| `databases.mmseqs2_taxonomy` | MMseqs2 contig taxonomy for TaxVAMB | 60-150 GB | `mmseqs databases GTDB` or `UniRef90` |
| `databases.anvio_scg` | anvi'o SCG taxonomy (optional; the workflow's anvi'o image includes it) | 1 GB | `anvi-setup-scg-taxonomy` |
| `ena.gtdb_metadata` | GTDB metadata for the NCBI lineages of the ENA tables (optional) | 0.3 GB | `bac120_metadata_r226.tsv.gz`, `ar53_metadata_r226.tsv.gz` |

All of them can be fetched with the pinned tool versions:

```bash
snakemake download_databases --cores 8 --software-deployment-method apptainer
```

The downloads run on the machine that runs Snakemake (internet access is needed) into
`databases.download_dir` (default `databases/`), and `databases/database_paths.yaml` lists the
values to copy into `config/config.yaml`. `databases.mmseqs2_download` selects the MMseqs2
database (default `GTDB`); with GPUs enabled a padded copy for MMseqs2 GPU searches is made.

### An existing MMseqs2 database and GPUs

MMseqs2 searches on GPUs need a padded copy of the target database. `makepaddedseqdb` does not copy
the taxonomy files, which `mmseqs taxonomy` needs, so link them and point
`databases.mmseqs2_taxonomy` to the padded copy:

```bash
mmseqs makepaddedseqdb db db_padded
for ext in _mapping _taxonomy; do ln -s db$ext db_padded$ext; done
```

Without a GPU (`gpu.enabled: false` or `binning.vamb.mmseqs_gpu: false`) use the database as it is.

### NCBI-based MMseqs2 databases built after March 2025

Since March 2025 the NCBI taxonomy calls the top rank of Bacteria, Archaea and Eukaryota
`domain` instead of `superkingdom` (Viruses: `acellular root`). Databases built by
`mmseqs databases` from NCBI sources (UniRef, Swiss-Prot, NR) download the current taxdump.
MMseqs2 18 rejects `domain` in `--lca-ranks`, so the first rank of every `--lca-ranks`
lineage is `unknown`. Without a correction every TaxVAMB lineage would be empty, and VAMB
stops with `NetworkXError: The node root is not in the digraph`. The step `vamb_taxonomy`
therefore reads the domain from the full lineage, the last column written with
`--tax-lineage 1`. Keep that option in `binning.vamb.mmseqs_extra`. GTDB databases built by
`mmseqs databases GTDB` still use `superkingdom` and are not affected.

Databases outside the working directory must be visible inside the containers: add their
parent directories to `--apptainer-args "--bind ..."` in your profile.
