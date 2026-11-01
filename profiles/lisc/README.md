# LiSC (University of Vienna) notes

Facts checked on 2026-10-31 for running this workflow on the Life Science Compute Cluster.

## Databases

| purpose | path |
|---|---|
| GTDB-Tk reference data (GTDB R226; `gtdbtk check_install` passed with GTDB-Tk 2.6.1) | `/lisc/opt/mirror/GTDBTk/r226` (also `r220`, `r232`) |
| CheckM2 DIAMOND database (CheckM2, Binette) | `/lisc/opt/sw/software/CheckM2/1.1.0-foss-2024a/CheckM2_database/uniref100.KO.1.dmnd` |
| anvi'o SCG taxonomy | built into the workflow's anvi'o image (leave `databases.anvio_scg` empty) |
| MMseqs2 taxonomy database for TaxVAMB | not provided centrally; build one with `snakemake download_databases` (`databases.mmseqs2_download`) |
| NCBI taxdump (ENA lineages) | `/lisc/opt/mirror/ncbi/current/taxonomy/taxdump.tar.gz` |

## GPUs and Slurm

* The GRES type is `gpu` (there is no `shard` GRES any more): the default `gpu.gres: "gpu:1"` works.
* A job-submit plugin rewrites partition requests, so the GPU model cannot be chosen with
  `--partition`. Test jobs ran on Tesla T4 GPUs (16 GB); GenomeFace (CUDA 11.7 user space) and
  PyTorch 2.8 (CUDA 12.9) both ran there.
* GPU nodes (scontrol): `zen2_1024_t4x1` (16 x 1 T4), `xeon_0768_t4x1`, `xeon_1024_t4x1`,
  `zen1_0512_t4x1`, `zen2_0512_l40sx2`, `zen4_0768_rtxx2`, `zen4_0768_a30x4`, `zen5_1536_h200x2`,
  `zen3_4096_l4x3`, `zen3_2048_l4x3`.
* `sacct -j <id>` reports jobs of the same day (memory and GPU use), useful to check the resource
  estimates.

## Apptainer

* Apptainer 1.5.5 (non-setuid). `apptainer build --fakeroot` works on compute nodes, so the
  images of this repository can be built on LiSC from `containers/<image>/<image>.def`
  (add `--notest` when the `%test` section needs a GPU).
* Compute nodes reach quay.io, ghcr.io, Docker Hub, community.wave.seqera.io and the conda channels.
* LiSC modules export `PYTHONPATH`; `profiles/lisc` keeps it out of the containers
  (`--env PYTHONPATH= --env PYTHONNOUSERSITE=1`).
* The login environment exports bash functions (`which`, `module`, `ml` and other Lmod
  helpers), and Apptainer passes them into the containers. The exported `which` broke barrnap's
  check for `nhmmer`; the workflow unsets it at the start of every shell command
  (`shell.prefix` in `workflow/Snakefile`).
* Set `apptainer-prefix` in the profile to a shared directory to keep one copy of the pulled
  images (about 30 GB) for all projects.
