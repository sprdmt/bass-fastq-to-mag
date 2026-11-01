# Changelog

All notable changes to this workflow are documented here. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and the project uses
[semantic versioning](https://semver.org/).

## [1.0.0] - 2026-11-01

A rewrite that makes the workflow portable, reproducible and public. Version 0.1.0 ran only on
the LiSC cluster of the University of Vienna (environment modules, personal conda environments,
hard-coded database paths).

### Added

* Standard Snakemake layout (`workflow/`, `config/`, `profiles/`) with a JSON-schema-validated
  configuration and sample sheet.
* Containers for every rule (BioContainers; images of this repository for GenomeFace, SemiBin2
  and VAMB with GPU support, anvi'o 9 and the helper scripts) and pinned conda environments as an
  alternative. Shell commands first unset an exported bash function `which`: clusters with RHEL's
  `which2.sh` or Lmod export it, Apptainer passes it into the containers, and there it broke
  barrnap's check for `nhmmer`.
* Read QC: cutadapt trimming with the original BASS parameters, FastQC before and after
  trimming, MultiQC report.
* Several co-assembly groups per run (`assembly_group` column), dereplicated together.
* Size-aware resource estimates for every rule, capped by user-defined cluster limits, with
  memory and run-time escalation on retries (`workflow/resources.yaml`, `docs/resources.md`).
  The defaults were checked against the requests and run times of the medea-II and ecotip2 runs
  of v0.1.0 (table in `docs/resources.md`).
* SLURM, LiSC and local profiles; GPU requests via `--gres` or `--gpus`.
* Automatic triage of the Binette bins for manual refinement in anvi'o (keep / refine / drop,
  collection `Refine`, `HOW_TO_REFINE.md` with one command per bin). Final candidates combine
  kept and refined bins (`mags.source: anvio`).
* Final MAG naming (`<project>_MAG_<n>`), GTDB-Tk classification, barrnap and tRNAscan-SE for
  the MIMAG rRNA/tRNA criteria, coverage per sample and one summary table.
* Quality classes MIMAG high-quality, near-complete (>= 90 % / <= 5 % without the rRNA/tRNA
  criteria), medium-quality and low-quality.
* ENA submission of the co-assemblies and MAGs with genome_uploader and Webin-CLI; the test
  server is used unless live submission is requested explicitly. A co-assembly references one
  registered co-assembly sample (the earlier separate scripts listed all source samples, which
  ENA rejected), and the first ENA job checks that the required `ena:` fields are set.
* `download_databases` target for CheckM2, GTDB-Tk (R226), MMseqs2 and anvi'o data.
* Continuous integration: snakefmt, lint, dry runs of all targets and configuration variants,
  unit tests of the helper and ENA scripts, and a GraphBin2 test suite that compares
  `graphbin2_fast.py` with upstream GraphBin2 1.3.3 carrying only the two bug fixes
  (`tests/graphbin2_fast/`).

### Changed

* GraphBin2 runs a corrected and faster implementation (`workflow/scripts/graphbin2_fast.py`)
  that keeps the label propagation of GraphBin2 1.3.3 (see `docs/graphbin2.md`). On the medea-II
  co-assembly (427,795 contigs) it needs 21-29 s and <= 1.6 GB per binner, compared with 3.7-4.9 h
  for upstream with the bugs fixed. By default Binette receives both the original and the
  GraphBin2 bins of every binner (`graphbin2.binette_input`).
* Tool versions (containers/README.md): SemiBin2 2.2.0 -> 2.4.1, Binette 1.1.2 -> 1.2.1 (final bins
  are now named `binette_bin<N>`), GTDB-Tk 2.5.2 -> 2.6.1 (still GTDB R226); MEGAHIT 1.2.9,
  VAMB 5.0.4, CheckM2 1.1.0, dRep 3.6.2 and genome_uploader 2.5.1 as before; MetaBAT2 2.18.
* GPUs are requested as `--gres=gpu:1` (the `shard` GRES of the old LiSC profile no longer exists).
* Containers run with `--env PYTHONPATH= --env PYTHONNOUSERSITE=1`, so a `PYTHONPATH` exported by
  environment modules (replacing the `module purge; unset PYTHONPATH` lines of v0.1.0) or
  `~/.local` packages cannot leak into them.
* Samples come from a sample sheet (`config/samples.tsv`). v0.1.0 took every
  `<name>_1.fastq.gz`/`<name>_2.fastq.gz` pair in `fastq/` as a sample, so each sequencing lane
  became a sample of its own and all samples formed one co-assembly.
* Prodigal runs in metagenome mode (`-p meta`) by default; v0.1.0 used single mode.
* MetaBAT2 runs with a fixed seed (`binning.metabat2.seed: 1`). v0.1.0 used MetaBAT2's default
  seed 0, which draws a random seed, so its bins could differ between runs.
* MAG completeness and contamination come from Binette's CheckM2 run and are passed to dRep
  (`--genomeInfo`); v0.1.0 ran CheckM2 again on the dRep output. A separate CheckM2 run remains
  for bins refined in anvi'o.
* anvi'o profiles and merges use the documented `-W` (`--overwrite-output-destinations`)
  instead of anvi'o's global `--force-overwrite` switch.
* The MEGAHIT memory limit is given in bytes from the job's allocation (v0.1.0 used a fraction of
  the node's memory, which can exceed the SLURM allocation).
* MAG filtering: Binette keeps bins >= 70 % complete (CheckM2) and dRep clusters at 99.5 % ANI,
  both as before. dRep now filters with Binette's CheckM2 values at >= 50 % completeness and
  <= 10 % contamination (MIMAG medium quality). v0.1.0 let dRep run CheckM (v1) and filter at
  dRep's defaults, >= 75 % completeness and <= 25 % contamination.

### Fixed

* The assembly graph was built with k = 127 while MEGAHIT's final k was 125 (`--k-max 125`).
  The graph had no edges and GraphBin2 changed no bins. The k is now read
  from the contig names, and the workflow stops when the graph has no links
  (`graph.require_edges`).
* GraphBin2 1.3.3: coverage looked up with the wrong (1-based) index, and an undefined variable
  (`contig_names`) on the output path.
* The GraphBin2 patch of v0.1.0 stopped label propagation early; the new implementation keeps
  the upstream algorithm (see `docs/graphbin2.md`).
* TaxVAMB with an NCBI-based MMseqs2 database built after March 2025: NCBI renamed the rank
  `superkingdom` to `domain`, which MMseqs2 18 does not accept in `--lca-ranks`. The lineage
  clean-up (the inline step `fix_taxonomy_for_vamb` of v0.1.0) then emptied every lineage and VAMB
  failed. The domain is now taken from the MMseqs2 full lineage (`docs/databases.md`).
* The CheckM2 database path of v0.1.0 (`/lisc/data/scratch/mirror/checkm2/...`) no longer exists;
  database paths are now configuration values.
