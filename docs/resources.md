# Automatic resource estimates

Every job asks the scheduler for threads, memory (`mem_mb`) and wall time (`runtime`, minutes)
computed from **its own inputs**, so a co-assembly of 100 samples requests far more than one of
8 samples, and a mapping job scales with the size of its read files.

## The model

For each step there is an entry in [`workflow/resources.yaml`](../workflow/resources.yaml):

```
value = base
      + per_gb        * GB          # total size of the inputs listed in size_of (all inputs if omitted)
      + per_sample    * N           # samples in the assembly group
      + per_gb_sample * GB * N      # e.g. SemiBin2: contigs x samples
      + per_file      * F           # number of input files, e.g. genomes for dRep/GTDB-Tk
```

The value is raised to `min`, multiplied by `retry_factor^(attempt - 1)` for memory and run
time when a job is retried (`--retries`, set to 2 in the provided profiles), and finally
clamped to the step's `max` and to the global limits of your machine:

```yaml
resources:
  limits:
    max_threads: 96        # cores of your largest node
    max_mem_mb: 1000000    # memory of your largest node (MB)
    max_runtime: 10080     # longest wall time your queue allows (minutes)
  retry_factor: 1.5
```

Sizes are measured when a job becomes ready, i.e. on the real files. In a dry run, memory and run
time of jobs whose inputs do not exist yet are shown as `<TBD>`; thread counts fall back to `base`.

Example (illustration with round numbers): a MEGAHIT co-assembly of 45 GB of compressed reads with
the default entry `threads: {base: 16, per_gb: 1.5}`, `mem_mb: {base: 64000, per_gb: 6000}`,
`runtime: {base: 240, per_gb: 30}` requests 84 threads, 334 GB and 26.5 h. MEGAHIT receives 90%
of the allocation as an absolute byte count (`assembly.memory_fraction`), not a fraction of the
node's memory.

## Tuning

* **Cluster limits** - set `resources.limits` once per cluster (or per project).
* **One step** - override only the keys you need; the rest keep their defaults:

  ```yaml
  resources:
    rules:
      semibin2: {mem_mb: {base: 64000, per_gb: 60000}}
      gtdbtk: {mem_mb: {base: 200000}}
  ```
* **Partitions** - if your cluster does not route jobs by memory itself, list partitions and the
  smallest one that fits each job is chosen (GPU jobs use `gpu_partitions`):

  ```yaml
  resources:
    partitions:
      - {name: small, max_mem_mb: 380000, max_runtime: 4320}
      - {name: bigmem, max_mem_mb: 2000000}
  ```
* **Fixed values** - a profile's `set-resources` still overrides everything for a rule, e.g.
  `set-resources: {megahit: {mem_mb: 900000}}`.

## Calibrating from your own runs

Every rule writes a Snakemake benchmark file (`benchmarks/.../*.tsv` with `max_rss` in MB and
wall-clock `s`). After a run, compare them with the inputs to refine the coefficients:

```bash
python workflow/scripts/summarize_benchmarks.py --benchmarks benchmarks --results results \
    --out benchmark_summary.tsv
```

`max_rss` adds up the resident memory of every process of a job. For tools whose worker
processes share memory, the value can be far above the real use. In the cluster test of
v1.0.0, GTDB-Tk reported 1,190,501 MB while it ran within its 120 GB allocation. For such
jobs, check the scheduler's record instead (`sacct -j <job ID> -o JobID,MaxRSS` or
`seff <job ID>`). It reflects the actual use when the cluster accounts memory by control group.

## How the defaults were calibrated

The defaults were checked against two co-assemblies that earlier versions of this workflow
completed on a SLURM cluster:

* **medea-II:** 7 samples, 28 read files (45.3 GB compressed), a 1.19 GB co-assembly with
  427,795 contigs, a 2.33 GB Bowtie2 index and 14 per-lane BAMs (22.4 GB). All sizes are measured.
* **ecotip2:** 41 samples and 2,833,342 contigs. The assembly has since been deleted. Its size,
  about 7.9 GB, is estimated from the contig count at medea-II's 2.78 kB per contig. Steps that
  scale with reads or BAMs are therefore evaluated for medea-II only.

The Snakemake logs of those runs record the resources each job requested and when it was
submitted and finished. Peak memory was not recorded, and the scheduler accounting has expired.
Two kinds of evidence follow from these logs:

* A request under which a job completed is an upper bound on its memory need.
* The submit-to-finish time, which includes queue wait, is an upper bound on its run time.

The v1.0.0 estimate below is the model evaluated at the same input sizes. Run times have a
margin of at least 1.2× over the observed upper bounds (TaxVAMB on ecotip2), mostly 2× or more.
Memory stays at or below the old requests wherever these were generous.

| step | data set | v0.1.0 request | v0.1.0 submit to finish (min) | v1.0.0 estimate |
|---|---|---|---|---|
| `megahit` | medea-II | 84 thr, 996 GB, 4880 min | ≤ 942 | 84 thr, 336 GB, 1599 min |
| `contig2fastg` | medea-II | 1 thr, 32 GB, 30 min | ≤ 0.7 | 1 thr, 9 GB, 38 min |
| `contig2fastg` | ecotip2 | 1 thr, 32 GB, 30 min | ≤ 2.2 | 1 thr, 35 GB, 139 min |
| `fastg2gfa` | medea-II | 1 thr, 32 GB, 10 min | ≤ 0.7 | 1 thr, 9 GB, 34 min |
| `fastg2gfa` | ecotip2 | 1 thr, 32 GB, 10 min | ≤ 1.3 | 1 thr, 35 GB, 168 min |
| `pprodigal` | medea-II | 40 thr, 60 GB, 90 min | ≤ 79 | 32 thr, 18 GB, 149 min |
| `pprodigal` | ecotip2 | 40 thr, 60 GB, 90 min | ≤ 7.2 | 32 thr, 71 GB, 817 min |
| `bowtie2_index` | medea-II | 24 thr, 96 GB, 120 min | ≤ 14 | 24 thr, 30 GB, 132 min |
| `bowtie2_index` | ecotip2 | 24 thr, 96 GB, 120 min | ≤ 61 | 24 thr, 110 GB, 533 min |
| `bowtie2_map` | medea-II | 21 thr, 72 GB, 360 min | ≤ 18 | 16 thr, 30 GB, 136 min |
| `coverm_contig` | medea-II | 24 thr, 64 GB, 180 min | ≤ 3.7 | 16 thr, 27 GB, 98 min |
| `jgi_depth` | medea-II | 20 thr, 18 GB, 135 min | ≤ 10 | 16 thr, 20 GB, 98 min |
| `metabat2` | medea-II | 24 thr, 256 GB, 600 min | ≤ 8.5 | 24 thr, 52 GB, 203 min |
| `metabat2` | ecotip2 | 24 thr, 256 GB, 600 min | ≤ 444 | 24 thr, 252 GB, 1005 min |
| `genomeface` | medea-II | 2 thr, 96 GB, 420 min | ≤ 6.5 | 4 thr, 44 GB, 132 min |
| `genomeface` | ecotip2 | 2 thr, 96 GB, 420 min | ≤ 144 | 4 thr, 111 GB, 533 min |
| `semibin2` | medea-II | 12 thr, 1500 GB, 900 min | ≤ 108 | 12 thr, 80 GB, 227 min |
| `semibin2` | ecotip2 | 12 thr, 500 GB, 900 min | ≤ 238 | 12 thr, 347 GB, 798 min |
| `mmseqs2_taxonomy` | medea-II | 22 thr, 600 GB, 450 min | ≤ 90 | 24 thr, 224 GB, 263 min |
| `mmseqs2_taxonomy` | ecotip2 | 22 thr, 600 GB, 450 min | ≤ n/a | 24 thr, 357 GB, 1065 min |
| `vamb` | medea-II | 3 thr, 64 GB, 4320 min | ≤ 284 | 4 thr, 26 GB, 596 min |
| `vamb` | ecotip2 | 3 thr, 64 GB, 4320 min | ≤ 2625 | 4 thr, 79 GB, 3268 min |
| `graphbin2` | medea-II | 4 thr, 64 GB, 60 min | ≤ 2.5 | 1 thr, 6 GB, 32 min |
| `graphbin2` | ecotip2 | 4 thr, 64 GB, 60 min | ≤ 25 | 1 thr, 20 GB, 99 min |
| `binette` | medea-II | 16 thr, 192 GB, 500 min | ≤ 64 | 16 thr, 56 GB, 303 min |
| `binette` | ecotip2 | 16 thr, 192 GB, 500 min | ≤ 200 | 16 thr, 189 GB, 1131 min |
| `anvio_contigs_db` | medea-II | 12 thr, 96 GB, 240 min | ≤ 32 | 12 thr, 30 GB, 203 min |
| `anvio_contigs_db` | ecotip2 | 12 thr, 96 GB, 240 min | ≤ 244 | 12 thr, 110 GB, 1005 min |
| `anvio_hmms` | medea-II | 24 thr, 128 GB, 600 min | ≤ 67 | 24 thr, 39 GB, 198 min |
| `anvio_hmms` | ecotip2 | 24 thr, 128 GB, 600 min | ≤ 346 | 24 thr, 168 GB, 974 min |
| `anvio_profile` | medea-II | 6 thr, 256 GB, 900 min | ≤ 440 | 6 thr, 58 GB, 1022 min |
| `anvio_merge` | medea-II | 8 thr, 128 GB, 120 min | ≤ 11 | 8 thr, 46 GB, 165 min |
| `anvio_merge` | ecotip2 | 8 thr, 128 GB, 120 min | ≤ 117 | 8 thr, 114 GB, 675 min |

Notes on the comparison:

* v0.1.0 treated every FASTQ pair as a sample, so medea-II was mapped and profiled per lane
  (14 BAMs of about 1.6 GB). With the two lanes of each sample concatenated into one file pair,
  as `config/README.md` recommends, v1.0.0 maps 7 BAMs of about 3.2 GB, and `bowtie2_map` and
  `anvio_profile` jobs are about twice as large as the old jobs.
* `graphbin2` was measured directly with `graphbin2_fast.py` on medea-II: 21-29 s and at most
  1.6 GB per binner ([graphbin2.md](graphbin2.md)). The v0.1.0 times are those of the previous
  patched GraphBin2.
* Tool logs give pure run times for some steps: SemiBin2 took 50 min on medea-II (GPU). GTDB-Tk
  took 57-156 min for 80-371 genomes, and CheckM2 3-83 min. The defaults request
  `120 + 2 × genomes` and `30 + genomes` minutes.
* `mmseqs2_taxonomy` used a UniRef90 database in v0.1.0 (600 GB requested, ≤ 90 min with a GPU).
  The v1.0.0 base of 200 GB assumes a database of that size. A small database needs far less.
  The cluster test with UniProtKB/Swiss-Prot completed with
  `resources: {rules: {mmseqs2_taxonomy: {mem_mb: {base: 48000}}}}`.
* Retries are not cause-aware: a job that fails for any reason is resubmitted with 1.5× the
  memory and run time, which can move it to a larger node and a longer queue. Fix the cause
  before resubmitting jobs that fail for reasons other than memory or time.

These figures are upper bounds and estimates, not measured peaks. After your first large run,
`summarize_benchmarks.py` gives the measured peak memory and run time of every job, and the
coefficients can be tightened from those.
