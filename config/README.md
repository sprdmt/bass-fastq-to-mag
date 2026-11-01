# Configuration

Two files describe a project:

* **`config/samples.tsv`** - one row per sample (tab-separated). Each row is one coverage
  profile for the binners, so give the lanes of one sample as a single file pair: `cat` the
  lane files (gzipped files can be concatenated) into one R1 and one R2 file.

  | column | required | meaning |
  |---|---|---|
  | `sample` | yes | unique name (letters, digits, `.`, `_`, `-`) |
  | `fq1`, `fq2` | yes | paths to the read files (gzipped or not) |
  | `assembly_group` | no | samples with the same group are co-assembled and co-binned; empty = one group named after `project` |
  | `run_accession` | for ENA | ENA/SRA run accession of the reads |
  | `sample_accession` | for ENA | ENA/BioSample accession of the sample |

* **`config/config.yaml`** - all settings with comments; every key and its default is defined in
  `workflow/schemas/config.schema.yaml`, so a minimal config only needs `project`, the database
  paths and anything you want to change.

Settings you will usually adapt:

* `qc.*` - adapters and trimming (defaults: TruSeq adapters, 9 bases clipped at the 5' ends,
  two-colour poly-G trimming), or `qc.enabled: false` for reads that are already trimmed
* `binning.binners` - any subset of `metabat2`, `genomeface`, `semibin2`, `vamb`
* `gpu.*` - `enabled: false` runs SemiBin2, TaxVAMB and MMseqs2 on CPUs. GenomeFace needs a GPU,
  so remove `genomeface` from `binning.binners` as well (the workflow stops at the start
  otherwise). `gres` is the SLURM request for one GPU on your cluster
* `databases.*` - CheckM2, GTDB-Tk and MMseqs2 databases (see `docs/databases.md`)
* `resources.limits` - largest node of your cluster; the per-job estimates never exceed it
  (see `docs/resources.md`)
* `mags.source` - `auto`, or `anvio` after manual refinement (see `docs/manual_refinement.md`)
* `ena.*` - submission metadata (see `docs/ena_submission.md`)
