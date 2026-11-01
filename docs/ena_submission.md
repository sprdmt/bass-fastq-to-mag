# Submitting assemblies and genomes to ENA

This module turns the workflow's final MAG set into ENA submissions with the
EBI tools **genome_uploader** (sample registration and manifests for bins/MAGs)
and **webin_cli_handler** / **Webin-CLI** (file upload), plus a small helper
for the primary co-assembly. Every command uses the **ENA TEST server** unless
`--live` is given.

## 1. ENA assembly tiers

| Tier | ENA object | Sample it points to | Produced here |
|---|---|---|---|
| Primary metagenome | one assembly per assembly group (`ASSEMBLY_TYPE primary metagenome`) | an **environmental** sample; for a co-assembly a new **co-assembly sample** whose `sample derived from` lists all source samples | `ena_coassembly.py` (sample XML, manifest) + `ena_submit.py coassembly` |
| Binned metagenome | one assembly per bin (`binned metagenome`), checklist ERC000050 | one **bin sample** per genome, `sample derived from` = source sample(s) | `ena_genome_table.py` + `ena_submit.py genomes --genome-type bins` (default) |
| MAG | one assembly per MAG (`Metagenome-Assembled Genome (MAG)`), checklist ERC000047 (GSC MIMAGS) | one **MAG sample** per genome | same with `--genome-type mags` |

* Binned metagenomes and MAGs are *derived* assemblies; ENA wants a separate
  organism-level sample for each, linked to the environmental samples through
  `sample derived from`. genome_uploader creates these samples from the run
  accessions you give it, copying collection date and location from the raw-read
  samples (all distinct sample accessions of a co-assembly are listed).
* ENA asks that only the highest-quality, **unique-taxon** genomes are
  submitted as MAGs (dereplicated to one representative per predicted species).
  The workflow's dRep step clusters at 99.5 % ANI (strain level), so the default
  tier is **bins**. Submit a species-level representative subset as MAGs if
  wanted (`ena_genome_table.py --mag-ids FILE --genome-type mags`).
* The primary tier is optional for the genome tiers: genome_uploader links bins
  and MAGs to the **runs** (not to the primary assembly). Passing the primary
  assembly accession (ERZ) instead of runs is possible in genome_uploader, but it
  inherits metadata from the **first** RUN_REF only, so the workflow always
  passes all runs of the co-assembly.

### Co-assemblies (primary tier)

Webin-CLI accepts exactly one `SAMPLE` per assembly. The earlier attempt with a
comma-separated list of the 24 source samples failed with
*"Unknown sample ERS…,ERS… or the sample cannot be referenced by your
submission account"*. ENA's FAQ "How do I register samples for co-assemblies?"
(read_docs `faq/assemblies-analyses.rst`, text unchanged since 2019) asks to

* reference **all source samples or runs** in the `sample derived from` field of
  the sample used for the assembly, as a comma-separated list **without spaces**
  (`ERSxxxxxx,ERSxxxxxx`) or as a range `ERSxxxxxx-ERSxxxxxx` that contains only
  accessions used in the assembly, with one accession format throughout;
* **inform the ENA helpdesk** before submitting a primary assembly that was
  co-assembled from raw reads.

`ena_coassembly.py sample-xml` builds that sample: environmental taxon
(`ena.metagenome_taxid`, e.g. 1561972 *seawater metagenome*), checklist
ERC000011, `sample derived from` = the group's `sample_accession` values (run
accessions if those are missing), and the mandatory `collection date` and
`geographic location (country and/or sea)`. A mandatory field without a single
true value gets the INSDC term `missing: sample group` ("the sample represents
a group of samples that do not have a single origin, e.g. for co-assembly").
`ena_submit.py coassembly --live` refuses to run without `--helpdesk-notified`.

## 2. Prerequisites

1. A Webin account; credentials only as environment variables
   (`export ENA_WEBIN=Webin-XXXXX`, `export ENA_WEBIN_PASSWORD='...'` in single
   quotes if the password has special characters). The scripts never read or
   write credential files; genome_uploader itself also loads
   `~/.genome_uploader.config.env` if that file exists (the wrapper runs it in
   the work directory, so a project `.env` is not picked up), and values already
   in the environment take precedence.
2. The **study** (PRJEB…) and the **raw reads** (runs ERR…, environmental
   samples ERS…) are already submitted by the same Webin account. Private data
   are fine (`ena.private: true` → genome_uploader `--private`, which reads the
   metadata through the authenticated Webin report service).
3. `config/samples.tsv` carries the accessions:

   | sample | fq1 | fq2 | assembly_group | run_accession | sample_accession |
   |---|---|---|---|---|---|
   | S1 | … | … | g1 | ERR0000001 | ERS0000001 |

   Several runs of one sample: list them comma-separated in `run_accession`.
4. For the default lineage method: the GTDB metadata files of the GTDB-Tk
   reference release (R226 for GTDB-Tk 2.5: `bac120_metadata_r226.tsv.gz`,
   `ar53_metadata_r226.tsv.gz` from data.gtdb.ecogenomic.org) and the GTDB-Tk
   `classify_wf` output directory (with `classify/*.tree`). Keep the metadata
   next to the GTDB-Tk reference data; when the workflow moves to a newer GTDB
   release (the GTDB R232 files are dated April 2026) the metadata and workbooks must be
   from that same release.

## 3. Configuration

Settings are in the `ena:` block of `config/config.yaml` (all keys and defaults:
`workflow/schemas/config.schema.yaml`). Key points:

| key | default | meaning |
|---|---|---|
| `study`, `centre_name` | — | genome_upload `-u`, `--centre_name` |
| `genome_type` | `bins` | ENA tier for the genome table |
| `private` | `true` | raw data still private |
| `metagenome`, `metagenome_taxid` | seawater metagenome / 1561972 | genome_uploader `metagenome` column (must be in its accepted list) and co-assembly sample taxon |
| `broad_environment`, `local_environment`, `environmental_medium` | — | ENVO terms; `groups:` overrides per assembly group |
| `assembly_software`, `binning_software`, `binning_parameters`, `stats_generation_software` | — | genome_uploader columns (`toolName_vX.X`) |
| `lineage_method`, `lineage_cross_check`, `lineage_max_rank` | `majority_vote`, `false`, `genus` | see §5 (`mapping_table` / `consensus` / cross-check need `gtdb_ncbi_workbooks` and `taxdump`) |
| `mags.coverage_method` | `sum_mean` | genome_coverage, see §5 (also used in `mag_summary.tsv`) |
| `mags.rna.*` | partial rRNA hits counted, `standard20`, 18 tRNAs | rRNA/tRNA criteria, see §5 (also used for the quality class) |
| `assembly_software`, `binning_software`, `stats_generation_software`, `coassembly.program` | empty | filled from the versions of the workflow's images; for local image files (`.sif`) the versions of the default images are used, so set these fields if your images differ |
| `coassembly.*` | — | alias prefix, checklist, shared attributes, PROGRAM/PLATFORM |
| `webin.java_heap_max_gb` | 4 | webin_cli_handler otherwise asks Java for 10 GB |

The first ENA job of a run (`ena_params`) checks that `metagenome`, the three environment
fields, the software fields and `binning_parameters` are set for every assembly group. If one is
empty, the run stops at its start with the names of the missing keys, not after binning and
GTDB-Tk.

There is deliberately **no config switch for live submission**.

## 4. Commands

### Workflow targets

```bash
snakemake ena_prepare  --workflow-profile profiles/slurm   # tables, FASTAs, lineages, co-assembly sample XML
snakemake ena_validate --workflow-profile profiles/slurm   # genome_uploader validator + ENA taxon preview
export ENA_WEBIN=Webin-XXXXX ENA_WEBIN_PASSWORD='...'
snakemake ena_submit   --workflow-profile profiles/slurm   # ENA TEST server only
```

`ena_prepare` writes `results/ena/genome_uploader_input.tsv` (+ `.provenance.tsv`),
`results/ena/fasta/`, `results/ena/ncbi_lineage.tsv` and `results/ena/coassembly/<group>/sample/`
for every assembly group with at least two samples. `ena_validate` and `ena_submit` run on the
machine that runs Snakemake (they need internet access). `ena_submit` registers the bin samples,
validates and submits all genomes and the co-assemblies to the test server and writes
`results/ena/submit/` (state, logs, test accessions). Check the test submission, then submit live
with the commands below (`--live`), which are deliberately not workflow targets.

A group with a single sample needs no co-assembly sample: its primary assembly references the
existing sample. Build the manifest with
`ena_coassembly.py manifest --sample-accession ERS... --assembly-fasta results/assembly/<group>/contigs.fa ...`
and submit it with Webin-CLI.

### Scripts


Workflow outputs used below: `results/mags/mags.tsv`, `results/mags/fasta/`,
`results/mags/rna/{barrnap,trnascan}/`, `results/mags/gtdbtk/`,
`results/mapping/<group>/coverm_contig_mean.tsv`, `config/samples.tsv`.

```bash
S=workflow/scripts
# tables (python container; Python >= 3.9, standard library only)
python $S/mimag_rna.py --barrnap-dir results/mags/rna/barrnap --trnascan-dir results/mags/rna/trnascan \
    --mags-table results/mags/mags.tsv --gtdbtk results/mags/gtdbtk/gtdbtk.summary.tsv --out results/mags/mimag_rna.tsv
python $S/mag_coverage.py --mags-table results/mags/mags.tsv --mags-dir results/mags/fasta \
    --coverm g1=results/mapping/g1/coverm_contig_mean.tsv --out results/mags/coverage.tsv
# GTDB-Tk container: GTDB-Tk's own majority vote
gtdb_to_ncbi_majority_vote.py --gtdbtk_output_dir results/mags/gtdbtk \
    --bac120_metadata_file bac120_metadata_r226.tsv.gz --ar53_metadata_file ar53_metadata_r226.tsv.gz \
    --output_file results/ena/gtdbtk_majority_vote.tsv
python $S/gtdb_to_ncbi.py --gtdbtk results/mags/gtdbtk/gtdbtk.summary.tsv \
    --mv-table results/ena/gtdbtk_majority_vote.tsv --out results/ena/ncbi_lineage.tsv
#   optional review of the majority vote against the GTDB workbooks (adds column cross_check)
python $S/gtdb_to_ncbi.py --gtdbtk results/mags/gtdbtk/gtdbtk.summary.tsv \
    --mv-table results/ena/gtdbtk_majority_vote.tsv --cross-check \
    --mapping-table gtdb_vs_ncbi_r226_bacteria.xlsx gtdb_vs_ncbi_r226_archaea.xlsx \
    --taxdump taxdmp_2025-12-01.zip --out results/ena/ncbi_lineage.tsv
python $S/ena_genome_table.py --mags-table results/mags/mags.tsv --mags-dir results/mags/fasta \
    --rna results/mags/mimag_rna.tsv --coverage results/mags/coverage.tsv --lineage results/ena/ncbi_lineage.tsv \
    --samples config/samples.tsv --params config/config.yaml \
    --fasta-outdir results/ena/fasta --out results/ena/genome_uploader_input.tsv
python $S/ena_coassembly.py sample-xml --samples config/samples.tsv --group g1 --params config/config.yaml \
    --outdir results/ena/coassembly/g1

# submission (genome-uploader container), always TEST first
export ENA_WEBIN=Webin-XXXXX ENA_WEBIN_PASSWORD='...'
python $S/ena_submit.py validate-table --genome-info results/ena/genome_uploader_input.tsv --genome-type bins
python $S/ena_submit.py preview-taxa --genome-info results/ena/genome_uploader_input.tsv --out results/ena/taxa_preview.tsv
python $S/ena_submit.py genomes --genome-info results/ena/genome_uploader_input.tsv --genome-type bins \
    --study PRJEB12345 --centre-name "University of Vienna" --workdir results/ena/submit          # TEST server
python $S/ena_submit.py genomes ... --workdir results/ena/submit --live                            # PRODUCTION
python $S/ena_submit.py coassembly --xml-dir results/ena/coassembly/g1 --samples config/samples.tsv --group g1 \
    --study PRJEB12345 --assembly-fasta results/assembly/g1/contigs.fa --assembly-name g1_coassembly \
    --program "MEGAHIT v1.2.9" --coverm results/mapping/g1/coverm_contig_mean.tsv \
    --outdir results/ena/coassembly/g1 --workdir results/ena/submit                               # TEST server
python $S/ena_submit.py coassembly ... --live --helpdesk-notified                                  # PRODUCTION
python $S/ena_submit.py status --workdir results/ena/submit
```

`genomes` runs `genome_upload` (sample registration, manifests), then for every
manifest `webin_cli_handler --mode validate` and `--mode submit` (`--test` in
test mode), each in its own directory `<workdir>/webin_<mode>/<alias>/`. Add
`--dry-run` to print the commands without credentials, `--stop-after register|validate`
to pause.

### Safety model

* TEST is the default everywhere; the test server is rebuilt daily (test objects
  disappear the next day; test aliases carry a timestamp suffix kept in the state file).
* `--live` requires the credentials in the environment **and** a completed TEST
  run of the identical table (MD5) or manifest in the same `--workdir`
  (`--skip-test-gate` overrides), and `--helpdesk-notified` for a co-assembly.
* Subprocesses run with the work directory as cwd (no project `.env` is
  picked up); all logs are written with the password redacted; webin-cli gets the
  password through `-passwordEnv`.
* Re-running a command resumes: registered samples, validated manifests (by MD5)
  and submitted assemblies are skipped. Accessions are in
  `<workdir>/ena_accessions_<test|live>.tsv` (genome_name, alias, sample_accession,
  analysis_accession) and the co-assembly accession TSV.

## 5. Definitions used in the tables

* **genome_coverage** (`mag_coverage.py`, default `sum_mean`): for each contig
  the CoverM `Mean` depths of all samples of the MAG's group are summed; the
  genome value is the mean of these sums over the MAG's contigs (each contig
  counts once). This is the definition used for the earlier submissions.
  `length_weighted` reports Σ_s Σ_c L_c·d_cs / Σ_c L_c instead. The output also has
  `mean_depth_<sample>` (length-weighted per sample). The primary-assembly
  `COVERAGE` uses the same definition over all contigs of the table.
* **rRNA_presence** (genome_uploader): True if barrnap finds 5S, 16S and 23S rRNA
  **and** at least 18 distinct standard tRNA isotypes are found. Partial barrnap
  hits count (`--exclude-partial-rrna` to change). tRNAs: tRNAscan-SE, non-pseudo,
  20 standard amino acids (iMet/fMet → Met, Ile2 → Ile; SeC, Pyl, Sup, Undet ignored).
* **mimag_quality**: `high-quality` (completeness > 90, contamination < 5, 23S+16S+5S,
  ≥ 18 tRNAs; MIMAG), `near-complete` (≥ 90 / ≤ 5 but not MIMAG high-quality;
  Almeida et al. 2019), `medium-quality` (≥ 50 / < 10), `low-quality` (rest).
* **ENA "assembly quality"** text is written by genome_uploader itself: the MIMAG
  high-quality text if completeness ≥ 90, contamination ≤ 5 and rRNA_presence,
  otherwise the medium-quality text. Because MIMAG requires > 90 and < 5,
  `ena_genome_table.py` stops on genomes at exactly 90.00 % or 5.00 % with
  rRNA_presence True (`--allow-assembly-quality-conflict` to accept them).
* **NCBI_lineage** (`gtdb_to_ncbi.py`): `d__…;p__…;c__…;o__…;f__…;g__…;s__` with NCBI
  names. genome_uploader needs the `x__` prefixes (plain names raise
  *Unrecognised taxonomy format*) and picks the first submittable ENA name from
  the most specific rank upwards (`Genus sp.`, `Family bacterium`, …).
  * `majority_vote` (default): GTDB-Tk's `gtdb_to_ncbi_majority_vote.py`; per
    rank, the NCBI taxon of ≥ 50 % of the GTDB reference genomes around the
    genome's tree placement (same GTDB release as the reference data).
  * `mapping_table` (previous route, kept as option): best NCBI taxon ≥ 50 % from
    the GTDB `gtdb_vs_ncbi` workbooks for the most specific mapped GTDB rank,
    else a direct name match in the NCBI taxdump, else the domain.
  * `--cross-check` (with `majority_vote`): also runs `mapping_table` and writes
    every rank where the two disagree to the column `cross_check` (a
    'Candidatus' prefix alone is not a disagreement); the lineage stays the
    majority vote. Review these genomes before submission.
  * `consensus`: keeps the ranks, from the domain down, on which both methods
    agree, in the current NCBI spelling; everything from the first disagreement
    on is left empty.
  * On a 172-MAG test set (marine co-assembly, GTDB R226) the two methods give the same lineage for
    97 genomes (one of them differs only in the spelling 'Candidatus
    Thermoplasmatota'), the majority vote is shallower for 59 and deeper for 4,
    and they conflict for 12 (first conflicting rank: phylum 2, class 2, order 2,
    family 3, genus 3). All 12 conflicts are genomes for which `mapping_table`
    used a direct GTDB-to-NCBI name match rather than a workbook mapping (134 of
    172 genomes); `consensus` is shallower than `majority_vote` for 17 genomes.
  * Lineages are cut after the genus (`--max-rank genus`), following ENA's
    guidance to name organism-level environmental taxa "up to genus level"; this
    keeps strain-specific NCBI species names (e.g. `X sp. ABC123`) off the samples.
  * `preview-taxa` shows the taxon genome_uploader will register (read-only ENA
    taxonomy service) before anything is registered.

## 6. Troubleshooting

| message | cause / fix |
|---|---|
| `Unknown sample ERS…,ERS… or the sample cannot be referenced` | several samples in `SAMPLE`; register the co-assembly sample (`ena_coassembly.py`) and use its accession |
| `Credentials for ENA are not set` / exit code 2 | export `ENA_WEBIN` and `ENA_WEBIN_PASSWORD` in the shell that runs the command |
| `Invalid submission account user name or password` | wrong credentials, or TEST/LIVE confusion |
| `refusing --live: no completed TEST submission` | run the same command without `--live` first (same `--workdir`, unchanged table) |
| `Unrecognised taxonomy format` (genome_uploader) | NCBI_lineage without `x__` prefixes; build it with `gtdb_to_ncbi.py` |
| `Metagenomes associated with each genome need to belong to ENA's approved metagenomes list` | `ena.metagenome` not in genome_uploader's list (e.g. "sea surface microlayer metagenome" is not; use "seawater metagenome" or "marine metagenome") |
| `discrepancy between number of runs involved and co-assembly status` | edit samples.tsv; the table builder derives co-assembly from the number of runs |
| genome name too long | genome_uploader appends `_<10-digit timestamp>` in test mode; keep names ≤ 39 characters |
| `missing: third party data` in genome_samples.xml (test run) | genome_uploader did not find country/latitude/longitude/date in the raw-read samples or the country is not in its INSDC country list; check the raw-read sample attributes before `--live` |
| `The object being added already exists` | already registered/submitted; the wrapper and webin_cli_handler reuse the accession |
| webin-cli reports its version is not supported | update `ena-webin-cli` in `workflow/envs/ena.yaml` |
| `java.lang.OutOfMemoryError` | raise `--java-heap-max-gb` |
| FASTA rule errors (`terminal N`, `< 20 bp`, invalid characters) | fix the FASTA; webin-cli applies the same rules |
| `MAG contigs missing from the CoverM table` | wrong group table or renamed contigs (e.g. anvi'o-renamed contigs vs MEGAHIT names) |
| MAG tier: `fails MAG criteria` | completeness < 50 or contamination ≥ 10: submit as bins or drop |

## 7. Differences to expect against a manual (Webin portal) preparation

Checked on a 172-MAG co-assembly set that had been prepared by hand before
(bin samples registered through the Webin portal spreadsheet, manifests written
by a script, co-assembly manifest with all source samples in `SAMPLE`):

| item | manual preparation | this module (genome_uploader) |
|---|---|---|
| bin sample registration | portal spreadsheet, alias `ena-SAMPLE-TAB-…` | genome_uploader, alias = `genome_name` (test: `_<timestamp>` added) |
| re-use of existing bin samples | — | only by identical alias in the same Webin account; portal-registered samples are **not** re-used, a new run registers new samples |
| organism name | `uncultured <taxon> bacterium` | genome_uploader's choice: `<Genus> sp.`, `<taxon> bacterium`, … (ENA recommendation) |
| HQ "assembly quality" text | ≥ 18 tRNA types incl. fMet/SeC as separate types | ≥ 18 of the 20 standard amino acids (fMet → Met; SeC, Pyl, Sup ignored) |
| `assembly software` | free text incl. parameters | `ena.assembly_software` (`MEGAHIT_v1.2.9`; parameters may be added) |
| `completeness software` | `CheckM2` | `CheckM2_v1.1.0` |
| `taxonomic identity marker` | `GTDB-Tk;…;GTDB release 226` | `multi-marker approach` (fixed by genome_uploader) |
| `isolation_source` | free text | the `metagenome` value (fixed by genome_uploader) |
| collection date, country, latitude, longitude | from the config, one value per data set | copied by genome_uploader from the raw-read samples; for a co-assembly `not provided` (date, coordinates) or `not applicable` (country) when the source samples differ |
| `sample derived from` | empty | the raw-read sample accessions |
| manifest `PROGRAM` / `PLATFORM` | binning software / `ILLUMINA` | `assembly_software` / instrument model of the runs (e.g. `Illumina NovaSeq 6000`) |
| manifest `MINGAPLENGTH` | 1 | not written |
| co-assembly `SAMPLE` | 24 comma-separated accessions (rejected by webin-cli: *Unknown sample …*) | one co-assembly sample with `sample derived from` (`ena_coassembly.py`) |

The TEST run writes `genome_samples.xml` and the manifests to the work
directory; check these (in particular date, location and taxon) before `--live`.
