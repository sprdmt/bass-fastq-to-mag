# Manual bin refinement with anvi'o (optional)

The automatic route already yields dereplicated MAGs (`mags.source: auto`). Manual refinement
in anvi'o can rescue bins that are nearly complete but contaminated. Reviewing every bin by
hand rarely fits into a project, so the workflow triages the Binette bins first and asks for
attention only where refinement can pay off.

## Triage

For each assembly group, every Binette bin gets one decision from its CheckM2 values
(`results/anvio/<group>/triage.tsv`):

| decision | rule (defaults) | what happens |
|---|---|---|
| refine | completeness >= `anvio.refine_min_completeness` (90) and contamination > `anvio.refine_contamination_above` (5) | goes into the collection **Refine** for manual cleaning |
| keep | otherwise, passes `mags.min_completeness` / `mags.max_contamination` (50 / 10) | enters the final MAG set as it is |
| drop | everything else | not used |

Why these defaults:

* Refinement removes contigs, so completeness can only stay the same or drop. A bin below 90 %
  completeness cannot become near-complete or MIMAG high-quality by refinement.
* A bin with 95 % completeness and 7 % contamination passes the 50/10 filter as medium quality,
  but cleaned below 5 % it becomes near-complete or high-quality. Such bins are therefore offered
  for refinement too (set `refine_contamination_above: 10` to refine only bins that would
  otherwise be lost).

With `anvio.review: all`, every bin goes into **Refine** and nothing is kept automatically.

## 1. Build the workspace

```yaml
anvio:
  enabled: true
  review: contaminated          # or: all
  refine_min_completeness: 90
  refine_contamination_above: 5
```

```bash
snakemake anvio --workflow-profile profiles/slurm
```

For every assembly group this creates `results/anvio/<group>/` with

* `contigs.db`: contigs, the workflow's Prodigal gene calls, HMM hits and SCG taxonomy
* `merged/PROFILE.db`: coverage profiles of all samples of the group
* collection **Binette**: all automatic bins, a read-only reference
* collection **Refine**: the bins to review (absent if there is nothing to refine)
* `triage.tsv`: the decision and reason for every bin
* `HOW_TO_REFINE.md`: the bins to review, most complete first, each with its `anvi-refine`
  command

## 2. Refine

```bash
anvi-refine -p results/anvio/<group>/merged/PROFILE.db -c results/anvio/<group>/contigs.db \
    -C Refine -b <bin_name>
anvi-estimate-genome-completeness -p results/anvio/<group>/merged/PROFILE.db \
    -c results/anvio/<group>/contigs.db -C Refine
```

Work only in the collection **Refine**. You can split a bin, remove contaminating contigs, or
delete bins that cannot be rescued. Bins you leave untouched are used as they are.

anvi-refine serves an interactive page. On a cluster, run it in an interactive job with
`--server-only -P 8080` and open an SSH tunnel from your computer
(`ssh -L 8080:localhost:8080 <user>@<node>`). Without a local anvi'o installation, run the
workflow's image: `apptainer exec --bind $PWD <anvio image> anvi-refine ...` (the image is
listed in `workflow/containers.yaml`).

Your edits are stored in `PROFILE.db`. The workflow never replaces an existing **Refine**
collection, but a rebuilt profile (for example after changing upstream rules) starts from
scratch. So:

* back up your work: `anvi-export-collection -p .../PROFILE.db -C Refine -O Refine_backup`
* run later steps with `--rerun-triggers mtime`, so that a workflow update or a changed
  parameter does not rebuild the anvi'o workspace

## 3. Finalize

```yaml
mags:
  source: anvio
```

```bash
snakemake mags --workflow-profile profiles/slurm --rerun-triggers mtime
```

For each group the candidates are the **kept** bins (with their Binette/CheckM2 values) plus all
bins of **Refine**. The latter are exported with `anvi-summarize`, re-checked with CheckM2 and
prefixed `refined_` (`results/anvio/<group>/candidates/`). A text copy of the refined collection
is written to `results/anvio/<group>/refined_collection.txt`. The candidates of all groups are
then filtered (`mags.min_completeness`, `mags.max_contamination`), dereplicated with dRep, named
`<project>_MAG_<n>` and annotated (GTDB-Tk, rRNA/tRNA genes, coverage), exactly as in the
automatic route. The column `origin` of `results/mags/mags.tsv` tells kept MAGs (Binette origin)
and refined MAGs (`anvio_refined`) apart.

Setting `mags.source: anvio` before any refinement is valid: the **Refine** bins are then used
unchanged.

## Quality classes in the final table

`results/mags/mag_summary.tsv` reports `mimag_quality`:

| class | definition |
|---|---|
| high-quality | MIMAG: > 90 % completeness, < 5 % contamination, 23S, 16S and 5S rRNA genes and >= 18 tRNAs |
| near-complete | >= 90 % completeness and <= 5 % contamination, without the rRNA/tRNA criteria (Almeida et al. 2019) |
| medium-quality | MIMAG: >= 50 % completeness, < 10 % contamination |
| low-quality | everything else |
