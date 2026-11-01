# GraphBin2 (MEGAHIT mode): corrected fast implementation and validation on medea-II

Track report for the bass-fastq-to-mag v1.0.0 restructuring, 2026-10-31. All numbers below come from runs on LiSC or local test runs; the comparison tables and job scripts named in Section 8 were delivered with the validation and are not part of this repository.

## Summary

- **Zero edges came from a k mismatch.** The medea-II co-assembly has final k = 125 (all 427,795 headers `k125_`), but the graph was built with k = 127. Rebuilt with k = 125 it has 48,310 edges; the k = 127 rebuild has none and is byte-identical to the GFA the old workflow used. GraphBin2 therefore left the bins of all four binners unchanged in the earlier runs.
- **Upstream GraphBin2 1.3.3 (MEGAHIT mode) cannot process this data.** On the k = 125 graph, the unmodified code ran 2.9–4.1 h per binner and then failed for all four binners: `NameError: contig_names` for MetaBAT2, GenomeFace and VAMB, `KeyError: 427795` (off-by-one) for SemiBin2. With only the NameError fixed, SemiBin2 still fails; for the other binners the off-by-one changes 88–169 contig assignments, mostly multi-bin calls.
- **The user's patch stops label propagation after one hop.** It marks the unbinned neighbours of each newly labelled contig as processed and then skips their fresh heap entries. On medea-II it propagates 21–31% fewer labels than upstream with the fixes and differs from it for 3,592–6,230 contigs per binner (2.7–7.1%).
- **`graphbin2_fast.py` reproduces upstream + fixes exactly.** On medea-II its CSV and all bin FASTAs are byte-identical to the fixed upstream code run with the same deterministic heap order, for all four binners. With upstream's own tie handling, 2–3 contigs per binner differ (≥ 0.99997 identical). It needs 21–29 s and ≤ 1.6 GB per binner, compared with 3.7–4.9 h and 6.2–7.1 GB for the fixed upstream code (600–840× faster) and 2.6–3.2 min and 6.5 GB for the patch.
- **Tests:** 317 pytest tests pass in about 15 s. A 600 + 600 instance sweep gives 600/600 identical results on tie-free graphs, and 600/600 identical with ties when the heap order matches (562/600 with upstream's tie order). Unmodified upstream raised KeyError in 318/400 tie-free random graphs and NameError in 6.
- **Refinement lowered MetaBAT2 bin quality on this assembly (CheckM2).** HQ bins (≥ 90% completeness, < 5% contamination) dropped from 26 to 12 and MQ-or-better from 95 to 77; the patch gives 12 and 87. In the stage ablations, label removal + refinement alone gives 17 HQ (completeness drops) and propagation alone gives 16 HQ (contamination rises). On this graph neither stage improves the bin set.
- **Recommendation for v1.0.0:** build the graph with the k read from the contig headers and run `graphbin2_fast.py` (no prepResult or shims needed). Make GraphBin2 optional, off by default. When it is enabled, give Binette the original and the refined bin sets so refinement can only add candidates. This is a design decision for the user.

## 1. Zero-edge assembly graphs

MEGAHIT ran with `--k-min 25 --k-max 125 --k-step 10 --min-contig-len 1000`. All 427,795 contig headers of `medea-II_coassembly_29_127.fasta` start with `k125_`, so the final k is 125 (the `29_127` in the file name does not describe the assembly). The old workflow ran `megahit_toolkit contig2fastg 127` and `agtools fastg2gfa -k 127`. Adjacent MEGAHIT contigs overlap by the final k, so with k = 127 no overlap is found and the graph has no edges.

| graph | L_lines | self_loop_links | unique_undirected_edges | isolated_contigs | components_ge2 | largest_component |
|---|---|---|---|---|---|---|
| k = 125 | 97,360 | 676 | 48,310 | 357,256 | 22,981 | 33 |
| k = 127 | 0 | 0 | 0 | 427,795 | 0 | 1 |

The k = 127 rebuild is byte-identical to the GFA used by the old workflow (md5 `6052c5d8a628ea2ee0dc27c71875941c`). The old GraphBin2 logs for all four binners report `Total number of edges in the assembly graph: 0` and `Number of isolated contigs: 427795`, so the previous GraphBin2 step returned the input bins unchanged. With k = 125, 70,539 contigs (16.5%) have at least one neighbour. Because `--min-contig-len 1000` removed the short contigs between them, components are small (largest 33 contigs). The same symptom was reported for ecotip2 (2,833,342 contigs) but was not re-checked here.

Contig2fastg + fastg2gfa took 14 s + 11 s (1.5 GB and 3.8 GB peak RSS) for each k.

## 2. Code review findings (verified)

- (a) User's patch, propagation: confirmed. After labelling contig x, `processed_contigs.update(unbinned_neighbours)` is followed by `if to_bin in processed_contigs: continue`, which also discards the fresh entries pushed for these neighbours. In a 5-contig chain with one binned end, upstream + fixes and the fast code label all 4 unbinned contigs; the patch labels 1 (test `test_chain_propagates_beyond_one_hop`).
- (b) Off-by-one: confirmed. `coverages` and `contig_lengths` are keyed by vertex index but read via `contigs_map[...]` (1-based NODE numbers) in 12 places. In contig2fastg/agtools graphs, segment i is `NODE_{i+1}` in FASTA order (checked on medea-II), so every read returns the next contig's value and vertex n−1 raises `KeyError: n`. The patch reads `coverages[node]`, which is correct.
- (c) `contig_names`: confirmed. NameError whenever at least one multi-labelled contig is found; raised after the full computation and before any output is written.
- (d) BFS order dependence: confirmed. Depths of queued nodes are overwritten on rediscovery, so the result depends on neighbour order. igraph `neighbors(v, mode=ALL)` returned sorted unique neighbours on 200 random graphs (igraph 1.0.0 locally); the medea-II results (identical to the reference run with igraph 0.11.9) confirm this order. The patch iterates Python sets, which is one more source of differences.
- (e) Zero edges: confirmed and explained (Section 1).
- Other: the old rule's `safe_drop` shim removed the row for `k125_<contig count>`, i.e. `k125_427795`. No such contig exists in medea-II (MEGAHIT IDs are not consecutive). The last graph vertex is `k125_2996973` (binned only by SemiBin2), so the shim could not prevent the KeyError. With the patch or the fixed code it is not needed.

## 3. Implementation: `workflow/scripts/graphbin2_fast.py`

- One file, standard library only, Python ≥ 3.9 (run on LiSC with `/usr/bin/python3` 3.9.25 and locally with 3.9.23 and 3.12). BSD-3 licence text of GraphBin2 kept, with a modification notice.
- CLI as GraphBin2 MEGAHIT mode: `--contigs --graph --binned --output --prefix --depth 5 --threshold 1.5 --delimiter , --nthreads` (accepted, single-threaded). `--assembler/--abundance/--paths` are accepted and ignored, for drop-in use. Extras: `--fasta-all-bins`, `--quiet`, and `--binned` may be a folder of bin FASTAs (read like prepResult).
- Outputs as upstream: `<output>/<prefix>graphbin2_output.csv` (CRLF rows, one per membership), `<output>/<prefix>bins/<prefix><bin>.fasta`, `<output>/<prefix>graphbin2.log`. A multi-binned contig goes into the FASTA of its last bin only, as upstream; `--fasta-all-bins` writes it into every bin.
- Deliberate differences: fix (b), fix (c), tolerant bins parsing (first header token, CRLF, blank rows, de-duplication, first bin kept for conflicting rows, contigs absent from the graph counted and ignored), sequence-digest fallback when GFA and FASTA order differ, output folder created and given a trailing `/`, WARNING with the k from the contig names when the graph has 0 edges.
- Speed: arrays/sets/dicts instead of list scans; neighbour tuples cached in ascending vertex order; deque BFS that reproduces upstream's visited-on-dequeue, depth-overwrite and duplicate-enqueue behaviour; components computed once; plain tuple heap with versioned lazy invalidation. The set of live heap entries always equals upstream's filter + heapify + re-push.
- Tie-breaking: heap key = (BFS depth, |coverage difference|, vertex of the contig to bin, bin index, vertex of the labelled contig). Upstream compares only the first two and pops ties in heap-position order. Results are identical when no ties occur, and always identical to upstream's algorithm run with this total order.
- `refine(..., stages=...)` lets a script switch stages off for diagnostics (used for the ablation); the CLI always runs all stages.
- Final file (md5 `d67272e6508fad529886e3e72c06f7b3`): the validation runs used an earlier revision without the bin-folder input and the diagnostic `stages` argument. Re-run on LiSC for MetaBAT2, the final file gave byte-identical CSV and bin FASTAs from the CSV (21 s, 1.0 GB) and directly from the MetaBAT2 bin folder (30 s, 0.9 GB). Locally it reproduced the validated CSVs of all four binners.

## 4. Tests

- `tests/graphbin2_fast/test_graphbin2_fast.py` (317 tests) with `tests/graphbin2_fast/gb2_testlib.py`. Implementations compared: fast; the reference (`tests/graphbin2_fast/reference/graphbin2_MEGAHIT_reference.py` = upstream 1.3.3 with only fixes (b) and (c), generated by text substitution, diff in `reference_vs_upstream.diff`); the reference with the total heap order; the user's patch; unmodified upstream.
- Targeted cases: chain (multi-hop propagation; the patch stops after one hop), star, a multi-labelled contig (fix (c), FASTA last-bin semantics, `--fasta-all-bins`), a 0-edge graph (no-op plus WARNING), tolerant bins parsing, shuffled GFA segment order, bin-folder input equal to prepResult's CSV, CLI file naming, and byte-identical CSV/FASTA output against the reference.
- Randomised: 150 tie-free graphs equal to the reference; 150 graphs with ties equal to the reference with total heap order; agreement with upstream's tie order (141/150 identical, mean 0.9969); the patch matches the reference in only 9/60 tie-free graphs; unmodified upstream raises KeyError and NameError on fixed seeds.
- Real subgraph fixture (685 medea-II contigs from 120 whole components plus 80 isolated contigs, 499 edges, 386 MetaBAT2-binned; synthetic sequences): fast output byte-identical to the reference (CSV and 114 bin FASTAs), 304 contigs changed, patch agreement 0.875, unmodified upstream `KeyError: 685`.
- `tests/graphbin2_fast/sweep_graphbin2.py` reproduces the larger sweep below (seeds 10000–10599 and 20000–20399):

| instances | fast_eq_reference | fast_eq_reference_total_order | mean_agreement_fast_reference | patch_eq_reference | mean_agreement_patch_reference |
|---|---|---|---|---|---|
| 600 tie-free (pairwise-distinct coverage differences) | 600/600 | 600/600 | 1.0000 | 101/600 | 0.676 |
| 600 with ties (small integer coverages) | 562/600 | 600/600 | 0.9942 (min 0.24) | 104/600 | 0.668 |

## 5. Validation on medea-II (k = 125 graph, four binners)

Inputs: `medea-II_coassembly_29_127.fasta`, the rebuilt k = 125 GFA, and the four `*_initial_contig_bins.csv` files (prepResult output, CRLF line endings; the VAMB file carries full FASTA headers). The CSVs were trimmed to the first header token and de-duplicated (no duplicates were present) and given to every implementation. The fast script was also run on the raw CSVs. Runs used one core each; wall time and peak RSS are from `/usr/bin/time -v`.

### 5.1 Run time and memory

| implementation | code | metabat2 | genomeface | semibin2 | vamb |
|---|---|---|---|---|---|
| fast | graphbin2_fast.py (Python 3.9.25) | 28 s, 1.5 GB | 21 s, 0.9 GB | 29 s, 1.1 GB | 21 s, 1.0 GB |
| patch | user's patch as installed (graphbin2-patched env) | 3.2 min, 6.5 GB | 3.2 min, 6.5 GB | 2.6 min, 6.5 GB | 2.6 min, 6.5 GB |
| reference | upstream 1.3.3 + fixes (b),(c) | 4.89 h, 6.3 GB | 3.67 h, 6.2 GB | 4.88 h, 7.1 GB | 4.90 h, 6.2 GB |
| reference_total | reference with total heap order | 4.88 h, 6.3 GB | 3.75 h, 6.2 GB | 4.99 h, 6.2 GB | 5.15 h, 6.2 GB |
| upstream_conly | upstream 1.3.3 + fix (c) only | 4.91 h, 6.2 GB | 3.31 h, 6.2 GB | 3.38 h, 6.1 GB (KeyError) | 5.06 h, 6.2 GB |
| upstream_unfixed | upstream 1.3.3 unmodified | 3.93 h, 6.1 GB (NameError) | 2.86 h, 6.1 GB (NameError) | 3.37 h, 6.1 GB (KeyError) | 4.07 h, 6.1 GB (NameError) |

The reference runs (and the other upstream-derived variants) ran in parallel in one 16-CPU job, one core each. Their run times include 48–63 min of writing outputs (a list scan over all bins for every contig). Fast and patch runs ran in parallel in a 32-CPU job. The fast run on the raw prepResult CSVs (full VAMB headers, CRLF) took 21–29 s and produced byte-identical CSVs and bin FASTAs. The fast run on the k = 127 graph (MetaBAT2) took 21 s, returned the input bins and logged the 0-edge WARNING with the hint to use k = 125.

### 5.2 Changes relative to the initial bins

Per contig: *removed* = binned before, unbinned after; *moved* = label replaced by a different bin; *propagated* = unbinned before, labelled after; *multi-binned* = more than one bin after. Stage counters from the fast log (labels removed in stage 1, including contigs re-labelled later) are in the CSV.

| binner | implementation | binned_initial | binned_final | removed | moved | propagated | multi_binned | changed |
|---|---|---|---|---|---|---|---|---|
| metabat2 | fast | 129,403 | 141,407 | 4,474 | 2,105 | 16,478 | 67 | 23,101 |
| metabat2 | patch | 129,403 | 136,369 | 5,097 | 1,673 | 12,063 | 46 | 18,867 |
| metabat2 | reference | 129,403 | 141,407 | 4,474 | 2,105 | 16,478 | 67 | 23,101 |
| metabat2 | reference_total | 129,403 | 141,407 | 4,474 | 2,105 | 16,478 | 67 | 23,101 |
| metabat2 | upstream_conly | 129,403 | 141,407 | 4,474 | 2,104 | 16,478 | 32 | 23,082 |
| genomeface | fast | 68,815 | 87,426 | 577 | 467 | 19,188 | 32 | 20,245 |
| genomeface | patch | 68,815 | 81,209 | 782 | 285 | 13,176 | 22 | 14,251 |
| genomeface | reference | 68,815 | 87,426 | 577 | 467 | 19,188 | 33 | 20,246 |
| genomeface | reference_total | 68,815 | 87,426 | 577 | 467 | 19,188 | 32 | 20,245 |
| genomeface | upstream_conly | 68,815 | 87,426 | 577 | 468 | 19,188 | 17 | 20,246 |
| semibin2 | fast | 103,086 | 120,483 | 3,685 | 2,043 | 21,082 | 95 | 26,856 |
| semibin2 | patch | 103,086 | 114,575 | 4,408 | 1,430 | 15,897 | 70 | 21,772 |
| semibin2 | reference | 103,086 | 120,483 | 3,685 | 2,043 | 21,082 | 94 | 26,855 |
| semibin2 | reference_total | 103,086 | 120,483 | 3,685 | 2,043 | 21,082 | 95 | 26,856 |
| vamb | fast | 120,343 | 133,288 | 2,309 | 1,421 | 15,254 | 67 | 19,030 |
| vamb | patch | 120,343 | 129,703 | 2,629 | 1,178 | 11,989 | 62 | 15,842 |
| vamb | reference | 120,343 | 133,288 | 2,309 | 1,421 | 15,254 | 66 | 19,029 |
| vamb | reference_total | 120,343 | 133,288 | 2,309 | 1,421 | 15,254 | 67 | 19,030 |
| vamb | upstream_conly | 120,343 | 133,288 | 2,309 | 1,423 | 15,254 | 33 | 19,015 |

### 5.3 Agreement between implementations

Fraction of contigs (labelled in either result) with identical label sets, number of differing contigs, and per-bin Jaccard index of the contig sets.

| binner | comparison | contigs_compared | contigs_differing | fraction_identical | bins_jaccard_lt1 | mean_bin_jaccard | min_bin_jaccard |
|---|---|---|---|---|---|---|---|
| metabat2 | fast_vs_reference | 141,407 | 3 | 0.99998 | 4 | 0.9999 | 0.986 |
| metabat2 | fast_vs_reference_total | 141,407 | 0 | 1.00000 | 0 | 1.0000 | 1.000 |
| metabat2 | patch_vs_reference | 141,407 | 5,059 | 0.96422 | 201 | 0.9538 | 0.448 |
| metabat2 | upstream_conly_vs_reference | 141,407 | 169 | 0.99880 | 56 | 0.9977 | 0.952 |
| genomeface | fast_vs_reference | 87,426 | 3 | 0.99997 | 4 | 0.9999 | 0.990 |
| genomeface | fast_vs_reference_total | 87,426 | 0 | 1.00000 | 0 | 1.0000 | 1.000 |
| genomeface | patch_vs_reference | 87,426 | 6,230 | 0.92874 | 233 | 0.9264 | 0.603 |
| genomeface | upstream_conly_vs_reference | 87,426 | 88 | 0.99899 | 42 | 0.9981 | 0.938 |
| semibin2 | fast_vs_reference | 120,483 | 2 | 0.99998 | 3 | 1.0000 | 0.991 |
| semibin2 | fast_vs_reference_total | 120,483 | 0 | 1.00000 | 0 | 1.0000 | 1.000 |
| semibin2 | patch_vs_reference | 120,483 | 5,934 | 0.95075 | 329 | 0.9520 | 0.781 |
| vamb | fast_vs_reference | 133,288 | 3 | 0.99998 | 5 | 0.9999 | 0.989 |
| vamb | fast_vs_reference_total | 133,288 | 0 | 1.00000 | 0 | 1.0000 | 1.000 |
| vamb | patch_vs_reference | 133,288 | 3,592 | 0.97305 | 255 | 0.9639 | 0.630 |
| vamb | upstream_conly_vs_reference | 133,288 | 122 | 0.99908 | 53 | 0.9985 | 0.942 |

The fast code and the reference with total heap order agree byte for byte (CSV and the md5 of every bin FASTA) for all four binners. The plain reference differs from it only in the comparison function of the heap. The 2–3 contigs per binner that differ from the plain reference are therefore tie-order effects; they show up as one more or one fewer multi-binned contig (e.g. GenomeFace 32 vs 33).

Upstream with only fix (c) (i.e. still off by one) finished for three binners and differs from the reference in 88–169 contigs, mostly multi-bin calls (MetaBAT2: 32 instead of 67 multi-binned contigs). For SemiBin2 it raised `KeyError: 427795` after 3.4 h. Unmodified upstream produced no output for any binner.

The user's patch differs from the reference for 3,592–6,230 contigs per binner. It propagates fewer labels (MetaBAT2 12,063 vs 16,478). Its net number of contigs that lose their label is higher (MetaBAT2 5,097 vs 4,474) because fewer of the contigs whose labels were removed in stage 1 are labelled again during the truncated propagation.

### 5.4 CheckM2 on MetaBAT2 bins

CheckM2 1.1.0 (LiSC module, uniref100.KO.1.dmnd). Bin FASTAs were written from the contig-bin CSVs (a multi-binned contig goes into each of its bins). Identical bins were evaluated once (689 + 488 unique bins). HQ = completeness ≥ 90 and contamination < 5; MQ = completeness ≥ 50 and contamination < 10, excluding HQ. Rows 2 and 3 are stage ablations of the corrected code (`refine(..., stages=...)`), not GraphBin2 modes.

| set | bins_changed | contigs | total_Mbp | HQ | MQ | MQ_or_better | HQ_lost | HQ_gained | median_completeness | median_contamination | mean_contamination |
|---|---|---|---|---|---|---|---|---|---|---|---|
| MetaBAT2 initial bins | 0 | 129403 | 573.3 | 26 | 69 | 95 | 0 | 0 | 66.04 | 3.01 | 9.12 |
| ablation: label removal + refinement only | 225 | 122698 | 541.5 | 17 | 76 | 93 | 10 | 1 | 63.17 | 2.59 | 8.39 |
| ablation: propagation only | 239 | 147839 | 621.7 | 16 | 60 | 76 | 10 | 0 | 70.6 | 4.8 | 11.24 |
| after graphbin2_fast.py (full) | 240 | 141474 | 594.6 | 12 | 65 | 77 | 14 | 0 | 66.33 | 4.72 | 10.92 |
| after user's patch (full) | 240 | 136415 | 574.5 | 12 | 75 | 87 | 14 | 0 | 63.72 | 4.37 | 9.78 |

Of the 14 HQ MetaBAT2 bins lost after full refinement, most gain contamination from propagated contigs. For example, bin.139 gains 152 contigs (1.28 Mbp; contamination 2.8 → 38.8) and bin.327 gains 99 propagated and 41 moved contigs (11.6% contamination, completeness 95.6 → 85.0). Label removal mostly lowers completeness by a few points. In total, refinement removed 6,579 contig memberships (30.7 Mbp) and added 16,497 propagated (42.6 Mbp), 2,101 moved (9.1 Mbp) and 52 multi-bin memberships (0.4 Mbp). Propagated contigs are short (median 1,391 bp vs 2,564 bp for initially binned contigs). Per-bin values are in `checkm2_metabat2_bins.csv` and `metabat2_bin_change_composition.csv`.

This assessment uses one binner, one assembly and CheckM2 estimates only. The reference run was not given to CheckM2; it differs from the fast run by 3 MetaBAT2 contigs.

## 6. Integration notes for v1.0.0

- Graph: read k from the contig headers and fail if they disagree. Pass the original MEGAHIT contigs (with `multi=` and `len=`) to contig2fastg and GraphBin2, not the header-trimmed `.fixed.fa`.
- GraphBin2: one rule per binner, `python workflow/scripts/graphbin2_fast.py --contigs ... --graph ... --binned <bin FASTA folder> --output <dir>`. The prepResult rule, the trim/dedupe and `safe_drop` shims and the `sitecustomize` patch are no longer needed. Any Python ≥ 3.9 image works (e.g. the stock Bioconda graphbin2 1.3.3 container or python:3.12-slim). Measured: 1 core, ≤ 1.6 GB, < 30 s for 428k contigs. Default resources of 1 CPU, 4 GB and 30 min leave a margin. For ecotip2 (2.8 M contigs) about 10 GB is an extrapolation, not a measurement.
- Empty bin FASTAs are still possible when a bin loses all its contigs (none did on medea-II). Keep a prune step before Binette.
- Binette: given the CheckM2 result, pass both the original and the GraphBin2-refined bin sets (8 sets) if GraphBin2 is enabled, or skip GraphBin2 by default.

```bash
# k from the MEGAHIT contig names (k125_123 -> 125); all headers must agree
K=$(grep '^>' {input.contigs} | sed -E 's/^>k([0-9]+)_.*/\1/' | sort -u)
[ "$(echo "$K" | wc -l)" -eq 1 ] || {{ echo "contig headers carry several k values: $K" >&2; exit 1; }}
megahit_toolkit contig2fastg "$K" {input.contigs} > {output.fastg}
agtools fastg2gfa -g {output.fastg} -k "$K" -o {params.outdir}/ && mv {params.outdir}/converted_graph.gfa {output.gfa}

# per binner (bin FASTA folder or prepResult-style CSV)
python workflow/scripts/graphbin2_fast.py --contigs {input.contigs} --graph {input.gfa} \
    --binned {input.bins_dir} --output {params.outdir}/ --quiet
```

## 7. Not assessed

- CheckM2 of GenomeFace, SemiBin2 and VAMB bins after refinement; effect on Binette/dRep output and on final MAGs (GTDB-Tk, MIMAG).
- ecotip2 (2,833,342 contigs): neither the graph k nor GraphBin2 run time/memory was checked; the ~10 GB memory figure is an extrapolation.
- Parameter sensitivity (`--depth`, `--threshold`; defaults 5 and 1.5 used throughout).
- Whether the extra contamination that CheckM2 reports for propagated contigs is real (no GUNC or read-based check).
- Other GraphBin2 modes (SPAdes, SGA, Flye) were not re-implemented; their source has the same list-based structure.
- Upstream's tie order was reproduced only in the reference runs; the fast code deliberately uses a total order.

## 8. Files and LiSC locations

- Artifacts: `graphbin2_fast.py`, `test_graphbin2_fast.py`, `gb2_testlib.py`, `sweep_graphbin2.py`, `graphbin2_MEGAHIT_reference.py`, `reference_vs_upstream.diff`, `graphbin2_MEGAHIT_fix_b_c.patch`, fixture files `medea2_subgraph_*`, comparison tables `graphbin2_medea2_*.csv`, `checkm2_metabat2_*.csv`, `metabat2_bin_change_composition.csv`, `graph_stats_medea2_k125_k127.json`, `sweep_graphbin2.json`, the component archive `graphbin2_component.tar.gz` and the issue draft.
- LiSC `<scratch>/bfm_v1_dev/graphbin2_validation/`: `graph_k125/` and `graph_k127/` (FASTG + GFA), `inputs/` (copies of the four initial-bins CSVs, `clean/` trimmed versions), `runs/<variant>/<binner>/` (fast, fast_rawcsv, fast_k127 and patch with bin FASTAs; reference variants without FASTAs), `runs/ablation/`, `checkm2/metabat2/` and `checkm2/metabat2_ablation/` (inputs, CheckM2 output, `bins_map.tsv`), `scripts/`, `skeleton/`, `graph_stats.json`, `component/` (final code).
