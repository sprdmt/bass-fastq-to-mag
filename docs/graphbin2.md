# Assembly graph and GraphBin2

GraphBin2 (Mallawaarachchi et al. 2021) refines the bins of a binner with the assembly graph:
it removes labels that disagree with the graph neighbourhood and propagates labels to unbinned
neighbours. This page describes how the workflow builds the graph, why it ships its own
implementation of GraphBin2, how that implementation was validated, and what refinement did to
bin quality on a real co-assembly. The full validation report is
[graphbin2/validation_report.md](graphbin2/validation_report.md).

## What the workflow runs

| step | rule | what it does |
|---|---|---|
| graph | `contig2fastg`, `fastg2gfa` | reads the final k from the MEGAHIT contig names (`>k125_…` → 125; all names must agree), runs `megahit_toolkit contig2fastg <k>` and `agtools fastg2gfa -k <k>`, and stops if the GFA has no edges (`graph.require_edges: true`) |
| refinement | `graphbin2` (one job per group and binner) | runs `workflow/scripts/graphbin2_fast.py` on the binner's contig-to-bin table and writes `results/graphbin2/<group>/<binner>/contig2bin.tsv` |
| consensus | `binette` | gets the original and the refined tables (`graphbin2.binette_input: both`, default) or only the refined ones (`refined`) |

```yaml
graph:
  require_edges: true     # stop on a graph without edges
graphbin2:
  enabled: true           # false: Binette gets only the binners' own bins
  binette_input: both     # both | refined
  depth: 5                # GraphBin2 --depth
  threshold: 1.5          # GraphBin2 --threshold
```

With `binette_input: both`, Binette sees eight bin sets (four binners, each before and after
refinement) and keeps the better version of each genome. Refinement can then only add
candidates; it cannot replace a good bin with a worse one. The `origin` column of the Binette
report shows which set each final bin came from.

## The zero-edge graph

Adjacent MEGAHIT contigs overlap by exactly the final k. If `contig2fastg` and `fastg2gfa` get
a different k, no overlap is found and the graph has no edges. GraphBin2 then leaves all bins
unchanged without any error. This happened in the earlier version of this workflow: the medea-II
co-assembly was run with `--k-max 125`, its contig names all start with `k125_`, but the graph
was built with k = 127.

| graph (medea-II, 427,795 contigs) | links in GFA | unique edges | isolated contigs | components ≥ 2 contigs | largest component |
|---|---|---|---|---|---|
| k = 125 (correct) | 97,360 | 48,310 | 357,256 | 22,981 | 33 |
| k = 127 (as used before) | 0 | 0 | 427,795 | 0 | 1 |

The k = 127 rebuild is byte-identical to the GFA of the old run, and the old GraphBin2 logs
report `Total number of edges in the assembly graph: 0` for all four binners. Even the correct
graph is sparse: `--min-contig-len 1000` removes the short contigs that connect the long ones,
so only 16.5% of the contigs have a neighbour. The workflow therefore reads k from the contig
names and fails on a graph without edges.

## Why a re-implementation

Upstream GraphBin2 1.3.3 in MEGAHIT mode could not process medea-II:

* **Off-by-one (bug b).** Coverages and lengths are stored by graph vertex index (0-based) but
  read through the 1-based NODE numbers in 12 places. Every lookup returns the next contig's
  value, and the last vertex raises `KeyError` (SemiBin2 on medea-II).
* **Undefined name (bug c).** As soon as one contig ends up in more than one bin,
  `NameError: contig_names` is raised after the full computation and before any output is
  written (MetaBAT2, GenomeFace and VAMB on medea-II).
* **Run time.** Lists are scanned where sets and dictionaries are needed: with both bugs fixed,
  one binner took 3.7–4.9 h and 6.2–7.1 GB on medea-II.

The patched module used previously avoided both errors, but it marks the unbinned neighbours of
each newly labelled contig as processed and then skips their queued entries, so labels travel
only one hop. On medea-II it propagated 21–31% fewer labels than upstream with the two fixes and
differed from it for 3,592–6,230 contigs per binner.

`workflow/scripts/graphbin2_fast.py` implements the upstream algorithm with fixes (b) and (c)
and uses only the Python standard library. It keeps the CLI and output files of GraphBin2
MEGAHIT mode (`--contigs --graph --binned --output --prefix --depth --threshold --delimiter`)
and adds a few practical changes:

* it accepts a folder of bin FASTAs as `--binned`;
* it tolerates CRLF line ends, full FASTA headers and duplicate rows in the bins file;
* it falls back to matching sequences when the GFA and FASTA orders differ;
* it warns, giving the k from the contig names, when the graph has no edges.

Heap ties are broken by a total order (depth, coverage difference, vertex, bin, source vertex),
so the result does not depend on heap layout.

## Validation

| implementation | result on medea-II (four binners) | time per binner | peak memory |
|---|---|---|---|
| upstream 1.3.3, unmodified | no output (NameError ×3, KeyError ×1) | 2.9–4.1 h until the error | 6.1 GB |
| upstream + fixes (b) and (c) | reference | 3.7–4.9 h | 6.2–7.1 GB |
| previous patch | 2.7–7.1% of contigs differ from the reference | 2.6–3.2 min | 6.5 GB |
| `graphbin2_fast.py` | byte-identical to the reference run with the same heap order; 2–3 contigs per binner differ under upstream's own tie order | 21–29 s | ≤ 1.6 GB |

The test suite (`tests/graphbin2_fast/`) compares the fast code with the reference implementation
(`tests/graphbin2_fast/reference/`, upstream 1.3.3 with only the two fixes, see
[graphbin2/reference_vs_upstream.diff](graphbin2/reference_vs_upstream.diff)). It covers
targeted cases (multi-hop chains, multi-bin contigs, a graph without edges, tolerant parsing,
shuffled GFA order), 300 random graphs, and a 685-contig subgraph of medea-II with synthetic
sequences:

```bash
pytest -q tests/graphbin2_fast                  # needs python-igraph for the reference comparisons
python tests/graphbin2_fast/sweep_graphbin2.py  # the 1,200-instance random sweep of the report
```

In the sweep, the fast code equalled the reference in 600 of 600 tie-free graphs and, with the
same heap order, in 600 of 600 graphs with ties. Tests that need the previous patch or an
installed upstream GraphBin2 skip when these are absent.

## Effect on bin quality

On medea-II, CheckM2 1.1.0 was run on the MetaBAT2 bins before and after refinement. HQ means
≥ 90% completeness and < 5% contamination; MQ means ≥ 50% and < 10%, excluding HQ. The two
middle rows switch off one stage of the algorithm.

| bin set | HQ | MQ or better | median completeness | median contamination |
|---|---|---|---|---|
| MetaBAT2 bins | 26 | 95 | 66.0 | 3.0 |
| label removal and refinement only | 17 | 93 | 63.2 | 2.6 |
| propagation only | 16 | 76 | 70.6 | 4.8 |
| full GraphBin2 (`graphbin2_fast.py`) | 12 | 77 | 66.3 | 4.7 |
| previous patch | 12 | 87 | 63.7 | 4.4 |

On this sparse graph, removal costs a few points of completeness, and propagation adds short
contigs that raise contamination. For example, one HQ bin gained 152 contigs (1.28 Mbp) and its
contamination rose from 2.8% to 38.8%. This is why the default gives Binette both versions of
every bin set. To skip refinement entirely, set `graphbin2.enabled: false`.

The workflow's cluster test points the same way. It used two small medea-II co-assemblies
(2 × 4 million read pairs each) with all four binners and the default `binette_input: both`.
All 8 final MAGs were closest to a bin of the binners' own sets (Jaccard index of contig sets
0.82-0.99); their GraphBin2-refined counterparts matched less well (0.73-0.95). GraphBin2 took
1-2 s per binner there, and Binette with eight bin sets took 98-216 s.

This assessment covers one assembly, one binner and CheckM2 estimates only. It does not cover:

* the other binners;
* the effect on the Binette and dRep results;
* other `--depth` and `--threshold` values;
* larger assemblies such as ecotip2 (2.8 million contigs), where memory of about 10 GB is an
  extrapolation.

## Upstream report

[graphbin2/upstream_issue_draft.md](graphbin2/upstream_issue_draft.md) describes bugs (b) and (c)
for the GraphBin2 maintainers, and
[graphbin2/graphbin2_MEGAHIT_fix_b_c.patch](graphbin2/graphbin2_MEGAHIT_fix_b_c.patch) fixes
both in the upstream code.
