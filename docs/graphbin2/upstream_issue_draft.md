# Draft issue for metagentools/GraphBin2

**Title:** MEGAHIT mode: coverage/length lookups are off by one vertex (KeyError or wrong values), NameError for multi-labelled contigs, and quadratic run time

**Version:** GraphBin2 1.3.3 (PyPI; `src/graphbin2/graphbin2_MEGAHIT.py` is identical on master), python-igraph 0.11.9, cogent3 2025.7. Graphs were built with `megahit_toolkit contig2fastg` (MEGAHIT 1.2.9) and `agtools fastg2gfa` (agtools 0.1.2), which name segments `NODE_1_length_…` to `NODE_n_length_…` in contig order.

## 1. Off-by-one in coverage and length lookups

`coverages` and `contig_lengths` are filled with the graph **vertex index** as key:

```python
contig_num = contigs_map_rev[graph_to_contig_map_rev[strings[0]]]   # vertex index 0..n-1
contig_lengths[contig_num] = length
coverages[contig_num] = coverage
```

All reads use the **NODE number** instead, e.g. `coverages[contigs_map[node]]` in `runBFS`, `contig_lengths[contigs_map[my_node]] < 10000` in the label-removal step, the bin coverage sums, and four lookups in `is_multi` (12 lookups in total). `contigs_map` maps vertex index → NODE number, which is 1-based, so every lookup returns the values of the next contig. Vertex n−1 raises `KeyError: n`.

What we observed:
- Of 400 random test graphs (MEGAHIT-style headers, distinct coverages), 318 runs raised `KeyError` and 6 raised `NameError` (item 2). 75 runs finished with the same result as the fixed code, and 1 finished with a different result.
- A 685-contig subgraph of a marine co-assembly (real headers, topology and MetaBAT2 bins) raised `KeyError: 685`.
- Full co-assembly (427,795 contigs ≥ 1 kb, 48,310 edges, four binners): on the graph built with the correct k (48,310 edges), the unmodified code ran 2.9–4.1 h per binner and then raised `NameError` (item 2) for three binners and `KeyError: 427795` for the fourth, whose last vertex is binned. With only item 2 fixed, that binner still raised `KeyError`, and for the other three the off-by-one changed 88–169 contig assignments, mostly multi-bin calls (e.g. 32 instead of 67 multi-labelled contigs).

Fix: index by vertex, `coverages[node]` and `contig_lengths[node]` (diff below).

## 2. `NameError: name 'contig_names' is not defined`

The multi-bin logging loop uses `contig_names[contig]`, which is never defined. The error is raised after all computation, as soon as at least one multi-labelled contig is found, so no output is written. Fix: `graph_to_contig_map[contigs_map[contig]]`.

## 3. Run time

Membership tests on lists (`binned_contigs`, `isolated`, `non_isolated`, `bins[n]`, `unbinned_contigs`), `queue.pop(0)`, the connected-component search by repeated list scans, and filtering the whole heap after each propagated label make the MEGAHIT mode roughly quadratic. On the 427,795-contig co-assembly the fixed code needed 3.7–4.9 h and 6.2–7.1 GB (one core) per binner. A set/dict re-implementation of the same algorithm (same BFS semantics and neighbour order) produced byte-identical output (CSV and every bin FASTA) for all four binners when heap ties are broken by the same total order (with the current tie handling, 2–3 contigs per binner differ; see item 4). It ran in 21–29 s with 0.9–1.6 GB peak memory.

One trap for anyone optimising the propagation step: after a contig is labelled, upstream removes the old heap entries of its unbinned neighbours and pushes fresh BFS entries for them. A "processed contigs" set that skips stale entries must not also skip the fresh ones. Otherwise propagation stops one hop away from the initial bins. A per-contig version counter checked at pop time reproduces upstream exactly.

## 4. Smaller observations

- `DataWrap.__lt__` compares only (BFS depth, coverage difference). Equal keys are popped in heap order, which depends on insertion order and set iteration. A total order (adding vertex index and bin index) would make results reproducible. In our tests, tied cases agreed with upstream on 94% of random graphs (562/600).
- A multi-labelled contig is written to the FASTA of only one bin (`final_bins` keeps the last bin), while the CSV lists all its bins. Tools that read the bin FASTAs (e.g. Binette) never see the overlap.
- Output paths are built as `f"{output_path}{prefix}bins/"` and `f"{output_path}{prefix}graphbin2_output.csv"`, so `--output out` (without a trailing slash) writes `outbins/` and `outgraphbin2_output.csv`. The log file uses `/`.
- With the default `prepResult`, contig names in the initial-bins CSV are full FASTA headers (e.g. `k125_1 flag=0 multi=3.0000 len=1139` for some binners), which raises `KeyError` in MEGAHIT mode. Taking the first whitespace-separated token would avoid this.
- When the graph is built with a k different from the assembly's final k (e.g. `contig2fastg 127` on a `--k-max 125` assembly), the GFA has no L lines and GraphBin2 silently returns the input binning. A warning when the graph has 0 edges would help. We ran into this ourselves.

## Patch for items 1 and 2

```diff
--- a/src/graphbin2/graphbin2_MEGAHIT.py
+++ b/src/graphbin2/graphbin2_MEGAHIT.py
@@ -335,4 +335,4 @@
                         abs(
-                            coverages[contigs_map[node]]
-                            - coverages[contigs_map[active_node]]
+                            coverages[node]
+                            - coverages[active_node]
                         ),
@@ -410,3 +410,3 @@
                             and BFS_labelled_bin_counts[max_index] > 1
-                            and contig_lengths[contigs_map[my_node]] < 10000
+                            and contig_lengths[my_node] < 10000
                         ):
@@ -657,6 +657,6 @@
                 bin_cov_sum[i] += (
-                    coverages[contigs_map[bins[i][j]]]
-                    * contig_lengths[contigs_map[bins[i][j]]]
+                    coverages[bins[i][j]]
+                    * contig_lengths[bins[i][j]]
                 )
-                bin_contig_len_total[i] += contig_lengths[contigs_map[bins[i][j]]]
+                bin_contig_len_total[i] += contig_lengths[bins[i][j]]
 
@@ -734,3 +734,3 @@
         logger.info(
-            contig_names[contig]
+            graph_to_contig_map[contigs_map[contig]]
             + " belongs to bins "
@@ -841,6 +841,6 @@
         bin_coverages[contig_bin] = bin_coverages[contig_bin] - (
-            coverages[contigs_map[contig]] * contig_lengths[contigs_map[contig]]
+            coverages[contig] * contig_lengths[contig]
         )
         bin_contig_lengths[contig_bin] = (
-            bin_contig_lengths[contig_bin] - contig_lengths[contigs_map[contig]]
+            bin_contig_lengths[contig_bin] - contig_lengths[contig]
         )
@@ -862,3 +862,3 @@
                     neighbour_bins[n].append(neighbour)
-                    neighbour_bin_coverages[n].append(coverages[contigs_map[neighbour]])
+                    neighbour_bin_coverages[n].append(coverages[neighbour])
                     break
@@ -891,3 +891,3 @@
 
-                cov_diff = abs(comb_cov_total - coverages[contigs_map[contig]])
+                cov_diff = abs(comb_cov_total - coverages[contig])
 
@@ -900,3 +900,3 @@
                 and len(min_diff_combination) > 1
-                and contig_lengths[contigs_map[contig]] > 1000
+                and contig_lengths[contig] > 1000
             ):
```

We can open a PR with these two fixes and regression tests (random graphs, a chain that needs multi-hop propagation, and a multi-labelled case). If useful, we can also contribute the faster implementation.
