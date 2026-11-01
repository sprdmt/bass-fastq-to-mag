# medea2_subgraph test fixture

Real topology, real contig headers and real initial bins, synthetic sequences.

- Source: medea-II MEGAHIT co-assembly (427,795 contigs >= 1 kb, k25-125, final k = 125),
  graph rebuilt with `megahit_toolkit contig2fastg 125` + `agtools fastg2gfa -k 125`
  (48,310 unique undirected edges), MetaBAT2 initial bins (prepResult CSV).
- Selection (random seed 2026): 60 connected components with at least 3 contigs that
  contain contigs of two or more MetaBAT2 bins, 60 further components with binned and
  unbinned contigs, and 80 isolated contigs. Whole components were taken, so no edge is cut.
- 685 contigs, 499 edges, 386 binned contigs (CRLF line endings as written by prepResult).
- Contig names, `multi=` and `len=` header fields and the edges are real; sequences are
  random unique 24-mers (identical in `contigs.fa` and `graph.gfa`) to keep the fixture small.
