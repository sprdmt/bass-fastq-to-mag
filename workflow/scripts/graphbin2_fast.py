#!/usr/bin/env python3
"""graphbin2_fast.py - corrected, faster GraphBin2 1.3.3 (MEGAHIT mode).

A self-contained re-implementation of the MEGAHIT mode of GraphBin2 1.3.3
(src/graphbin2/graphbin2_MEGAHIT.py, https://github.com/metagentools/GraphBin2).
Standard library only (Python >= 3.9); it does not import graphbin2, igraph or
cogent3, so it runs inside the stock Bioconda graphbin2 container or any Python 3.

Algorithm (unchanged from upstream): (1) iteratively remove labels of binned
contigs whose depth-limited BFS neighbourhood contains no contig of their own
bin; (2) iteratively move inconsistent contigs to the bin with the highest
distance-weighted support (1/2^depth) if it exceeds `threshold` x own support;
(3) propagate labels to unbinned contigs in components that contain binned
contigs, closest (BFS depth, then |coverage difference|) first; (4) assign
contigs to several bins when the summed mean coverage of a combination of
neighbouring bins matches the contig coverage best.

Deliberate differences from upstream 1.3.3
------------------------------------------
* FIX coverage/length lookup: upstream stores coverages and lengths keyed by
  graph node index but reads them as coverages[contigs_map[node]], where
  contigs_map maps node index -> NODE number (1-based in contig2fastg/agtools
  graphs). Every coverage difference, the 10 kb length test, the bin mean
  coverages and the multi-bin test therefore used the NEXT contig's values, and
  the last contig raised KeyError. Here they are looked up by node index.
* FIX multi-bin logging: upstream logs contig_names[contig], an undefined name
  (NameError as soon as a multi-labelled contig is found). The contig name is
  logged instead.
* Tolerant parsing of the initial bins file: the contig field is cut at the
  first whitespace (prepResult copies full FASTA headers), CR/LF line endings,
  blank and one-column rows are accepted, identical rows are de-duplicated, a
  contig listed under several bins keeps its first bin (WARNING), and contigs
  absent from the contigs file/graph are counted and ignored (upstream crashes
  with KeyError). Bin labels are stripped of surrounding whitespace; a bin
  none of whose contigs is in the graph gets no output FASTA. --binned may
  also be a folder of bin FASTAs (labels = file names, as prepResult).
* Segment-to-contig mapping: as upstream, segments are paired with contigs by
  position and accepted when the sequences are identical (compared through a
  128-bit BLAKE2 digest instead of keeping every sequence twice in memory);
  segments that do not match by position are then matched by sequence digest.
  Unmatched segments abort the run; contigs missing from the graph are
  reported and ignored (upstream crashes).
* The output folder gets a trailing '/' and is created if needed (upstream
  concatenates '<output><prefix>graphbin2_output.csv' without a separator).
* Logging: node, link, edge, isolated-node and component counts, the MEGAHIT
  k implied by the contig names, a WARNING when the graph has 0 edges (then
  GraphBin2 cannot change anything; usually contig2fastg/fastg2gfa were run
  with a k different from the assembly's final k), and per-step summaries.
* --nthreads is accepted for command-line compatibility but not used: the
  run is single-threaded (upstream's thread pool gives no speed-up under the
  GIL either).

Performance changes (no change in results)
------------------------------------------
Set/array membership instead of list scans, neighbour lists cached once in
ascending vertex order (the order igraph's Graph.neighbors() returns, which
the BFS depth bookkeeping depends on), deque-based BFS, connected components
computed once, plain-tuple heap entries, and versioned lazy invalidation in
label propagation: when a contig is labelled, the version of each still
unbinned neighbour is incremented and fresh BFS candidates are pushed with the
new version; older entries for that contig are skipped when popped. The set of
live heap entries is therefore always the one upstream obtains by filtering the
heap, heapifying and re-pushing, so labels keep propagating beyond one hop.

Tie-breaking
------------
Upstream orders propagation candidates only by (BFS depth, |coverage
difference|); equal keys pop in an order determined by heap positions (set
iteration and push order). Here the heap key is the total order
(depth, |coverage difference|, node index of the contig to bin, bin index,
node index of the labelled contig), so results are deterministic and do not
depend on set iteration order. With no ties the result equals upstream's;
with ties it equals upstream's algorithm run with this total order.

Outputs (names as upstream)
---------------------------
<output>/<prefix>graphbin2_output.csv     contig,bin (one row per membership)
<output>/<prefix>bins/<prefix><bin>.fasta one FASTA per bin (empty bins included)
<output>/<prefix>graphbin2.log
As upstream, a contig assigned to several bins is written to the FASTA of the
last of these bins only (sorted bin order) unless --fasta-all-bins is given;
the CSV always lists every membership.

-------------------------------------------------------------------------------
Derived from GraphBin2 (https://github.com/metagentools/GraphBin2), licensed
under the BSD 3-Clause License:

Copyright (c) 2022, Vijini Mallawaarachchi
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Modification notice: this file is a modified re-implementation of
graphbin2_MEGAHIT.py from GraphBin2 1.3.3 (authors: Vijini Mallawaarachchi,
Anuradha Wickramarachchi, Yu Lin), written for the bass-fastq-to-mag workflow
in 2026. Changes are listed above. It is not an official GraphBin2 release.
-------------------------------------------------------------------------------
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import heapq
import itertools as it
import logging
import os
import re
import sys
import time
from collections import Counter, deque

__version__ = "1.3.3+bfm.1"
UPSTREAM_VERSION = "1.3.3"

NODE_RE = re.compile(r"NODE_(.*)_length")
K_PREFIX_RE = re.compile(r"^k(\d+)_")

DELIMITERS = {
    ",": ",", "comma": ",",
    ";": ";", "semicolon": ";",
    "|": "|", "pipe": "|",
    "\t": "\t", "\\t": "\t", "tab": "\t", "$'\\t'": "\t",
    " ": " ", "space": " ", '" "': " ",
}


class GraphBin2Error(RuntimeError):
    """Raised for unusable inputs (reported as an ERROR and exit code 1)."""


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

def _open_text(path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def read_fasta(path):
    """Yield (label, sequence) like cogent3's MinimalFastaParser.

    The label is the header line without '>' (stripped); all whitespace is
    removed from sequence lines, which are concatenated; blank lines and lines
    starting with '#' are ignored.
    """
    label = None
    chunks = []
    with _open_text(path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line[0] == "#":
                continue
            if line[0] == ">":
                if label is not None:
                    yield label, "".join(chunks)
                label = line[1:].strip()
                chunks = []
            else:
                if " " in line or "\t" in line:
                    line = "".join(line.split())
                chunks.append(line)
    if label is not None:
        yield label, "".join(chunks)


def _digest(seq):
    return hashlib.blake2b(seq.encode("ascii", "replace"), digest_size=16).digest()


class Contigs:
    """Contig names, lengths and coverages in FASTA order (MEGAHIT headers)."""

    def __init__(self):
        self.names = []
        self.lengths = []
        self.coverages = []
        self.digests = []
        self.index = {}
        self.len_from_seq = 0


def load_contigs(path):
    c = Contigs()
    for label, seq in read_fasta(path):
        toks = label.split()
        if not toks:
            raise GraphBin2Error(f"empty FASTA header in {path}")
        name = toks[0]
        length = cov = None
        for t in toks[1:]:
            if t.startswith("len="):
                length = int(t[4:])
            elif t.startswith("multi="):
                # upstream: coverage = int(float(strings[2][6:]))
                cov = int(float(t[6:]))
        if cov is None:
            raise GraphBin2Error(
                f"contig header '{label}' has no 'multi=' field; MEGAHIT mode needs "
                "the original MEGAHIT contig headers (k<K>_<id> flag= multi= len=)")
        if length is None:
            length = len(seq)
            c.len_from_seq += 1
        if name in c.index:
            raise GraphBin2Error(f"duplicate contig name '{name}' in {path}")
        c.index[name] = len(c.names)
        c.names.append(name)
        c.lengths.append(length)
        c.coverages.append(cov)
        c.digests.append(_digest(seq))
    return c


class Graph:
    """Undirected, simple assembly graph (self-loops and multi-edges removed)."""

    def __init__(self):
        self.seg_numbers = []
        self.seg_digests = []
        self.n_links = 0
        self.n_self_loops = 0
        self.n_links_unknown = 0
        self.adj = []

    @property
    def n(self):
        return len(self.seg_numbers)

    @property
    def n_edges(self):
        return sum(len(a) for a in self.adj) // 2


def load_gfa(path):
    g = Graph()
    links = []
    with _open_text(path) as fh:
        for line in fh:
            if line.startswith("S"):
                strings = line.split()
                m = NODE_RE.search(strings[1])
                if m is None:
                    raise GraphBin2Error(
                        f"segment name '{strings[1]}' does not follow the NODE_<n>_length_ "
                        "pattern of megahit_toolkit contig2fastg + agtools fastg2gfa graphs")
                g.seg_numbers.append(int(m.group(1)))
                g.seg_digests.append(_digest(strings[2]) if len(strings) > 2 else b"")
            elif line.startswith("L"):
                strings = line.rstrip("\r\n").split("\t")
                m1 = NODE_RE.search(strings[1])
                m2 = NODE_RE.search(strings[3])
                if m1 is None or m2 is None:
                    raise GraphBin2Error(f"cannot parse link line: {line[:200]!r}")
                links.append((int(m1.group(1)), int(m2.group(1))))
    node_of = {}
    for idx, num in enumerate(g.seg_numbers):
        if num in node_of:
            raise GraphBin2Error(f"duplicate segment number NODE_{num} in {path}")
        node_of[num] = idx
    nb = {}
    for a, b in links:
        g.n_links += 1
        if a == b:
            g.n_self_loops += 1
            continue
        ia = node_of.get(a)
        ib = node_of.get(b)
        if ia is None or ib is None:
            g.n_links_unknown += 1
            continue
        nb.setdefault(ia, set()).add(ib)
        nb.setdefault(ib, set()).add(ia)
    adj = [()] * g.n
    for v, s in nb.items():
        adj[v] = tuple(sorted(s))  # ascending vertex ids = igraph neighbors() order
    g.adj = adj
    return g


def map_segments_to_contigs(g, c):
    """Return node -> contig index; position first (as upstream), then digest."""
    n, m = g.n, len(c.names)
    node_contig = [-1] * n
    contig_node = [-1] * m
    for i in range(min(n, m)):
        if g.seg_digests[i] == c.digests[i]:
            node_contig[i] = i
            contig_node[i] = i
    n_positional = sum(1 for x in node_contig if x >= 0)
    unmatched = [i for i in range(n) if node_contig[i] < 0]
    n_by_digest = 0
    if unmatched:
        pool = {}
        for j in range(m):
            if contig_node[j] < 0:
                pool.setdefault(c.digests[j], deque()).append(j)
        for i in unmatched:
            q = pool.get(g.seg_digests[i])
            if q:
                j = q.popleft()
                node_contig[i] = j
                contig_node[j] = i
                n_by_digest += 1
    still = [i for i in range(n) if node_contig[i] < 0]
    if still:
        raise GraphBin2Error(
            f"{len(still)} of {n} graph segments have no identical sequence in the contigs "
            f"file (first: NODE_{g.seg_numbers[still[0]]}); the graph was not built from "
            "these contigs")
    return node_contig, contig_node, n_positional, n_by_digest


BIN_EXTENSIONS = (".fasta", ".fa", ".fna")


def bins_dir_rows(path):
    """Rows (contig header, bin label) from a folder of bin FASTAs, as prepResult.

    prepResult uses the file name (with extension) as bin label and the full
    FASTA header as contig field; files are read in sorted name order here.
    """
    for fname in sorted(os.listdir(path)):
        if fname.lower().endswith(BIN_EXTENSIONS):
            for label, _ in read_fasta(os.path.join(path, fname)):
                yield [label, fname]


def load_bins(path, delimiter, name_to_node):
    """Tolerant reader for the GraphBin2 initial-bins CSV (contig,bin).

    `path` may also be a folder of bin FASTA files (read like prepResult).
    """
    stats = Counter()
    absent_examples = []
    first_bin = {}
    if os.path.isdir(path):
        rows = bins_dir_rows(path)
        fh = None
    else:
        fh = open(path, newline="")
        rows = csv.reader(fh, delimiter=delimiter)
    try:
        for row in rows:
            if not row or all(not x.strip() for x in row):
                stats["blank_rows"] += 1
                continue
            if len(row) < 2:
                stats["rows_without_bin"] += 1
                continue
            stats["rows"] += 1
            raw = row[0]
            toks = raw.split()
            if not toks:
                stats["rows_without_contig"] += 1
                continue
            name = toks[0]
            if name != raw:
                stats["contig_names_trimmed"] += 1
            label = row[1].strip()
            if not label:
                stats["rows_without_bin"] += 1
                continue
            node = name_to_node.get(name)
            if node is None:
                stats["contigs_absent_from_graph"] += 1
                if len(absent_examples) < 5:
                    absent_examples.append(name)
                continue
            prev = first_bin.get(node)
            if prev is None:
                first_bin[node] = label
            elif prev == label:
                stats["duplicate_rows"] += 1
            else:
                stats["contigs_in_several_bins_kept_first"] += 1
    finally:
        if fh is not None:
            fh.close()
    return first_bin, stats, absent_examples


# ---------------------------------------------------------------------------
# GraphBin2 core
# ---------------------------------------------------------------------------

def make_bfs(adj, bin_of, max_depth):
    """Depth-limited BFS of upstream runBFS, reproduced exactly.

    Upstream marks nodes visited when they are dequeued, overwrites the depth
    of queued-but-unvisited nodes on rediscovery (even beyond the depth limit),
    may enqueue a node several times, and stops at labelled nodes. Results
    depend on neighbour order, so adjacency lists are in ascending vertex order
    as returned by igraph. Returns {labelled node: (recorded depth, its bin)}.
    """

    def run_bfs(node):
        queue = deque((node,))
        pop = queue.popleft
        push = queue.append
        visited = set()
        vadd = visited.add
        depth = {node: 0}
        found = {}
        while queue:
            a = pop()
            vadd(a)
            b = bin_of[a]
            if b >= 0 and len(visited) > 1:
                found[a] = (depth[a], b)
            else:
                da = depth[a] + 1
                for nb in adj[a]:
                    if nb not in visited:
                        depth[nb] = da
                        if da > max_depth:
                            continue
                        push(nb)
        return found

    return run_bfs


ALL_STAGES = ("remove", "refine", "propagate", "multibin")


def refine(adj, coverages, lengths, bin_of, n_bins, depth=5, threshold=1.5, log=None,
           stages=ALL_STAGES):
    """Run the four GraphBin2 stages in place on bin_of (node -> bin index or -1).

    Returns (memberships, info): memberships[node] is a sorted tuple of bin
    indices (empty if unbinned); info holds per-stage counts. `stages` selects
    stages for diagnostics only (default: all, as GraphBin2).
    """
    log = log or logging.getLogger("graphbin2_fast.null")
    n = len(adj)
    info = {}
    run_bfs = make_bfs(adj, bin_of, depth)
    t_stage = time.time()

    # -- 1. remove labels of unsupported vertices -------------------------
    removed = 0
    passes = 0
    while "remove" in stages:
        passes += 1
        remove = []
        for v in range(n):
            my_bin = bin_of[v]
            if my_bin < 0 or not adj[v]:
                continue
            found = run_bfs(v)
            if not found:
                continue
            counts = {}
            for _, b in found.values():
                counts[b] = counts.get(b, 0) + 1
            if counts.get(my_bin, 0) == 0:
                remove.append(v)
            elif len(counts) == 1:
                # upstream branch "all labelled BFS nodes in one bin"; with the own
                # bin present it is that bin, so this never removes (kept verbatim)
                (max_index, max_count), = counts.items()
                if max_index != my_bin and max_count > 1 and lengths[v] < 10000:
                    remove.append(v)
        log.debug(f"Remove-labels pass {passes}: {len(remove)} labels removed")
        if not remove:
            break
        for v in remove:
            bin_of[v] = -1
        removed += len(remove)
    info["removed"] = removed
    info["remove_passes"] = passes
    log.info(f"Removed labels of {removed} unsupported contigs in {passes} pass(es) "
             f"({time.time() - t_stage:.1f} s)")

    # -- 2. refine labels of inconsistent vertices -------------------------
    t_stage = time.time()
    once_moved = set()
    moved = 0
    passes = 0
    while "refine" in stages:
        passes += 1
        moves = []
        for v in range(n):
            my_bin = bin_of[v]
            if my_bin < 0 or not adj[v] or v in once_moved:
                continue
            found = run_bfs(v)
            if not found:
                continue
            weights = {}
            for d, b in found.values():
                weights[b] = weights.get(b, 0) + 1 / (2 ** d)
            limit = weights.get(my_bin, 0) * threshold
            max_weight = -1
            max_bin = -1
            for b in sorted(weights):
                w = weights[b]
                if b != my_bin and w > 0 and w > limit and max_weight < w:
                    max_weight = w
                    max_bin = b
            if max_bin != -1:
                moves.append((v, max_bin))
                once_moved.add(v)
        log.debug(f"Refine pass {passes}: {len(moves)} contigs moved")
        if not moves:
            break
        for v, b in moves:
            bin_of[v] = b
        moved += len(moves)
    info["moved"] = moved
    info["refine_passes"] = passes
    log.info(f"Moved {moved} inconsistent contigs in {passes} pass(es) "
             f"({time.time() - t_stage:.1f} s)")

    # -- 3. components with at least one labelled contig ("non-isolated") --
    t_stage = time.time()
    comp = [-1] * n
    n_comp = 0
    for s in range(n):
        if comp[s] != -1:
            continue
        comp[s] = n_comp
        stack = [s]
        while stack:
            v = stack.pop()
            for w in adj[v]:
                if comp[w] == -1:
                    comp[w] = n_comp
                    stack.append(w)
        n_comp += 1
    labelled_comp = set(comp[v] for v in range(n) if bin_of[v] >= 0)
    non_isolated = bytearray(1 if comp[v] in labelled_comp else 0 for v in range(n))
    to_label = bytearray(1 if (non_isolated[v] and bin_of[v] < 0) else 0 for v in range(n))
    info["components"] = n_comp
    info["non_isolated"] = sum(non_isolated)
    info["non_isolated_unbinned"] = sum(to_label)
    log.info(f"Connected components (including single contigs): {n_comp}")
    log.info(f"Number of non-isolated contigs: {info['non_isolated']}")
    log.info(f"Number of non-isolated unbinned contigs: {info['non_isolated_unbinned']}")

    # -- 4. propagate labels (versioned lazy invalidation) -----------------
    heap = []
    version = {}

    def candidates(x, ver):
        cx = coverages[x]
        return [(d, abs(cx - coverages[a]), x, b, a, ver)
                for a, (d, b) in run_bfs(x).items()]

    frontier = set()
    for v in range(n):
        if bin_of[v] >= 0 and non_isolated[v]:
            for w in adj[v]:
                if bin_of[w] < 0:
                    frontier.add(w)
    for x in (sorted(frontier) if "propagate" in stages else ()):
        heap.extend(candidates(x, 0))
    heapq.heapify(heap)
    info["initial_heap"] = len(heap)
    propagated = 0
    stale = 0
    heappop, heappush = heapq.heappop, heapq.heappush
    while heap:
        d, cd, x, b, a, ver = heappop(heap)
        if not to_label[x]:
            continue
        if ver != version.get(x, 0):
            stale += 1
            continue
        bin_of[x] = b
        to_label[x] = 0
        propagated += 1
        fresh = [w for w in adj[x] if bin_of[w] < 0]
        for w in fresh:
            version[w] = version.get(w, 0) + 1
        for w in fresh:
            for entry in candidates(w, version[w]):
                heappush(heap, entry)
    info["propagated"] = propagated
    info["stale_entries_skipped"] = stale
    log.info(f"Propagated labels to {propagated} contigs ({time.time() - t_stage:.1f} s)")

    # -- 5. multi-binned contigs ------------------------------------------
    t_stage = time.time()
    bin_cov_sum = [0] * n_bins
    bin_len_total = [0] * n_bins
    for v in range(n):
        b = bin_of[v]
        if b >= 0 and non_isolated[v]:
            bin_cov_sum[b] += coverages[v] * lengths[v]
            bin_len_total[b] += lengths[v]
    multi = []
    for v in (range(n) if "multibin" in stages else ()):
        cb = bin_of[v]
        if cb < 0 or not non_isolated[v]:
            continue
        nz = sorted(set(bin_of[w] for w in adj[v] if bin_of[w] >= 0))
        if len(nz) < 2:  # upstream: zero_bin_count < n_bins - 1
            continue
        if len(nz) > 20:
            log.warning(f"contig node {v} has neighbours in {len(nz)} bins; "
                        f"evaluating 2^{len(nz)} bin combinations")
        cov_v = coverages[v]
        len_v = lengths[v]
        mean = {}
        for i in nz:
            s = bin_cov_sum[i]
            ln = bin_len_total[i]
            if i == cb:
                s = s - cov_v * len_v
                ln = ln - len_v
            mean[i] = s / ln if ln != 0 else s
        min_diff = sys.maxsize
        best = None
        for r in range(1, len(nz) + 1):
            for combination in it.combinations(nz, r):
                total = 0
                for i in combination:
                    total += mean[i]
                diff = abs(total - cov_v)
                if diff < min_diff:
                    min_diff = diff
                    best = combination
        if best is not None and len(best) > 1 and len_v > 1000:
            multi.append((v, best))
    memberships = [((bin_of[v],) if bin_of[v] >= 0 else ()) for v in range(n)]
    for v, combination in multi:
        memberships[v] = tuple(sorted(set(memberships[v]) | set(combination)))
    info["multi_binned"] = len(multi)
    info["multi"] = multi
    log.info(f"Multi-bin stage: {len(multi)} contigs assigned to several bins "
             f"({time.time() - t_stage:.1f} s)")
    return memberships, info


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def _setup_logger(log_path, quiet=False):
    logger = logging.getLogger(f"graphbin2_fast.{id(log_path)}.{time.time()}")
    logger.setLevel(logging.DEBUG)
    logger.propagate = False
    fmt = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")
    if not quiet:
        ch = logging.StreamHandler()
        ch.setFormatter(fmt)
        ch.setLevel(logging.INFO)
        logger.addHandler(ch)
    fh = logging.FileHandler(log_path)
    fh.setFormatter(fmt)
    fh.setLevel(logging.DEBUG)
    logger.addHandler(fh)
    return logger


def _close_logger(logger):
    for h in list(logger.handlers):
        h.close()
        logger.removeHandler(h)


def run(contigs, graph, binned, output, prefix="", depth=5, threshold=1.5,
        delimiter=",", nthreads=1, fasta_all_bins=False, quiet=False):
    """Run GraphBin2 (MEGAHIT mode). Returns a dict of summary counts."""
    if prefix and not prefix.endswith("_"):
        prefix = prefix + "_"
    prefix = prefix or ""
    if not output.endswith("/"):
        output = output + "/"
    os.makedirs(output, exist_ok=True)
    logger = _setup_logger(f"{output}{prefix}graphbin2.log", quiet=quiet)
    try:
        return _run(contigs, graph, binned, output, prefix, depth, threshold,
                    delimiter, nthreads, fasta_all_bins, logger)
    except GraphBin2Error as err:
        logger.error(str(err))
        raise
    finally:
        _close_logger(logger)


def _run(contigs_file, graph_file, bins_file, output, prefix, depth, threshold,
         delimiter, nthreads, fasta_all_bins, logger):
    start = time.time()
    logger.info(f"graphbin2_fast {__version__}: corrected re-implementation of "
                f"GraphBin2 {UPSTREAM_VERSION} (MEGAHIT mode)")
    logger.info(f"Contigs file: {contigs_file}")
    logger.info(f"Assembly graph file: {graph_file}")
    logger.info(f"Existing binning output file: {bins_file}")
    logger.info(f"Final binning output folder: {output}")
    logger.info(f"Depth: {depth}  Threshold: {threshold}  Delimiter: {delimiter!r}")
    if nthreads and nthreads > 1:
        logger.info(f"Number of threads: {nthreads} (accepted for compatibility; "
                    "this implementation is single-threaded)")

    t = time.time()
    c = load_contigs(contigs_file)
    logger.info(f"Contigs read: {len(c.names)} ({time.time() - t:.1f} s)")
    ks = Counter()
    for name in c.names:
        m = K_PREFIX_RE.match(name)
        ks[int(m.group(1)) if m else None] += 1
    if ks and None not in ks:
        logger.info("Contig names indicate MEGAHIT k: " + ", ".join(
            f"k{k}: {v}" for k, v in sorted(ks.items())))
    if c.len_from_seq:
        logger.warning(f"{c.len_from_seq} contig headers had no 'len=' field; "
                       "the sequence length was used")

    t = time.time()
    g = load_gfa(graph_file)
    n_edges = g.n_edges
    n_isolated = sum(1 for a in g.adj if not a)
    logger.info(f"Total number of contigs available: {g.n}")
    logger.info(f"Total number of edges in the assembly graph: "
                f"{g.n_links - g.n_self_loops} (links without self-loops, as upstream)")
    logger.info(f"Unique undirected edges after simplification: {n_edges}; "
                f"self-loop links: {g.n_self_loops}; isolated contigs: {n_isolated} "
                f"({time.time() - t:.1f} s)")
    if g.n_links_unknown:
        logger.warning(f"{g.n_links_unknown} links refer to segments missing from the "
                       "graph and were ignored")
    if n_edges == 0:
        k_note = ""
        if ks and None not in ks and len(ks) == 1:
            k_note = (f" The contig names indicate a final MEGAHIT k of {next(iter(ks))};"
                      " use exactly this k for 'megahit_toolkit contig2fastg' and "
                      "'agtools fastg2gfa -k'.")
        logger.warning("The assembly graph has 0 edges: every contig is isolated and "
                       "GraphBin2 cannot change the binning (output = input bins). This "
                       "usually means the graph was built with a k-mer size different "
                       "from the assembly's final k." + k_note)

    node_contig, contig_node, n_pos, n_dig = map_segments_to_contigs(g, c)
    if n_dig:
        logger.info(f"Segments matched to contigs by position: {n_pos}; by sequence: {n_dig}")
    n_missing = sum(1 for x in contig_node if x < 0)
    if n_missing:
        logger.warning(f"{n_missing} contigs of the contigs file are not in the graph; "
                       "they cannot be refined and are ignored")
    names = [c.names[node_contig[v]] for v in range(g.n)]
    coverages = [c.coverages[node_contig[v]] for v in range(g.n)]
    lengths = [c.lengths[node_contig[v]] for v in range(g.n)]
    name_to_node = {names[v]: v for v in range(g.n)}

    first_bin, bstats, absent = load_bins(bins_file, delimiter, name_to_node)
    bins_list = sorted(set(first_bin.values()))
    bin_index = {b: i for i, b in enumerate(bins_list)}
    n_bins = len(bins_list)
    bin_of = [-1] * g.n
    for v, label in first_bin.items():
        bin_of[v] = bin_index[label]
    logger.info(f"Number of bins available in binning result: {n_bins}")
    logger.info("Initial bins file: " + ", ".join(f"{k}={v}" for k, v in sorted(bstats.items())))
    if bstats.get("contigs_absent_from_graph"):
        logger.warning(f"{bstats['contigs_absent_from_graph']} contigs of the bins file are "
                       f"not in the contigs file/graph and were ignored (e.g. {', '.join(absent)})")
    if bstats.get("contigs_in_several_bins_kept_first"):
        logger.warning(f"{bstats['contigs_in_several_bins_kept_first']} rows assign an "
                       "already binned contig to another bin; the first bin was kept")
    initial = list(bin_of)
    n_binned = sum(1 for b in bin_of if b >= 0)
    logger.info(f"Number of binned contigs: {n_binned}")
    logger.info(f"Total number of unbinned contigs: {g.n - n_binned}")
    logger.info(f"Number of isolated contigs: {n_isolated}")

    memberships, info = refine(g.adj, coverages, lengths, bin_of, n_bins,
                               depth=depth, threshold=threshold, log=logger)

    if not info["multi"]:
        logger.info("No multi-labelled contigs were found")
    else:
        logger.info(f"Found {len(info['multi'])} multi-labelled contigs")
        for v, combination in info["multi"]:
            logger.info(names[v] + " belongs to bins "
                        + ", ".join(bins_list[s] for s in combination))

    # -- outputs ------------------------------------------------------------
    t = time.time()
    out_csv = f"{output}{prefix}graphbin2_output.csv"
    with open(out_csv, mode="w") as fh:
        w = csv.writer(fh, delimiter=delimiter, quotechar='"', quoting=csv.QUOTE_MINIMAL)
        for v in range(g.n):
            for k in memberships[v]:
                w.writerow([names[v], bins_list[k]])

    bins_dir = f"{output}{prefix}bins/"
    os.makedirs(bins_dir, exist_ok=True)
    targets = {}
    for v in range(g.n):
        if memberships[v]:
            targets[names[v]] = memberships[v] if fasta_all_bins else (memberships[v][-1],)
    chunk = 500
    for first in range(0, n_bins, chunk):
        sel = range(first, min(n_bins, first + chunk))
        handles = {k: open(f"{bins_dir}{prefix}{bins_list[k]}.fasta", "w") for k in sel}
        try:
            for label, seq in read_fasta(contigs_file):
                ks_v = targets.get(label.split()[0])
                if ks_v:
                    for k in ks_v:
                        h = handles.get(k)
                        if h is not None:
                            h.write(f">{label}\n{seq}\n")
        finally:
            for h in handles.values():
                h.close()

    final_binned = sum(1 for m in memberships if m)
    changed = sum(1 for v in range(g.n)
                  if memberships[v] != ((initial[v],) if initial[v] >= 0 else ()))
    summary = {
        "version": __version__, "contigs": len(c.names), "nodes": g.n,
        "links_without_self_loops": g.n_links - g.n_self_loops, "edges": n_edges,
        "isolated": n_isolated, "components": info["components"], "bins": n_bins,
        "binned_initial": n_binned,
        "removed": info["removed"], "remove_passes": info["remove_passes"],
        "moved": info["moved"], "refine_passes": info["refine_passes"],
        "non_isolated": info["non_isolated"],
        "non_isolated_unbinned": info["non_isolated_unbinned"],
        "propagated": info["propagated"], "multi_binned": info["multi_binned"],
        "binned_final": final_binned, "contigs_changed": changed,
        "bins_file": dict(bstats), "seconds": round(time.time() - start, 2),
    }
    logger.info(f"Outputs written ({time.time() - t:.1f} s)")
    logger.info("Summary: " + ", ".join(f"{k}={v}" for k, v in summary.items()
                                       if k not in ("bins_file", "version")))
    logger.info(f"Elapsed time: {time.time() - start:.1f} seconds")
    logger.info(f"Final binning results can be found at {out_csv}")
    return summary


def build_parser():
    ap = argparse.ArgumentParser(
        prog="graphbin2_fast.py",
        description="Corrected, faster GraphBin2 1.3.3 (MEGAHIT mode). Refines an initial "
                    "binning with the assembly graph (megahit_toolkit contig2fastg + "
                    "agtools fastg2gfa, built with the assembly's final k).")
    ap.add_argument("--contigs", required=True, help="MEGAHIT contigs FASTA (original headers)")
    ap.add_argument("--graph", required=True, help="assembly graph in GFA format")
    ap.add_argument("--binned", required=True,
                    help="initial binning CSV (contig,bin), or a folder of bin FASTA files "
                         "(.fasta/.fa/.fna; read like graphbin2 prepResult)")
    ap.add_argument("--output", required=True, help="output folder")
    ap.add_argument("--prefix", default="", help="prefix for the output files")
    ap.add_argument("--depth", type=int, default=5,
                    help="maximum depth for the breadth-first search [5]")
    ap.add_argument("--threshold", type=float, default=1.5,
                    help="threshold for determining inconsistent vertices [1.5]")
    ap.add_argument("--delimiter", default=",",
                    help="delimiter of input and output CSV: , ; | tab space [,]")
    ap.add_argument("--nthreads", type=int, default=1,
                    help="accepted for compatibility; the run is single-threaded")
    ap.add_argument("--fasta-all-bins", action="store_true",
                    help="write multi-binned contigs into every bin FASTA they belong to "
                         "(upstream writes them only to the last bin)")
    ap.add_argument("--assembler", default="megahit",
                    help="compatibility option; only 'megahit' is supported")
    ap.add_argument("--abundance", default=None, help="compatibility option; ignored")
    ap.add_argument("--paths", default=None, help="compatibility option; ignored")
    ap.add_argument("--quiet", action="store_true", help="log to file only")
    ap.add_argument("-v", "--version", action="version", version=f"%(prog)s {__version__}")
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.assembler.lower() != "megahit":
        print("graphbin2_fast.py supports only --assembler megahit", file=sys.stderr)
        return 2
    if a.depth < 1:
        print("Please enter a valid number for depth (>= 1)", file=sys.stderr)
        return 2
    if a.threshold < 1.0:
        print("Please enter a valid number for threshold (>= 1.0)", file=sys.stderr)
        return 2
    if a.nthreads <= 0:
        print("Please enter a valid number for the number of threads", file=sys.stderr)
        return 2
    delim = DELIMITERS.get(a.delimiter)
    if delim is None:
        print(f"Unsupported delimiter {a.delimiter!r}", file=sys.stderr)
        return 2
    try:
        run(a.contigs, a.graph, a.binned, a.output, prefix=a.prefix, depth=a.depth,
            threshold=a.threshold, delimiter=delim, nthreads=a.nthreads,
            fasta_all_bins=a.fasta_all_bins, quiet=a.quiet)
    except GraphBin2Error as err:
        print(f"ERROR: {err}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
