"""Helpers for the GraphBin2 tests: synthetic instances and runners.

Implementations compared:
  fast        workflow/scripts/graphbin2_fast.py (stdlib only)
  reference   tests/reference/graphbin2_MEGAHIT_reference.py = upstream 1.3.3
              with only the coverage-index fix (b) and the contig_names fix (c)
  ref_total   reference with the propagation heap ordered by the total key of
              graphbin2_fast (generated from the reference source at runtime)
  patch       the previously used patched module (not shipped; fixtures/...user_patch.py or
              $GB2_USER_PATCH); needs the graphbin2 package for its relative import
  upstream    unmodified graphbin2 1.3.3 from the installed package
"""
from __future__ import annotations

import csv
import importlib
import importlib.util
import logging
import os
import random
import sys
import types
from pathlib import Path

os.environ.setdefault("TQDM_DISABLE", "1")

TESTS = Path(__file__).resolve().parent
REPO = TESTS.parents[1]  # tests/graphbin2_fast -> repository root
FAST_PATH = REPO / "workflow" / "scripts" / "graphbin2_fast.py"
REF_PATH = TESTS / "reference" / "graphbin2_MEGAHIT_reference.py"
PATCH_PATH = Path(os.environ.get("GB2_USER_PATCH",
                                 TESTS / "fixtures" / "graphbin2_MEGAHIT_user_patch.py"))

ORIG_LT = "            return (self.data[3], self.data[-1]) < (other.data[3], other.data[-1])\n"
TOTAL_LT = ("            return (self.data[3], self.data[4], self.data[0], self.data[2], self.data[1]) "
            "< (other.data[3], other.data[4], other.data[0], other.data[2], other.data[1])\n")


def _load(name, path, source=None):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    if source is None:
        spec.loader.exec_module(mod)
    else:
        mod.__file__ = str(path)
        exec(compile(source, str(path), "exec"), mod.__dict__)
    sys.modules[name] = mod
    return mod


def load_fast():
    return _load("graphbin2_fast", FAST_PATH)


def load_reference(total_order=False):
    src = REF_PATH.read_text()
    if total_order:
        assert src.count(ORIG_LT) == 1
        src = src.replace(ORIG_LT, TOTAL_LT)
        return _load("graphbin2_ref_total", REF_PATH, src)
    return _load("graphbin2_ref", REF_PATH, src)


def load_patch():
    if not PATCH_PATH.exists():
        return None
    try:
        importlib.import_module("graphbin2.bidirectionalmap.bidirectionalmap")
    except ImportError:
        return None
    return _load("graphbin2._user_patch_under_test", PATCH_PATH)


def load_upstream():
    try:
        mod = importlib.import_module("graphbin2.graphbin2_MEGAHIT")
    except ImportError:
        return None
    return mod if getattr(mod, "__version__", "") == "1.3.3" else None


# ---------------------------------------------------------------------------
# instances
# ---------------------------------------------------------------------------

BASES = "ACGT"


class Instance:
    """A MEGAHIT-like test instance; node order = FASTA order = GFA order."""

    def __init__(self, names, coverages, lengths, edges, bins, seqs=None, multi_str=None):
        self.names = list(names)
        self.coverages = list(coverages)
        self.lengths = list(lengths)
        self.edges = list(edges)
        self.bins = dict(bins)            # node -> bin label
        self.multi_str = multi_str        # optional exact 'multi=' strings
        rnd = random.Random(len(names) * 7919 + len(edges))
        if seqs is None:
            seqs, seen = [], set()
            while len(seqs) < len(names):
                s = "".join(rnd.choice(BASES) for _ in range(24))
                if s not in seen:
                    seen.add(s)
                    seqs.append(s)
        self.seqs = seqs

    @property
    def n(self):
        return len(self.names)

    def write(self, d, bins_rows=None, gfa_extra_links=True, no_edges=False):
        d = Path(d)
        d.mkdir(parents=True, exist_ok=True)
        with open(d / "contigs.fa", "w") as fh:
            for i in range(self.n):
                multi = self.multi_str[i] if self.multi_str else f"{self.coverages[i]}.0000"
                fh.write(f">{self.names[i]} flag=0 multi={multi} len={self.lengths[i]}\n"
                         f"{self.seqs[i]}\n")
        seg = [f"NODE_{i + 1}_length_{self.lengths[i]}_cov_{self.coverages[i]}_ID_{2 * i + 1}"
               for i in range(self.n)]
        with open(d / "graph.gfa", "w") as fh:
            fh.write("H\tVN:Z:1.0\n")
            for i in range(self.n):
                fh.write(f"S\t{seg[i]}\t{self.seqs[i]}\n")
            if not no_edges:
                for k, (a, b) in enumerate(self.edges):
                    fh.write(f"L\t{seg[a]}\t+\t{seg[b]}\t-\t125M\n")
                    if gfa_extra_links and k % 3 == 0:  # reverse duplicate, as in real GFAs
                        fh.write(f"L\t{seg[b]}\t+\t{seg[a]}\t-\t125M\n")
                if gfa_extra_links and self.n:
                    fh.write(f"L\t{seg[0]}\t+\t{seg[0]}\t-\t125M\n")  # self-loop
        with open(d / "bins.csv", "w", newline="") as fh:
            if bins_rows is not None:
                fh.write(bins_rows)
            else:
                w = csv.writer(fh)
                for v in sorted(self.bins):
                    w.writerow([self.names[v], self.bins[v]])
        return d


def sidon_coverages(n, rnd, bits=40):
    """n distinct integers < 2**bits with pairwise-distinct differences."""
    while True:
        vals = rnd.sample(range(1, 2 ** bits), n)
        diffs = set()
        ok = True
        for i in range(n):
            for j in range(i + 1, n):
                dlt = abs(vals[i] - vals[j])
                if dlt in diffs:
                    ok = False
                    break
                diffs.add(dlt)
            if not ok:
                break
        if ok:
            return vals


def random_instance(seed, n=None, kind=None, tie_free=True):
    rnd = random.Random(seed)
    n = n or rnd.randint(6, 70)
    kind = kind or rnd.choice(["er", "planted", "tree", "chains", "grid"])
    k = rnd.randint(2, 5)
    comm = [rnd.randrange(k) for _ in range(n)]
    edges = set()
    if kind == "er":
        m = int(n * rnd.uniform(0.5, 1.6))
        for _ in range(m):
            a, b = rnd.randrange(n), rnd.randrange(n)
            if a != b:
                edges.add((min(a, b), max(a, b)))
    elif kind == "planted":
        for a in range(n):
            for b in range(a + 1, n):
                p = 0.25 if comm[a] == comm[b] else 0.02
                if rnd.random() < p * 6 / max(n, 6):
                    edges.add((a, b))
    elif kind == "tree":
        for v in range(1, n):
            u = rnd.randrange(v)
            edges.add((u, v))
    elif kind == "chains":
        order = list(range(n))
        rnd.shuffle(order)
        for i in range(n - 1):
            if rnd.random() < 0.85:
                a, b = order[i], order[i + 1]
                edges.add((min(a, b), max(a, b)))
    elif kind == "grid":
        w = max(2, int(n ** 0.5))
        for v in range(n):
            if (v + 1) % w and v + 1 < n:
                edges.add((v, v + 1))
            if v + w < n and rnd.random() < 0.7:
                edges.add((v, v + w))
    edges = sorted(edges)
    if tie_free:
        cov = sidon_coverages(n, rnd)
    else:
        base = [rnd.randint(2, 40) for _ in range(k)]
        cov = [max(1, base[comm[v]] + rnd.randint(-3, 3)) for v in range(n)]
    lengths = [rnd.choice([1000, 1500, 2500, 4000, 9000, 12000, 30000]) + rnd.randint(0, 50)
               if rnd.random() < 0.9 else 1000 for _ in range(n)]
    frac = rnd.uniform(0.15, 0.7)
    bins = {}
    for v in range(n):
        if rnd.random() < frac:
            b = comm[v] if rnd.random() > 0.2 else rnd.randrange(k)
            bins[v] = f"bin.{b + 1}.fa"
    names = rnd.sample(range(1, 10 ** 7), n)
    names = [f"k125_{x}" for x in names]
    return Instance(names, cov, lengths, edges, bins)


# ---------------------------------------------------------------------------
# runners
# ---------------------------------------------------------------------------

class Args:
    def __init__(self, d, out, prefix="", depth=5, threshold=1.5, delimiter=",", nthreads=1,
                 binned=None):
        d = Path(d)
        self.contigs = str(d / "contigs.fa")
        self.graph = str(d / "graph.gfa")
        self.binned = str(binned or d / "bins.csv")
        self.paths = None
        self.abundance = None
        self.output = str(out).rstrip("/") + "/"
        self.prefix = prefix
        self.depth = depth
        self.threshold = threshold
        self.delimiter = delimiter
        self.nthreads = nthreads


def _clear_loggers():
    for name in list(logging.Logger.manager.loggerDict):
        if name.startswith("GraphBin2"):
            lg = logging.getLogger(name)
            for h in list(lg.handlers):
                h.close()
                lg.removeHandler(h)


def run_module(mod, d, out, **kw):
    Path(out).mkdir(parents=True, exist_ok=True)
    args = Args(d, out, **kw)
    try:
        mod.run(args)
    finally:
        _clear_loggers()
    return read_result(Path(args.output) / f"{args.prefix}graphbin2_output.csv")


def run_fast(fast, d, out, **kw):
    a = Args(d, out, **kw)
    fast.run(a.contigs, a.graph, a.binned, a.output, prefix=a.prefix, depth=a.depth,
             threshold=a.threshold, delimiter=a.delimiter, quiet=True)
    prefix = a.prefix if (not a.prefix or a.prefix.endswith("_")) else a.prefix + "_"
    return read_result(Path(a.output) / f"{prefix}graphbin2_output.csv")


def read_result(path):
    res = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh):
            if row:
                res.setdefault(row[0], set()).add(row[1])
    return res


def initial_result(inst):
    return {inst.names[v]: {b} for v, b in inst.bins.items()}


def agreement(a, b):
    """Fraction of contigs (labelled in a or b) with identical label sets."""
    keys = set(a) | set(b)
    if not keys:
        return 1.0
    same = sum(1 for k in keys if a.get(k, set()) == b.get(k, set()))
    return same / len(keys)
