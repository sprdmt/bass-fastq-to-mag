"""Tests for workflow/scripts/graphbin2_fast.py.

Run:  pytest -q tests/graphbin2_fast
Comparisons with the reference need python-igraph (and optionally cogent3);
comparisons with the previously used patch (not shipped) and with unmodified upstream need the
graphbin2 1.3.3 package. Missing dependencies skip those tests.
"""
from __future__ import annotations

import ast
import random
import shutil
import statistics
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gb2_testlib as T  # noqa: E402

fast = T.load_fast()
FIXTURE = T.TESTS / "fixtures" / "medea2_subgraph"
N_RANDOM = 150


@pytest.fixture(scope="module")
def ref():
    pytest.importorskip("igraph")
    return T.load_reference()


@pytest.fixture(scope="module")
def ref_total():
    pytest.importorskip("igraph")
    return T.load_reference(total_order=True)


@pytest.fixture(scope="module")
def patch():
    pytest.importorskip("igraph")
    mod = T.load_patch()
    if mod is None:
        pytest.skip("previous patch or graphbin2 package not available")
    return mod


@pytest.fixture(scope="module")
def upstream():
    pytest.importorskip("igraph")
    mod = T.load_upstream()
    if mod is None:
        pytest.skip("unmodified graphbin2 1.3.3 not installed")
    return mod


def chain_instance():
    # 0-1-2-3-4 chain, only contig 0 binned; contig 5 isolated
    return T.Instance([f"k125_{i}" for i in (11, 7, 99, 3, 42, 5)], [10, 11, 12, 13, 14, 15],
                      [5000] * 6, [(0, 1), (1, 2), (2, 3), (3, 4)], {0: "bin.A.fa"})


def multibin_instance():
    # c (bin A, cov 30) touches l1 (A, cov 10) and l2 (B, cov 20); l3 (B, cov 20) supports l2.
    # Mean coverages A\{c} = 10 and B = 20 sum to 30 -> c belongs to A and B.
    names = ["k125_100", "k125_200", "k125_300", "k125_400", "k125_500"]
    return T.Instance(names, [30, 10, 20, 20, 3], [5000, 5000, 5000, 5000, 5000],
                      [(0, 1), (0, 2), (2, 3)],
                      {0: "bin.A.fa", 1: "bin.A.fa", 2: "bin.B.fa", 3: "bin.B.fa"})


# ---------------------------------------------------------------------------
# static properties
# ---------------------------------------------------------------------------

STDLIB = {"argparse", "csv", "gzip", "hashlib", "heapq", "itertools", "logging", "os", "re",
          "sys", "time", "collections", "__future__"}


def test_stdlib_only_and_python39_syntax():
    src = T.FAST_PATH.read_text()
    tree = ast.parse(src, feature_version=(3, 9))
    mods = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name.split(".")[0] for a in node.names)
        elif isinstance(node, ast.ImportFrom):
            mods.add((node.module or "").split(".")[0])
    assert mods <= STDLIB, mods - STDLIB


def test_licence_and_modification_notice():
    src = T.FAST_PATH.read_text()
    assert "Copyright (c) 2022, Vijini Mallawaarachchi" in src
    assert "BSD 3-Clause" in src and "Modification notice" in src


def test_neighbour_order_matches_igraph(tmp_path):
    igraph = pytest.importorskip("igraph")
    for seed in range(40):
        inst = T.random_instance(seed)
        d = inst.write(tmp_path / f"g{seed}")
        g = fast.load_gfa(str(d / "graph.gfa"))
        ig = igraph.Graph()
        ig.add_vertices(inst.n)
        ig.add_edges(inst.edges + [(b, a) for a, b in inst.edges[::3]])
        ig.simplify(multiple=True, loops=False, combine_edges=None)
        assert [tuple(ig.neighbors(v, mode="all")) for v in range(inst.n)] == list(g.adj)


# ---------------------------------------------------------------------------
# targeted cases
# ---------------------------------------------------------------------------

def test_chain_propagates_beyond_one_hop(tmp_path, ref, patch):
    inst = chain_instance()
    d = inst.write(tmp_path / "in")
    expected = {n: {"bin.A.fa"} for n in inst.names[:5]}
    assert T.run_fast(fast, d, tmp_path / "fast") == expected
    assert T.run_module(ref, d, tmp_path / "ref") == expected
    # the previous patch marks unbinned neighbours as processed and then skips their
    # fresh candidates, so propagation stops after one hop
    assert T.run_module(patch, d, tmp_path / "patch") == {"k125_11": {"bin.A.fa"},
                                                          "k125_7": {"bin.A.fa"}}


def test_star_all_leaves_labelled(tmp_path, ref, patch):
    names = [f"k125_{i}" for i in range(1, 8)]
    inst = T.Instance(names, [20, 21, 19, 22, 18, 23, 17], [3000] * 7,
                      [(0, i) for i in range(1, 7)], {0: "bin.1.fa"})
    d = inst.write(tmp_path / "in")
    expected = {n: {"bin.1.fa"} for n in names}
    assert T.run_fast(fast, d, tmp_path / "fast") == expected
    assert T.run_module(ref, d, tmp_path / "ref") == expected
    assert T.run_module(patch, d, tmp_path / "patch") == expected  # one hop suffices here


def test_multibin_and_fasta_semantics(tmp_path, ref):
    inst = multibin_instance()
    d = inst.write(tmp_path / "in")
    rf = T.run_fast(fast, d, tmp_path / "fast")
    assert rf["k125_100"] == {"bin.A.fa", "bin.B.fa"}
    assert rf == T.run_module(ref, d, tmp_path / "ref")  # reference logs the name (fix c)
    log = (tmp_path / "fast" / "graphbin2.log").read_text()
    assert "k125_100 belongs to bins bin.A.fa, bin.B.fa" in log
    # upstream semantics: a multi-binned contig goes into the FASTA of its last bin only
    a = (tmp_path / "fast" / "bins" / "bin.A.fa.fasta").read_text()
    b = (tmp_path / "fast" / "bins" / "bin.B.fa.fasta").read_text()
    assert ">k125_100 " not in a and ">k125_100 " in b
    fast.run(str(d / "contigs.fa"), str(d / "graph.gfa"), str(d / "bins.csv"),
             str(tmp_path / "all"), fasta_all_bins=True, quiet=True)
    assert ">k125_100 " in (tmp_path / "all" / "bins" / "bin.A.fa.fasta").read_text()


def test_zero_edge_graph_is_noop_with_warning(tmp_path, ref):
    inst = T.random_instance(7, n=40)
    d = inst.write(tmp_path / "in", no_edges=True)
    rf = T.run_fast(fast, d, tmp_path / "fast")
    assert rf == T.initial_result(inst)
    assert T.run_module(ref, d, tmp_path / "ref") == rf
    log = (tmp_path / "fast" / "graphbin2.log").read_text()
    assert "WARNING" in log and "0 edges" in log and "final MEGAHIT k of 125" in log


def test_tolerant_bins_parsing(tmp_path):
    inst = T.random_instance(3, n=50)
    clean = inst.write(tmp_path / "clean")
    expected = T.run_fast(fast, clean, tmp_path / "out_clean")
    rows = []
    for k, v in enumerate(sorted(inst.bins)):
        full = f"{inst.names[v]} flag=0 multi={inst.coverages[v]}.0000 len={inst.lengths[v]}"
        rows.append(f"{full},{inst.bins[v]}\r\n")
        if k % 2 == 0:
            rows.append(f"{inst.names[v]},{inst.bins[v]}\r\n")       # duplicate row
    rows.insert(1, "\r\n")                                            # blank line
    rows.append("k125_absent_1,bin.9.fa\r\n")                         # contig not in graph
    first = min(inst.bins)
    rows.append(f"{inst.names[first]},bin.other.fa\r\n")              # conflicting second bin
    messy = tmp_path / "messy"
    shutil.copytree(clean, messy)
    with open(messy / "bins.csv", "w", newline="") as fh:
        fh.write("".join(rows))
    s = fast.run(str(messy / "contigs.fa"), str(messy / "graph.gfa"), str(messy / "bins.csv"),
                 str(tmp_path / "out_messy"), quiet=True)
    assert T.read_result(tmp_path / "out_messy" / "graphbin2_output.csv") == expected
    st = s["bins_file"]
    assert st.get("contig_names_trimmed", 0) == len(inst.bins)
    assert st.get("duplicate_rows", 0) == (len(inst.bins) + 1) // 2 > 0
    assert st.get("contigs_absent_from_graph", 0) == 1
    assert st.get("contigs_in_several_bins_kept_first", 0) == 1


def test_bins_folder_input_equals_prepresult_csv(tmp_path):
    """--binned DIR (bin FASTAs) gives the result of upstream prepResult's CSV."""
    import subprocess
    pytest.importorskip("graphbin2")
    inst = T.random_instance(9, n=40)
    d = inst.write(tmp_path / "in")
    bdir = tmp_path / "bin_fastas"
    bdir.mkdir()
    by_bin = {}
    for v, b in inst.bins.items():
        by_bin.setdefault(b, []).append(v)
    for b, vs in by_bin.items():
        with open(bdir / b, "w") as fh:
            for v in vs:
                fh.write(f">{inst.names[v]} flag=0 multi={inst.coverages[v]}.0000 "
                         f"len={inst.lengths[v]}\n{inst.seqs[v]}\n")
    subprocess.run([sys.executable, "-m", "graphbin2.support.prepResult", "--binned", str(bdir),
                    "--output", str(tmp_path / "prep")], check=True, capture_output=True)
    prep_csv = tmp_path / "prep" / "initial_contig_bins.csv"
    via_csv = fast.run(str(d / "contigs.fa"), str(d / "graph.gfa"), str(prep_csv),
                       str(tmp_path / "a"), quiet=True)
    via_dir = fast.run(str(d / "contigs.fa"), str(d / "graph.gfa"), str(bdir),
                       str(tmp_path / "b"), quiet=True)
    assert via_csv["bins_file"]["contig_names_trimmed"] == len(inst.bins)
    assert ((tmp_path / "a" / "graphbin2_output.csv").read_bytes()
            == (tmp_path / "b" / "graphbin2_output.csv").read_bytes())
    assert T.read_result(tmp_path / "a" / "graphbin2_output.csv") == T.run_fast(
        fast, d, tmp_path / "c")


def test_gfa_segment_order_differs_from_fasta(tmp_path):
    inst = T.random_instance(11, n=45)
    d = inst.write(tmp_path / "in")
    expected = T.run_fast(fast, d, tmp_path / "a")
    lines = (d / "graph.gfa").read_text().splitlines(keepends=True)
    s_lines = [x for x in lines if x.startswith("S")]
    random.Random(1).shuffle(s_lines)
    other = [x for x in lines if not x.startswith("S")]
    (d / "graph.gfa").write_text(other[0] + "".join(s_lines) + "".join(other[1:]))
    assert T.run_fast(fast, d, tmp_path / "b") == expected


def test_cli_naming_and_bytes_match_reference(tmp_path, ref):
    inst = T.random_instance(5, tie_free=True)
    d = inst.write(tmp_path / "in")
    rc = fast.main(["--contigs", str(d / "contigs.fa"), "--graph", str(d / "graph.gfa"),
                    "--binned", str(d / "bins.csv"), "--output", str(tmp_path / "fast"),
                    "--prefix", "mb2", "--abundance", str(d / "bins.csv"),
                    "--nthreads", "4", "--quiet"])
    assert rc == 0
    T.run_module(ref, d, tmp_path / "ref", prefix="mb2_")
    for sub in ("fast", "ref"):
        assert (tmp_path / sub / "mb2_graphbin2_output.csv").exists()
        assert (tmp_path / sub / "mb2_graphbin2.log").exists()
    fa = sorted(p.name for p in (tmp_path / "fast" / "mb2_bins").iterdir())
    fr = sorted(p.name for p in (tmp_path / "ref" / "mb2_bins").iterdir())
    assert fa == fr and all(x.startswith("mb2_") for x in fa)
    assert ((tmp_path / "fast" / "mb2_graphbin2_output.csv").read_bytes()
            == (tmp_path / "ref" / "mb2_graphbin2_output.csv").read_bytes())
    for x in fa:
        assert ((tmp_path / "fast" / "mb2_bins" / x).read_bytes()
                == (tmp_path / "ref" / "mb2_bins" / x).read_bytes())


def test_cli_rejects_invalid_arguments(tmp_path):
    base = ["--contigs", "x", "--graph", "y", "--binned", "z", "--output", str(tmp_path)]
    assert fast.main(base + ["--depth", "0"]) == 2
    assert fast.main(base + ["--threshold", "0.5"]) == 2
    assert fast.main(base + ["--assembler", "spades"]) == 2


# ---------------------------------------------------------------------------
# randomised equivalence
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("seed", range(N_RANDOM))
def test_random_tie_free_equals_reference(tmp_path, ref, seed):
    """Coverages with pairwise-distinct differences -> no heap ties -> identical."""
    inst = T.random_instance(1000 + seed, tie_free=True)
    d = inst.write(tmp_path / "in")
    assert T.run_fast(fast, d, tmp_path / "fast") == T.run_module(ref, d, tmp_path / "ref")


@pytest.mark.parametrize("seed", range(N_RANDOM))
def test_random_with_ties_equals_reference_total_order(tmp_path, ref_total, seed):
    """With ties, results equal the reference run with the same total heap order."""
    inst = T.random_instance(5000 + seed, tie_free=False)
    d = inst.write(tmp_path / "in")
    assert T.run_fast(fast, d, tmp_path / "fast") == T.run_module(ref_total, d, tmp_path / "ref")


def test_random_with_ties_agreement_with_upstream_heap_order(tmp_path, ref):
    """Upstream breaks heap ties by heap position; agreement is high but not total."""
    agree = []
    for seed in range(N_RANDOM):
        inst = T.random_instance(5000 + seed, tie_free=False)
        d = inst.write(tmp_path / f"in{seed}")
        rf = T.run_fast(fast, d, tmp_path / f"f{seed}")
        rr = T.run_module(ref, d, tmp_path / f"r{seed}")
        agree.append(T.agreement(rf, rr))
    exact = sum(1 for a in agree if a == 1.0)
    print(f"\n[ties] {exact}/{len(agree)} identical; mean contig agreement "
          f"{statistics.mean(agree):.4f}; min {min(agree):.3f}")
    assert exact >= 0.85 * len(agree)
    assert statistics.mean(agree) >= 0.97


def test_random_patch_differs_from_reference(tmp_path, ref, patch):
    same = 0
    for seed in range(60):
        inst = T.random_instance(1000 + seed, tie_free=True)
        d = inst.write(tmp_path / f"in{seed}")
        same += T.run_module(patch, d, tmp_path / f"p{seed}") == T.run_module(
            ref, d, tmp_path / f"r{seed}")
    print(f"\n[patch] identical to reference in {same}/60 tie-free instances")
    assert same < 45


# ---------------------------------------------------------------------------
# unmodified upstream 1.3.3: documented failures
# ---------------------------------------------------------------------------

def test_upstream_off_by_one_keyerror(tmp_path, upstream):
    inst = T.random_instance(20000, tie_free=True)  # 61 contigs; lookup of NODE n+1 fails
    d = inst.write(tmp_path / "in")
    with pytest.raises(KeyError):
        T.run_module(upstream, d, tmp_path / "up")


def test_upstream_contig_names_nameerror(tmp_path, upstream, ref):
    inst = T.random_instance(20014, tie_free=True)  # multi-labelled contig found
    d = inst.write(tmp_path / "in")
    with pytest.raises(NameError):
        T.run_module(upstream, d, tmp_path / "up")
    rr = T.run_module(ref, d, tmp_path / "ref")
    assert any(len(v) > 1 for v in rr.values())


# ---------------------------------------------------------------------------
# real medea-II subgraph (real headers, topology and MetaBAT2 bins)
# ---------------------------------------------------------------------------

@pytest.mark.skipif(not (FIXTURE / "graph.gfa").exists(), reason="fixture missing")
def test_medea2_subgraph(tmp_path, ref, ref_total, patch):
    rf = T.run_fast(fast, FIXTURE, tmp_path / "fast")
    assert rf == T.run_module(ref_total, FIXTURE, tmp_path / "ref_total")
    rr = T.run_module(ref, FIXTURE, tmp_path / "ref")
    rp = T.run_module(patch, FIXTURE, tmp_path / "patch")
    init = T.read_result(FIXTURE / "bins.csv")
    print(f"\n[medea2 fixture] changed vs initial: {sum(1 for k in set(rf) | set(init) if rf.get(k) != init.get(k))}; "
          f"agreement fast/reference {T.agreement(rf, rr):.4f}; patch/reference {T.agreement(rp, rr):.4f}")
    assert rf != init
    assert T.agreement(rf, rr) >= 0.95
