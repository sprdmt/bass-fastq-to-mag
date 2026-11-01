#!/usr/bin/env python3
"""Per-MAG (and per-assembly) sequencing depth from CoverM contig tables.

Input tables
------------
``coverm contig --methods mean`` output, one table per assembly group:
``Contig<TAB><sample> Mean<TAB>...``. Headerless tables (for example the old
GraphBin2 ``abundance.tsv`` with the header line removed) are detected and
their columns named sample_1..n. Only columns ending in `` Mean`` are used
when a header is present.

Definitions (documented in docs/ena_submission.md)
--------------------------------------------------
For contig c with length L_c and per-sample CoverM mean depth d_cs:

* ``mean_depth_<sample>`` = sum_c(L_c * d_cs) / sum_c(L_c) over the MAG's
  contigs (length-weighted mean depth of the genome in that sample).
* ``genome_coverage`` (ENA COVERAGE, default ``--method sum_mean``, the
  definition used for the earlier ENA submissions) = mean over the MAG's
  contigs of the per-contig sum over samples, sum_s d_cs. It is the depth the
  co-assembly had from all reads pooled; each contig counts once regardless
  of its length.
* ``--method length_weighted`` instead reports sum_s mean_depth_s (length-
  weighted). The two agree when depth is uniform along contigs of all lengths.

The primary-assembly coverage (``assembly`` sub-command, used by
ena_coassembly.py) applies the same definitions to every contig in the table.

Usage
-----
  mag_coverage.py --mags-table results/mags/mags.tsv --mags-dir results/mags/fasta \
      --coverm g1=results/mapping/g1/coverm_contig_mean.tsv [g2=...] --out results/mags/coverage.tsv
  mag_coverage.py assembly --coverm results/mapping/g1/coverm_contig_mean.tsv [--assembly-fasta F]
"""
from __future__ import annotations

import argparse
import json
import logging
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import (EnaError, fasta_lengths, fmt_float, genome_id_from_path,  # noqa: E402
                        list_fasta_files, open_text, read_tsv, setup_logging, write_tsv)

LOG = logging.getLogger("bfm.mag_coverage")
METHODS = ("sum_mean", "length_weighted")
_MEAN_COL = re.compile(r"^(?P<sample>.+?)\s+Mean$")


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_header(path) -> Tuple[bool, List[str], List[int]]:
    """Return (has_header, sample_names, column_indices) for a CoverM contig table."""
    with open_text(path) as fh:
        first = fh.readline().rstrip("\r\n")
    if not first:
        raise EnaError(f"{path}: empty coverage table")
    cols = first.split("\t")
    if len(cols) < 2:
        raise EnaError(f"{path}: expected tab-separated contig and depth columns")
    if all(_is_number(c) for c in cols[1:]):
        return False, [f"sample_{i}" for i in range(1, len(cols))], list(range(1, len(cols)))
    names, idx = [], []
    for i, c in enumerate(cols[1:], start=1):
        m = _MEAN_COL.match(c.strip())
        if m:
            names.append(m.group("sample"))
            idx.append(i)
    if not idx:
        raise EnaError(f"{path}: no '<sample> Mean' columns in header {cols[:4]}...; "
                       "run 'coverm contig --methods mean'")
    return True, names, idx


def iter_depths(path):
    """Yield (contig, [depth per sample]) for every data row."""
    has_header, _names, idx = read_header(path)
    with open_text(path) as fh:
        if has_header:
            fh.readline()
        for n, line in enumerate(fh, start=2 if has_header else 1):
            line = line.rstrip("\r\n")
            if not line:
                continue
            cols = line.split("\t")
            try:
                yield cols[0], [float(cols[i]) for i in idx]
            except (IndexError, ValueError) as exc:
                raise EnaError(f"{path}:{n}: malformed row ({exc})") from exc


def load_depths(path, contigs: Optional[set] = None) -> Tuple[List[str], Dict[str, List[float]]]:
    """Return sample names and depths for ``contigs`` (all contigs if None)."""
    _h, names, _i = read_header(path)
    depths = {}
    for contig, vals in iter_depths(path):
        if contigs is None or contig in contigs:
            depths[contig] = vals
    return names, depths


def genome_depth(contig_lengths: Dict[str, int], depths: Dict[str, List[float]], n_samples: int,
                 method: str = "sum_mean") -> Dict[str, object]:
    """Per-sample length-weighted depth and ENA genome coverage for one genome."""
    found = [c for c in contig_lengths if c in depths]
    if not found:
        return {"n_found": 0, "per_sample": [math.nan] * n_samples, "genome_coverage": math.nan}
    total_len = sum(contig_lengths[c] for c in found)
    per_sample = [sum(contig_lengths[c] * depths[c][s] for c in found) / total_len for s in range(n_samples)]
    if method == "sum_mean":
        cov = sum(sum(depths[c]) for c in found) / len(found)
    elif method == "length_weighted":
        cov = sum(per_sample)
    else:
        raise EnaError(f"unknown method {method}")
    return {"n_found": len(found), "per_sample": per_sample, "genome_coverage": cov}


def assembly_coverage(path, method: str = "sum_mean", assembly_fasta=None) -> Dict[str, object]:
    """ENA COVERAGE of a whole (co-)assembly from its CoverM contig table (streamed)."""
    lengths = fasta_lengths(assembly_fasta) if (method == "length_weighted") else None
    if method == "length_weighted" and not assembly_fasta:
        raise EnaError("length_weighted assembly coverage needs --assembly-fasta for contig lengths")
    n = 0
    total = 0.0
    wsum = 0.0
    wlen = 0
    missing_len = 0
    _h, names, _i = read_header(path)
    for contig, vals in iter_depths(path):
        s = sum(vals)
        n += 1
        total += s
        if lengths is not None:
            L = lengths.get(contig)
            if L is None:
                missing_len += 1
                continue
            wsum += L * s
            wlen += L
    if n == 0:
        raise EnaError(f"{path}: no contigs")
    cov = total / n if method == "sum_mean" else (wsum / wlen if wlen else math.nan)
    return {"coverage": cov, "n_contigs": n, "n_samples": len(names), "method": method,
            "contigs_without_length": missing_len}


def parse_coverm_args(values: Sequence[str]) -> Dict[str, str]:
    """['g1=path', 'g2=path'] -> {g1: path}; a bare path maps to group ''."""
    out: Dict[str, str] = {}
    for v in values:
        if "=" in v and not Path(v).exists():
            g, p = v.split("=", 1)
            out[g] = p
        else:
            out[""] = v
    return out


def mag_coverage_table(mags: List[Dict[str, str]], mags_dir, coverm: Dict[str, str], method: str = "sum_mean",
                       min_found_fraction: float = 0.95, fasta_ext: Optional[str] = None) -> Tuple[List[Dict], List[str]]:
    """Compute coverage rows for MAGs given as dicts with mag_id and group."""
    by_group: Dict[str, List[Dict[str, str]]] = {}
    for m in mags:
        g = m.get("group", "") or ""
        if g not in coverm:
            if "" in coverm and len(coverm) == 1:
                g = ""
            else:
                raise EnaError(f"no --coverm table for group {g!r} (MAG {m['mag_id']})")
        by_group.setdefault(g, []).append(m)
    fasta_index = {genome_id_from_path(p): p for p in list_fasta_files(mags_dir, [fasta_ext] if fasta_ext else None)}
    rows, all_samples, problems = [], [], []
    for g, members in by_group.items():
        lengths = {}
        for m in members:
            p = fasta_index.get(m["mag_id"])
            if p is None:
                raise EnaError(f"FASTA for {m['mag_id']} not found in {mags_dir}")
            lengths[m["mag_id"]] = fasta_lengths(p)
        wanted = {c for L in lengths.values() for c in L}
        names, depths = load_depths(coverm[g], wanted)
        for s in names:
            if s not in all_samples:
                all_samples.append(s)
        for m in members:
            L = lengths[m["mag_id"]]
            res = genome_depth(L, depths, len(names), method)
            frac = res["n_found"] / len(L) if L else 0.0
            if frac < min_found_fraction:
                problems.append(f"{m['mag_id']}: {res['n_found']}/{len(L)} contigs found in {coverm[g]}")
            row = {"mag_id": m["mag_id"], "group": m.get("group", g), "n_contigs": len(L),
                   "n_contigs_with_depth": res["n_found"],
                   "genome_coverage": "" if math.isnan(res["genome_coverage"]) else fmt_float(res["genome_coverage"], 4),
                   "coverage_method": method}
            for s, v in zip(names, res["per_sample"]):
                row[f"mean_depth_{s}"] = "" if math.isnan(v) else fmt_float(v, 4)
            rows.append(row)
    if problems:
        raise EnaError("MAG contigs missing from the CoverM table (wrong group or renamed contigs?):\n  "
                       + "\n  ".join(problems[:20]))
    return rows, all_samples


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("mode", nargs="?", choices=("mags", "assembly"), default="mags")
    ap.add_argument("--coverm", nargs="+", required=True, help="<group>=<coverm contig table> (or one bare path)")
    ap.add_argument("--mags-table", help="workflow mags.tsv (mag_id, group); default: every FASTA in --mags-dir")
    ap.add_argument("--mags-dir", help="directory with <mag_id>.fna")
    ap.add_argument("--fasta-ext", default=None, help="restrict --mags-dir to this extension (e.g. .fna)")
    ap.add_argument("--method", choices=METHODS, default="sum_mean")
    ap.add_argument("--min-found-fraction", type=float, default=0.95,
                    help="fail if fewer MAG contigs than this fraction are in the table (default 0.95)")
    ap.add_argument("--assembly-fasta", help="assembly FASTA (lengths for length_weighted assembly coverage)")
    ap.add_argument("--out", help="output TSV (mags) or JSON (assembly; default stdout)")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        coverm = parse_coverm_args(args.coverm)
        if args.mode == "assembly":
            if len(coverm) != 1:
                raise EnaError("assembly mode takes exactly one --coverm table")
            res = assembly_coverage(next(iter(coverm.values())), args.method, args.assembly_fasta)
            text = json.dumps(res, indent=2)
            if args.out:
                Path(args.out).write_text(text + "\n")
            else:
                print(text)
            return 0
        if not args.mags_dir or not args.out:
            raise EnaError("mags mode needs --mags-dir and --out")
        if args.mags_table:
            mags = read_tsv(args.mags_table, required=("mag_id",))
            if mags and "group" not in mags[0]:
                LOG.warning("mags table has no 'group' column; all MAGs use the single --coverm table")
        else:
            mags = [{"mag_id": genome_id_from_path(p), "group": next(iter(coverm))}
                    for p in list_fasta_files(args.mags_dir, [args.fasta_ext] if args.fasta_ext else None)]
        rows, samples = mag_coverage_table(mags, args.mags_dir, coverm, args.method,
                                           args.min_found_fraction, args.fasta_ext)
    except EnaError as exc:
        LOG.error("%s", exc)
        return 1
    cols = ["mag_id", "group", "n_contigs", "n_contigs_with_depth", "genome_coverage", "coverage_method"] + \
        [f"mean_depth_{s}" for s in samples]
    write_tsv(args.out, rows, cols)
    LOG.info("wrote %d MAG coverage rows (%d samples) to %s", len(rows), len(samples), args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
