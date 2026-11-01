#!/usr/bin/env python3
"""Collect MAG candidates from all assembly groups for dereplication.

Each candidate bin is copied as <group>__<bin>.fa. Quality values come from Binette's
final_bins_quality_reports.tsv (automatic route; Binette 1.1 and 1.2 layouts), a CheckM2 quality_report.tsv, or a simple
table bin/completeness/contamination/origin (anvi'o route, written by anvio_candidates.py).
Writes:
  * genome_info.csv  - dRep --genomeInfo table (genome, completeness, contamination)
  * candidates.tsv   - provenance and basic statistics of every candidate
"""

import argparse
import csv
import os
import shutil
import sys

FASTA_EXTENSIONS = (".fasta", ".fna", ".fa")


def stem(filename):
    for ext in FASTA_EXTENSIONS:
        if filename.endswith(ext):
            return filename[: -len(ext)]
    return filename


def fasta_stats(path):
    lengths, current = [], 0
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if current:
                    lengths.append(current)
                current = 0
            else:
                current += len(line.strip())
    if current:
        lengths.append(current)
    total = sum(lengths)
    n50, acc = 0, 0
    for length in sorted(lengths, reverse=True):
        acc += length
        if acc >= total / 2:
            n50 = length
            break
    return total, len(lengths), n50


def binette_bin_name(row):
    """File stem of a Binette bin: 'bin_<bin_id>' (Binette 1.1) or the column 'name' (1.2)."""
    if row.get("bin_id"):
        return f"bin_{row['bin_id']}"
    return row["name"]


def read_quality(report, fmt):
    """Map bin file stem -> (completeness, contamination, origin)."""
    quality = {}
    with open(report, newline="") as fh:
        rows = csv.DictReader(fh, delimiter="\t")
        for row in rows:
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            if fmt == "binette":
                quality[binette_bin_name(row)] = (float(row["completeness"]), float(row["contamination"]),
                                                  row.get("origin", ""))
            elif fmt == "simple":
                quality[row["bin"]] = (float(row["completeness"]), float(row["contamination"]), row.get("origin", ""))
            else:
                quality[row["Name"]] = (float(row["Completeness"]), float(row["Contamination"]), "checkm2")
    return quality


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--groups", nargs="+", required=True)
    parser.add_argument("--bins", nargs="+", required=True, help="bin directory per group")
    parser.add_argument("--reports", nargs="+", required=True, help="quality report per group")
    parser.add_argument("--report-format", choices=["binette", "checkm2", "simple"], required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--genome-info", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args(argv)

    if not (len(args.groups) == len(args.bins) == len(args.reports)):
        sys.exit("--groups, --bins and --reports must have the same length")
    os.makedirs(args.out_dir, exist_ok=True)

    rows = []
    for group, bins_dir, report in zip(args.groups, args.bins, args.reports):
        quality = read_quality(report, args.report_format)
        files = sorted(f for f in os.listdir(bins_dir) if f.endswith(FASTA_EXTENSIONS))
        for filename in files:
            name = stem(filename)
            if name not in quality:
                print(f"WARNING: no quality values for {group}/{filename}; skipped", file=sys.stderr)
                continue
            completeness, contamination, origin = quality[name]
            candidate = f"{group}__{name}.fa"
            shutil.copyfile(os.path.join(bins_dir, filename), os.path.join(args.out_dir, candidate))
            size, n_contigs, n50 = fasta_stats(os.path.join(args.out_dir, candidate))
            rows.append(
                {
                    "candidate": candidate,
                    "group": group,
                    "source_bin": filename,
                    "origin": origin,
                    "completeness": completeness,
                    "contamination": contamination,
                    "size_bp": size,
                    "n_contigs": n_contigs,
                    "n50": n50,
                }
            )

    with open(args.genome_info, "w", newline="") as fh:
        writer = csv.writer(fh, lineterminator="\n")
        writer.writerow(["genome", "completeness", "contamination"])
        for row in rows:
            writer.writerow([row["candidate"], row["completeness"], row["contamination"]])
    fields = ["candidate", "group", "source_bin", "origin", "completeness", "contamination",
              "size_bp", "n_contigs", "n50"]
    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"{len(rows)} candidates from {len(args.groups)} group(s)")


if __name__ == "__main__":
    main()
