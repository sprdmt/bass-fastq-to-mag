#!/usr/bin/env python3
"""Combine kept Binette bins with the manually refined bins of one assembly group.

Kept bins (decision 'keep' in triage.tsv) are taken with their Binette/CheckM2 values and
Binette origin; bins of
the refined anvi'o collection are taken with the CheckM2 values computed after refinement and
prefixed 'refined_'. Output: a directory of candidate FASTA files and a quality table
(bin, completeness, contamination, origin).
"""

import argparse
import csv
import os
import shutil


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--triage", required=True)
    p.add_argument("--binette-bins", required=True)
    p.add_argument("--refined-bins", required=True)
    p.add_argument("--refined-checkm2", required=True, help="CheckM2 quality_report.tsv of the refined bins")
    p.add_argument("--out-dir", required=True)
    p.add_argument("--quality", required=True)
    a = p.parse_args(argv)

    os.makedirs(a.out_dir, exist_ok=True)
    rows = []
    for t in read_tsv(a.triage):
        if t["decision"] != "keep":
            continue
        shutil.copyfile(os.path.join(a.binette_bins, t["bin"] + ".fa"), os.path.join(a.out_dir, t["bin"] + ".fa"))
        rows.append({"bin": t["bin"], "completeness": t["completeness"],
                     "contamination": t["contamination"], "origin": t.get("origin") or "binette"})
    refined = {r["Name"]: r for r in read_tsv(a.refined_checkm2)} if os.path.getsize(a.refined_checkm2) else {}
    for fname in sorted(os.listdir(a.refined_bins)):
        name = os.path.splitext(fname)[0]
        if name not in refined:
            continue
        new = f"refined_{name}"
        shutil.copyfile(os.path.join(a.refined_bins, fname), os.path.join(a.out_dir, new + ".fa"))
        rows.append({"bin": new, "completeness": refined[name]["Completeness"],
                     "contamination": refined[name]["Contamination"], "origin": "anvio_refined"})
    with open(a.quality, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["bin", "completeness", "contamination", "origin"], delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    n_refined = sum(r["origin"] == "anvio_refined" for r in rows)
    print(f"{len(rows) - n_refined} kept + {n_refined} refined bins")


if __name__ == "__main__":
    main()
