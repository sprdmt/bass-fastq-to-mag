#!/usr/bin/env python3
"""Name the dereplicated MAGs and write the MAG table.

dRep winners that pass the quality thresholds are sorted by assembly group and source bin
(natural order) and renamed <project>_MAG_<n>. Contig names are kept unchanged, so every
MAG can be traced back to its co-assembly. Output: one <mag>.fna per MAG and mags.tsv.
"""

import argparse
import csv
import os
import re
import shutil
import sys
from collections import Counter


def natural_key(text):
    return [int(tok) if tok.isdigit() else tok.lower() for tok in re.split(r"(\d+)", text)]


def read_tsv(path, delimiter="\t"):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter=delimiter))


def quality_class(completeness, contamination):
    """Class from CheckM2 values only. 'near-complete' (>= 90% / <= 5%, Almeida et al. 2019)
    becomes MIMAG 'high-quality' in mag_summary.tsv when the 23S, 16S and 5S rRNA genes and
    >= 18 tRNAs are present (mimag_rna.py)."""
    if completeness >= 90 and contamination <= 5:
        return "near-complete"
    if completeness >= 50 and contamination < 10:
        return "medium-quality"
    return "low-quality"


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--winners", required=True, help="dRep data_tables/Wdb.csv")
    parser.add_argument("--clusters", required=True, help="dRep data_tables/Cdb.csv")
    parser.add_argument("--candidates", required=True, help="candidates.tsv")
    parser.add_argument("--genomes-dir", required=True, help="dRep dereplicated_genomes/")
    parser.add_argument("--project", required=True)
    parser.add_argument("--min-completeness", type=float, default=50.0)
    parser.add_argument("--max-contamination", type=float, default=10.0)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--table", required=True)
    args = parser.parse_args(argv)

    candidates = {row["candidate"]: row for row in read_tsv(args.candidates)}
    winners = [row["genome"] for row in read_tsv(args.winners, ",")]
    cluster_rows = read_tsv(args.clusters, ",")
    cluster_of = {row["genome"]: row["secondary_cluster"] for row in cluster_rows}
    cluster_size = Counter(cluster_of.values())

    kept = []
    for genome in winners:
        cand = candidates.get(genome)
        if cand is None:
            sys.exit(f"dRep winner {genome} is not in {args.candidates}")
        comp, cont = float(cand["completeness"]), float(cand["contamination"])
        if comp >= args.min_completeness and cont <= args.max_contamination:
            kept.append(cand)
    kept.sort(key=lambda c: (natural_key(c["group"]), natural_key(c["source_bin"])))

    os.makedirs(args.out_dir, exist_ok=True)
    fields = ["mag_id", "candidate", "group", "source_bin", "origin", "completeness",
              "contamination", "quality_checkm2", "size_bp", "n_contigs", "n50",
              "drep_cluster", "drep_cluster_size"]
    with open(args.table, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for i, cand in enumerate(kept, start=1):
            mag_id = f"{args.project}_MAG_{i}"
            shutil.copyfile(os.path.join(args.genomes_dir, cand["candidate"]),
                            os.path.join(args.out_dir, f"{mag_id}.fna"))
            cluster = cluster_of.get(cand["candidate"], "")
            writer.writerow({
                "mag_id": mag_id,
                "candidate": cand["candidate"],
                "group": cand["group"],
                "source_bin": cand["source_bin"],
                "origin": cand.get("origin", ""),
                "completeness": cand["completeness"],
                "contamination": cand["contamination"],
                "quality_checkm2": quality_class(float(cand["completeness"]),
                                                 float(cand["contamination"])),
                "size_bp": cand["size_bp"],
                "n_contigs": cand["n_contigs"],
                "n50": cand["n50"],
                "drep_cluster": cluster,
                "drep_cluster_size": cluster_size.get(cluster, 1),
            })
    print(f"{len(kept)} MAGs written ({len(winners)} dRep representatives before filtering)")


if __name__ == "__main__":
    main()
