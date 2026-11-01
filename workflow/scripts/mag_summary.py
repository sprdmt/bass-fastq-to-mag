#!/usr/bin/env python3
"""Join the per-MAG tables into results/mags/mag_summary.tsv.

Inputs: mags.tsv (names, CheckM2 quality, dRep cluster), GTDB-Tk summary, MIMAG rRNA/tRNA
table and coverage table. Missing optional inputs (e.g. GTDB-Tk disabled) leave empty columns.
"""

import argparse

import pandas as pd


def read_optional(path, **kwargs):
    if not path:
        return None
    try:
        return pd.read_csv(path, sep="\t", **kwargs)
    except (FileNotFoundError, pd.errors.EmptyDataError):
        return None


def join(table, other):
    """Left-join `other` on mag_id; columns that `table` already has are kept from `table`."""
    other = other.drop(columns=[c for c in other.columns if c != "mag_id" and c in table.columns])
    return table.merge(other, on="mag_id", how="left")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mags", required=True)
    parser.add_argument("--gtdbtk", default="")
    parser.add_argument("--mimag", required=True)
    parser.add_argument("--coverage", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    mags = pd.read_csv(args.mags, sep="\t")
    table = mags.copy()

    gtdb = read_optional(args.gtdbtk)
    if gtdb is not None and "user_genome" in gtdb.columns:
        keep = [c for c in ("user_genome", "classification", "classification_method", "note") if c in gtdb.columns]
        gtdb = gtdb[keep].rename(columns={"user_genome": "mag_id", "classification": "gtdb_classification",
                                          "classification_method": "gtdb_method", "note": "gtdb_note"})
        table = join(table, gtdb)
    else:
        table["gtdb_classification"] = ""

    mimag = pd.read_csv(args.mimag, sep="\t")
    mimag = mimag.drop(columns=[c for c in ("completeness", "contamination", "group") if c in mimag.columns])
    table = join(table, mimag)

    coverage = pd.read_csv(args.coverage, sep="\t")
    coverage = coverage.drop(columns=[c for c in ("group",) if c in coverage.columns])
    table = join(table, coverage)

    first = ["mag_id", "group", "completeness", "contamination", "mimag_quality",
             "gtdb_classification", "size_bp", "n_contigs", "n50", "genome_coverage"]
    ordered = [c for c in first if c in table.columns] + [c for c in table.columns if c not in first]
    table[ordered].to_csv(args.out, sep="\t", index=False)
    counts = table["mimag_quality"].value_counts().to_dict() if "mimag_quality" in table else {}
    print(f"{len(table)} MAGs; MIMAG classes: {counts}")


if __name__ == "__main__":
    main()
