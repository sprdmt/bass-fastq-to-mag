#!/usr/bin/env python3
"""Summarise Snakemake benchmark files next to the input sizes that drive the resource model.

For every benchmarks/<step...>/<name>.tsv the assembly group or sample is recognised from the
path (names from the sample sheet). Read volume (GB of FASTQ), assembly size (GB of contigs) and
the number of samples are attached, together with peak memory per GB, so that the coefficients
in workflow/resources.yaml can be checked against real runs.
"""

import argparse
import csv
import os


def file_gb(path):
    try:
        return os.path.getsize(path) / 1e9
    except OSError:
        return None


def read_samples(path):
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    for row in rows:
        row["assembly_group"] = row.get("assembly_group") or ""
    return rows


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--benchmarks", default="benchmarks")
    parser.add_argument("--samples", default="config/samples.tsv")
    parser.add_argument("--results", default="results")
    parser.add_argument("--project", default="", help="default group name for samples without one")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    samples = read_samples(args.samples) if os.path.exists(args.samples) else []
    for row in samples:
        if not row["assembly_group"]:
            row["assembly_group"] = args.project
    by_sample = {r["sample"]: r for r in samples}
    groups = sorted({r["assembly_group"] for r in samples})

    def reads_gb(sample):
        row = by_sample[sample]
        trimmed = [os.path.join(args.results, "qc", "trimmed", f"{sample}_R{m}.fastq.gz") for m in (1, 2)]
        paths = trimmed if all(os.path.exists(p) for p in trimmed) else [row["fq1"], row["fq2"]]
        sizes = [file_gb(p) for p in paths]
        return None if None in sizes else sum(sizes)

    out_rows = []
    for root, _dirs, files in os.walk(args.benchmarks):
        for fname in sorted(files):
            if not fname.endswith(".tsv"):
                continue
            path = os.path.join(root, fname)
            rel = os.path.relpath(path, args.benchmarks)[: -len(".tsv")]
            tokens = rel.replace(os.sep, "/").replace(".", "/").replace("_", "/").split("/")
            group = next((g for g in groups if g in rel.split(os.sep)), "")
            sample = next((s for s in by_sample if s in tokens or os.path.basename(rel).startswith(s)), "")
            if sample and not group:
                group = by_sample[sample]["assembly_group"]
            with open(path, newline="") as fh:
                bench = list(csv.DictReader(fh, delimiter="\t"))
            if not bench:
                continue
            b = bench[-1]
            group_samples = [s for s, r in by_sample.items() if r["assembly_group"] == group]
            reads = reads_gb(sample) if sample else (
                sum(filter(None, (reads_gb(s) for s in group_samples))) if group_samples else None)
            contigs = file_gb(os.path.join(args.results, "assembly", group, "contigs.fa")) if group else None
            rss = float(b.get("max_rss") or 0)
            out_rows.append({
                "step": rel,
                "group": group,
                "sample": sample,
                "n_samples": len(group_samples) if group else "",
                "reads_gb": round(reads, 3) if reads else "",
                "contigs_gb": round(contigs, 3) if contigs else "",
                "wall_min": round(float(b.get("s") or 0) / 60, 2),
                "max_rss_mb": round(rss, 1),
                "cpu_time_s": b.get("cpu_time", ""),
                "mean_load": b.get("mean_load", ""),
                "rss_mb_per_contigs_gb": round(rss / contigs, 1) if contigs else "",
                "rss_mb_per_reads_gb": round(rss / reads, 1) if reads else "",
            })

    fields = ["step", "group", "sample", "n_samples", "reads_gb", "contigs_gb", "wall_min",
              "max_rss_mb", "cpu_time_s", "mean_load", "rss_mb_per_contigs_gb", "rss_mb_per_reads_gb"]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(out_rows, key=lambda r: r["step"]))
    print(f"{len(out_rows)} benchmark files summarised in {args.out}")


if __name__ == "__main__":
    main()
