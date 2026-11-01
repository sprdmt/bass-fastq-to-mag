#!/usr/bin/env python3
"""Triage the Binette bins of one assembly group for manual refinement in anvi'o.

Decisions (CheckM2 values from the Binette report):
  refine  completeness >= --refine-min-completeness and contamination > --refine-contamination-above
          (worth cleaning by hand; refinement removes contigs, so completeness cannot rise)
  keep    otherwise, if the bin passes the final MAG thresholds: enters the MAG set as it is
  drop    everything else
With --review all, every bin is put up for review (decision 'review') and nothing is kept
automatically.

Outputs: triage table, anvi'o collection tables (all bins -> 'Binette', bins to review ->
'Refine') and a markdown guide with one anvi-refine command per bin.
"""

import argparse
import csv
import os

REVIEW = ("refine", "review")


def read_binette(report):
    bins = {}
    with open(report, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            row = {k.strip(): (v or "").strip() for k, v in row.items() if k}
            name = f"bin_{row['bin_id']}" if row.get("bin_id") else row["name"]  # Binette 1.1 / 1.2
            bins[name] = {
                "completeness": float(row["completeness"]),
                "contamination": float(row["contamination"]),
                "size": int(float(row.get("size") or 0)),
                "origin": row.get("origin", ""),
            }
    return bins


def bin_files(bins_dir):
    files = {}
    for fname in os.listdir(bins_dir):
        stem, ext = os.path.splitext(fname)
        if ext in (".fa", ".fna", ".fasta"):
            files[stem] = os.path.join(bins_dir, fname)
    return files


def contigs(path):
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split()[0]


def decide(q, a):
    if a.review == "all":
        return "review", "anvio.review: all"
    if q["completeness"] >= a.refine_min_completeness and q["contamination"] > a.refine_contamination_above:
        return "refine", (f"completeness >= {a.refine_min_completeness:g}%, "
                          f"contamination > {a.refine_contamination_above:g}%")
    if q["completeness"] >= a.keep_min_completeness and q["contamination"] <= a.keep_max_contamination:
        return "keep", "passes the MAG thresholds"
    return "drop", "below the MAG thresholds"


def write_howto(a, rows, path):
    todo = [r for r in rows if r["decision"] in REVIEW]
    counts = {d: sum(r["decision"] == d for r in rows) for d in ("keep", "refine", "review", "drop")}
    profile, db = os.path.abspath(a.profile), os.path.abspath(a.contigs_db)
    title = f"# Manual refinement - group {a.group}" if a.group else "# Manual refinement"
    lines = [
        title,
        "",
        f"{len(rows)} Binette bins: {counts['keep']} kept as they are, "
        f"{counts['refine'] + counts['review']} to review, {counts['drop']} dropped "
        f"(details: triage.tsv).",
        "",
    ]
    if a.review == "all":
        lines.append("All bins are in the collection **Refine** (`anvio.review: all`).")
    else:
        lines.append(
            f"The collection **Refine** holds the bins with >= {a.refine_min_completeness:g}% completeness "
            f"and > {a.refine_contamination_above:g}% contamination, including bins that already pass the "
            f"MAG thresholds (>= {a.keep_min_completeness:g}% / <= {a.keep_max_contamination:g}%) but could "
            "become near-complete or high-quality. All other bins that pass the MAG thresholds enter the "
            "final set unchanged; the collection **Binette** keeps all automatic bins for reference.")
    lines += ["", "| bin | completeness | contamination | size (bp) | command |", "|---|---|---|---|---|"]
    for r in todo:
        cmd = f"anvi-refine -p {profile} -c {db} -C Refine -b {r['bin']}"
        lines.append(f"| {r['bin']} | {r['completeness']:.1f} | {r['contamination']:.1f} | {r['size']} | `{cmd}` |")
    if not todo:
        lines.append("| - | | | | nothing to refine |")
    lines += [
        "",
        "## Notes",
        "",
        "* Work only in **Refine**: split a bin, remove contaminating contigs, or delete bins that cannot",
        "  be rescued. Bins left untouched are used as they are.",
        "* Check progress with",
        f"  `anvi-estimate-genome-completeness -p {profile} -c {db} -C Refine`.",
        "* On a cluster, add `--server-only -P 8080` to anvi-refine and open a tunnel from your computer:",
        "  `ssh -L 8080:localhost:8080 <user>@<node>`, then browse to http://localhost:8080.",
    ]
    if a.image:
        lines.append(f"* Without a local anvi'o: `apptainer exec --bind $PWD {a.image} anvi-refine ...`.")
    lines += [
        "* Your edits live in PROFILE.db. Back them up with",
        f"  `anvi-export-collection -p {profile} -C Refine -O Refine_backup`, and run later steps with",
        "  `--rerun-triggers mtime` so that a workflow update does not rebuild the anvi'o workspace.",
        "",
        "## Finish",
        "",
        "Set `mags: {source: anvio}` in config.yaml and run `snakemake mags ...`. The final candidates",
        "are the kept bins plus the bins of **Refine** (re-checked with CheckM2), dereplicated with dRep",
        "across all groups.",
        "",
    ]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    return counts


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bins-dir", required=True)
    p.add_argument("--report", required=True, help="Binette final_bins_quality_reports.tsv")
    p.add_argument("--keep-min-completeness", type=float, required=True)
    p.add_argument("--keep-max-contamination", type=float, required=True)
    p.add_argument("--refine-min-completeness", type=float, default=90.0)
    p.add_argument("--refine-contamination-above", type=float, default=5.0)
    p.add_argument("--review", choices=["contaminated", "all"], default="contaminated")
    p.add_argument("--group", default="")
    p.add_argument("--profile", default="PROFILE.db")
    p.add_argument("--contigs-db", default="contigs.db")
    p.add_argument("--image", default="", help="anvi'o container image, for the guide")
    p.add_argument("--triage", required=True)
    p.add_argument("--all-table", required=True)
    p.add_argument("--refine-table", required=True)
    p.add_argument("--howto", required=True)
    a = p.parse_args(argv)

    quality = read_binette(a.report)
    files = bin_files(a.bins_dir)
    missing = sorted(set(quality) - set(files))
    if missing:
        print(f"WARNING: {len(missing)} bins of the report have no FASTA file: {', '.join(missing[:5])}")
    order = sorted((n for n in quality if n in files),
                   key=lambda n: (-quality[n]["completeness"], quality[n]["contamination"]))
    rows = []
    for name in order:
        decision, reason = decide(quality[name], a)
        rows.append({"bin": name, **quality[name], "decision": decision, "reason": reason})

    with open(a.triage, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=["bin", "completeness", "contamination", "size",
                                                           "origin", "decision", "reason"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    with open(a.all_table, "w") as all_fh, open(a.refine_table, "w") as ref_fh:
        for row in rows:
            for contig in contigs(files[row["bin"]]):
                all_fh.write(f"{contig}\t{row['bin']}\n")
                if row["decision"] in REVIEW:
                    ref_fh.write(f"{contig}\t{row['bin']}\n")
    counts = write_howto(a, rows, a.howto)
    print(f"{len(rows)} bins: " + ", ".join(f"{k} {v}" for k, v in counts.items() if v))


if __name__ == "__main__":
    main()
