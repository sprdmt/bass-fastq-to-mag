#!/usr/bin/env python3
"""Write a contig-to-bin table (contig<TAB>bin, no header) from a directory of bin FASTA files.

Bin files are recognised by content (first non-empty line starts with '>'), so files with
or without .fa/.fna/.fasta extensions and gzip-compressed files are all accepted. The bin
name is the file name without FASTA/gzip extensions, optionally prefixed. Contig names are
the first whitespace-delimited token of each header.
"""

import argparse
import gzip
import os
import sys

FASTA_EXTENSIONS = (".fasta", ".fas", ".fna", ".fa")


def bin_name(filename, prefix=""):
    name = filename[:-3] if filename.endswith(".gz") else filename
    for ext in FASTA_EXTENSIONS:
        if name.endswith(ext):
            name = name[: -len(ext)]
            break
    return prefix + name


def open_text(path):
    with open(path, "rb") as fh:
        magic = fh.read(2)
    return gzip.open(path, "rt") if magic == b"\x1f\x8b" else open(path, "rt")


def is_fasta(path):
    try:
        with open_text(path) as fh:
            for line in fh:
                if line.strip():
                    return line.startswith(">")
    except (OSError, UnicodeDecodeError, EOFError):
        return False
    return False


def table_rows(bins_dir, prefix="", exclude=()):
    """Yield (contig, bin) pairs in a deterministic order."""
    for filename in sorted(os.listdir(bins_dir)):
        path = os.path.join(bins_dir, filename)
        if filename.startswith(".") or filename in exclude or not os.path.isfile(path):
            continue
        if not is_fasta(path):
            continue
        name = bin_name(filename, prefix)
        with open_text(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    yield line[1:].split()[0], name


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bins-dir", required=True, help="directory with one FASTA file per bin")
    parser.add_argument("--out", required=True, help="output table (contig<TAB>bin)")
    parser.add_argument("--prefix", default="", help="prefix added to every bin name")
    parser.add_argument("--exclude", nargs="*", default=[], help="file names to ignore")
    args = parser.parse_args(argv)

    if not os.path.isdir(args.bins_dir):
        sys.exit(f"Bin directory not found: {args.bins_dir}")
    seen = {}
    bins = set()
    with open(args.out, "w") as out:
        for contig, name in table_rows(args.bins_dir, args.prefix, set(args.exclude)):
            out.write(f"{contig}\t{name}\n")
            seen[contig] = seen.get(contig, 0) + 1
            bins.add(name)
    multi = sum(1 for n in seen.values() if n > 1)
    print(f"{len(bins)} bins, {len(seen)} contigs, {multi} contigs in more than one bin")
    if not bins:
        print(f"WARNING: no bin FASTA files found in {args.bins_dir}", file=sys.stderr)


if __name__ == "__main__":
    main()
