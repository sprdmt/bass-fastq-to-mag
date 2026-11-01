#!/usr/bin/env python3
"""Convert Prodigal/pprodigal GFF + protein FASTA into an anvi'o external-gene-calls table.

Ported from the rule `build_external_gene_calls` of the original workflow. Coordinates
become 0-based starts, strand becomes f/r, partial genes are flagged, and the amino-acid
sequence is attached by matching '<contig>_<gene number>' protein names.
"""

import argparse
import re


def simplify(text):
    return re.split(r"[ \t|#]", text)[0]


def read_proteins(faa):
    seqs, name, chunks = {}, None, []
    with open(faa) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks).rstrip("*")
                name, chunks = simplify(line[1:].strip()), []
            else:
                chunks.append(line.strip())
    if name is not None:
        seqs[name] = "".join(chunks).rstrip("*")
    return seqs


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gff", required=True)
    parser.add_argument("--faa", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--source", default="prodigal")
    parser.add_argument("--version", default="v2.6.3")
    args = parser.parse_args(argv)

    proteins = read_proteins(args.faa)
    kept = 0
    with open(args.gff) as gff, open(args.out, "w") as out:
        out.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\t"
                  "source\tversion\taa_sequence\n")
        for line in gff:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "CDS":
                continue
            contig, _src, _feat, start, end, _score, strand, _phase, attributes = parts
            contig = simplify(contig)
            attrs = dict(kv.split("=", 1) for kv in attributes.split(";") if "=" in kv)
            gene_id = attrs.get("ID", "")
            gene_number = gene_id.split("_")[-1] if "_" in gene_id else gene_id
            partial = 0 if attrs.get("partial", "00") == "00" else 1
            kept += 1
            out.write(
                f"{kept}\t{contig}\t{int(start) - 1}\t{int(end)}\t{'f' if strand == '+' else 'r'}\t"
                f"{partial}\t1\t{args.source}\t{args.version}\t"
                f"{proteins.get(f'{contig}_{gene_number}', '')}\n"
            )
    print(f"Wrote {kept} external gene calls to {args.out}")


if __name__ == "__main__":
    main()
