#!/usr/bin/env python3
"""Prepare a taxconverter table for `vamb bin taxvamb`.

1. Lineages are cut at the first 'unknown' or 'uc_*' rank.
2. A taxon name that occurs under more than one parent is renamed '<parent>_<name>' after
   its first occurrence, because TaxVAMB requires each name to have a single parent.
3. Every contig of the assembly is written, in FASTA order; contigs without a prediction
   get an empty lineage.

Since March 2025 the NCBI taxonomy names the top rank of Bacteria, Archaea and Eukaryota
"domain" (Viruses: "acellular root") instead of "superkingdom". MMseqs2 18 does not accept
"domain" in --lca-ranks, so with a database built from a current taxdump the first rank of
every --lca-ranks lineage is "unknown" and step 1 would empty all lineages. With --mmseqs, the
domain is read from the full lineage (--tax-lineage 1, last column): 'd_<name>' (superkingdom
in older NCBI and GTDB databases) or '-_<name>' for Bacteria, Archaea, Eukaryota and Viruses.

Ported from the inline step `fix_taxonomy_for_vamb` of the original workflow.
"""

import argparse
import sys


DOMAINS = {"Bacteria", "Archaea", "Eukaryota", "Viruses"}


def domain_from_full_lineage(full):
    """Domain from an MMseqs2 full lineage such as '-_cellular organisms;-_Bacteria;p_...'."""
    for tok in (full or "").split(";"):
        code, sep, name = tok.strip().partition("_")
        if not sep:
            continue
        if code == "d" or (code == "-" and name in DOMAINS):
            return name
    return None


def read_mmseqs_domains(path):
    """contig -> domain from the last column of an MMseqs2 taxonomy TSV (--tax-lineage 1)."""
    domains = {}
    with open(path) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            dom = domain_from_full_lineage(parts[-1])
            if dom:
                domains[parts[0].strip()] = dom
    return domains


def fill_domain(raw, domain):
    """Replace an 'unknown' first rank by the domain; returns (lineage, replaced)."""
    tokens = [t.strip() for t in (raw or "").split(";")]
    if domain and tokens and tokens[0] in ("", "unknown"):
        tokens[0] = domain
        return ";".join(tokens), True
    return raw, False


def read_contig_ids(fasta):
    ids = []
    with open(fasta) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


class LineageCleaner:
    def __init__(self):
        self.first_parent = {}
        self.n_ambiguous = 0

    def clean(self, raw):
        tokens = [t.strip() for t in (raw or "").split(";") if t.strip()]
        kept = []
        for tok in tokens:
            if tok == "unknown" or tok.startswith("uc_"):
                break
            kept.append(tok)
        out = []
        for i, tok in enumerate(kept):
            parent = kept[i - 1] if i > 0 else "root"
            previous = self.first_parent.get(tok)
            if previous is None:
                self.first_parent[tok] = parent
                out.append(tok)
            elif previous == parent:
                out.append(tok)
            else:
                out.append(f"{parent}_{tok}")
                self.n_ambiguous += 1
        return ";".join(out)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--contigs", required=True)
    parser.add_argument("--taxonomy", required=True, help="taxconverter output (contigs, predictions)")
    parser.add_argument("--mmseqs", help="MMseqs2 taxonomy TSV with the full lineage (--tax-lineage 1)")
    parser.add_argument("--out", required=True)
    args = parser.parse_args(argv)

    contigs = read_contig_ids(args.contigs)
    domains = read_mmseqs_domains(args.mmseqs) if args.mmseqs else {}
    cleaner = LineageCleaner()
    lineages = {}
    n_rows = n_filled = 0
    with open(args.taxonomy) as fh:
        fh.readline()  # header: contigs<TAB>predictions
        for line in fh:
            if not line.strip():
                continue
            n_rows += 1
            parts = line.rstrip("\n").split("\t")
            name = parts[0].strip()
            raw, filled = fill_domain(parts[1] if len(parts) > 1 else "", domains.get(name))
            n_filled += filled
            lineages[name] = cleaner.clean(raw)

    missing = with_lineage = 0
    with open(args.out, "w") as out:
        out.write("contigs\tpredictions\n")
        for name in contigs:
            if name not in lineages:
                missing += 1
            lineage = lineages.get(name, "")
            with_lineage += bool(lineage)
            out.write(f"{name}\t{lineage}\n")

    print(f"FASTA contigs: {len(contigs)}")
    print(f"Taxonomy rows read: {n_rows}")
    print(f"Contigs missing from the taxonomy table: {missing}")
    print(f"Contigs with a lineage: {with_lineage}")
    print(f"Domains taken from the MMseqs2 full lineage: {n_filled}")
    print(f"Ambiguous names disambiguated: {cleaner.n_ambiguous}")
    if with_lineage == 0:
        sys.exit(
            "ERROR: no contig has a lineage, so TaxVAMB cannot build its taxonomy tree. Check the "
            "MMseqs2 output: with a database built from an NCBI taxdump of March 2025 or later, "
            "pass --mmseqs with an MMseqs2 table that has the full lineage (--tax-lineage 1)."
        )


if __name__ == "__main__":
    main()
