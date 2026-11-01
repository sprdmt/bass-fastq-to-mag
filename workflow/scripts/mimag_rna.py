#!/usr/bin/env python3
"""rRNA / tRNA summary and MIMAG quality class per MAG.

Inputs
------
* barrnap GFF3 per MAG (``<mag_id>.gff``; ``<mag_id>.fna.gff`` etc. also
  accepted). The workflow runs barrnap with ``--kingdom bac`` or ``arc``
  chosen from the GTDB-Tk domain; both model sets name their features
  ``5S_rRNA`` / ``16S_rRNA`` / ``23S_rRNA``, so one parser serves both.
  Eukaryotic/organellar models (5_8S, 18S, 28S, 12S) are counted separately
  and never satisfy the prokaryotic MIMAG rRNA criterion.
* tRNAscan-SE tabular output per MAG (``-o`` file of ``tRNAscan-SE -B``/``-A``).
* completeness / contamination: ``--mags-table`` (workflow mags.tsv:
  mag_id, completeness, contamination; CheckM2 values) or ``--checkm2``
  (CheckM2 quality_report.tsv: Name, Completeness, Contamination).

Definitions
-----------
* rrna_5S / rrna_16S / rrna_23S: 1 if barrnap reports at least one hit for the
  gene (partial hits count unless ``--exclude-partial-rrna``), else 0.
* trna_unique: number of distinct standard amino-acid tRNA isotypes (of 20)
  among non-pseudo tRNAscan-SE predictions. iMet/fMet count as Met and Ile2
  (tRNA-Ile CAU) as Ile; SeC, Pyl, Sup (suppressor) and Undet are not counted.
  ``--trna-mode legacy`` reproduces the earlier ENA-prep Snakefile instead
  (trailing digits stripped, iMet->Met, every other type name counted,
  pseudogenes not excluded).
* rRNA_presence (genome_uploader definition): True if 5S, 16S and 23S rRNA
  are all detected AND at least 18 tRNAs (trna_unique >= 18).
* mimag_quality:
  high-quality   = MIMAG high-quality draft (Bowers et al. 2017): completeness > 90,
                   contamination < 5, 23S+16S+5S rRNA and >= 18 tRNAs
  near-complete  = completeness >= 90 and contamination <= 5 but not MIMAG
                   high-quality (rRNA/tRNA criteria not met; Almeida et al. 2019)
  medium-quality = completeness >= 50 and contamination < 10 (and not above)
  low-quality    = everything else (contamination >= 10 is flagged in ``note``)
* genome_uploader_quality: what genome_uploader 2.5.1 will write as the
  "assembly quality" attribute: HQ if completeness >= 90, contamination <= 5
  (after rounding to 2 decimals) and rRNA_presence, else MQ.
* ena_quality_conflict: True when genome_uploader would write the HQ text for a
  genome that is not MIMAG high-quality (only possible at exactly 90.00 %
  completeness or 5.00 % contamination); ena_genome_table.py stops on these.

Usage
-----
  mimag_rna.py --barrnap-dir results/mags/rna/barrnap \
               --trnascan-dir results/mags/rna/trnascan \
               --mags-table results/mags/mags.tsv \
               --out results/mags/mimag_rna.tsv [--gtdbtk results/mags/gtdbtk/gtdbtk.summary.tsv]
"""
from __future__ import annotations

import argparse
import logging
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import (EnaError, fmt_bool, fmt_float, read_tsv, setup_logging,  # noqa: E402
                        strip_suffixes, write_tsv)

LOG = logging.getLogger("bfm.mimag_rna")

STANDARD_AA = ("Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile",
               "Leu", "Lys", "Met", "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val")
_STD_LOWER = {a.lower(): a for a in STANDARD_AA}
_ALIASES = {"imet": "Met", "fmet": "Met", "ile2": "Ile"}
PROK_RRNA = ("5S", "16S", "23S")
GFF_SUFFIXES = (".gff3", ".gff")
TRNA_SUFFIXES = (".txt", ".tsv", ".out", ".trnascan")
SEQ_SUFFIXES = (".fna", ".fa", ".fasta")
MIN_TRNA = 18

OUT_COLUMNS = ["mag_id", "domain", "rrna_5S", "rrna_16S", "rrna_23S", "n_5S", "n_16S", "n_23S",
               "n_rrna_partial", "n_rrna_other", "trna_unique", "trna_isotypes", "n_trna",
               "n_trna_pseudo", "rRNA_presence", "completeness", "contamination",
               "mimag_quality", "genome_uploader_quality", "ena_quality_conflict", "note"]


# --------------------------------------------------------------------------
# barrnap
# --------------------------------------------------------------------------
def _gff_attributes(col9: str) -> Dict[str, str]:
    attrs = {}
    for part in col9.strip().split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            attrs[k.strip()] = v.strip()
    return attrs


def classify_rrna(attrs: Dict[str, str]) -> Optional[str]:
    """Return '5S', '16S', '23S', '5.8S', '18S', '28S', '12S' ... from barrnap attributes."""
    name = attrs.get("Name", "")
    m = re.match(r"^(\d+(?:[._]\d+)?)S_rRNA", name)
    if not m:
        m = re.search(r"(\d+(?:\.\d+)?)S ribosomal RNA", attrs.get("product", ""), re.I)
    if not m:
        return None
    return m.group(1).replace("_", ".") + "S"


def parse_barrnap_gff(path, exclude_partial: bool = False) -> Dict[str, int]:
    """Count rRNA genes in a barrnap GFF3 file.

    Returns counts for '5S', '16S', '23S', 'partial' (hits flagged partial,
    counted or not depending on ``exclude_partial``) and 'other' (non-
    prokaryotic models such as 18S/28S/5.8S/12S).
    """
    counts = {"5S": 0, "16S": 0, "23S": 0, "partial": 0, "other": 0}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2] != "rRNA":
                continue
            attrs = _gff_attributes(cols[8])
            gene = classify_rrna(attrs)
            if gene is None:
                continue
            partial = "partial" in attrs.get("product", "").lower() or "partial" in attrs.get("note", "").lower() \
                or "aligned only" in attrs.get("note", "").lower()
            if partial:
                counts["partial"] += 1
                if exclude_partial:
                    continue
            if gene in PROK_RRNA:
                counts[gene] += 1
            else:
                counts["other"] += 1
    return counts


# --------------------------------------------------------------------------
# tRNAscan-SE
# --------------------------------------------------------------------------
def normalise_isotype(raw: str, mode: str = "standard20") -> Optional[str]:
    s = (raw or "").strip()
    if not s:
        return None
    if mode == "legacy":
        if s.lower() == "imet":
            return "Met"
        s = re.sub(r"(\d+)$", "", s)
        if not s or s.lower() in ("pseudo", "undet"):
            return None
        return s[0].upper() + s[1:].lower()
    low = s.lower()
    if low in _ALIASES:
        return _ALIASES[low]
    return _STD_LOWER.get(low)


def parse_trnascan(path, mode: str = "standard20") -> Dict[str, object]:
    """Parse tRNAscan-SE tabular output (``-o``)."""
    isotypes = set()
    n_trna = 0
    n_pseudo = 0
    with open(path) as fh:
        for line in fh:
            if not line.strip():
                continue
            cols = [c.strip() for c in line.rstrip("\n").split("\t")]
            if len(cols) < 5:
                cols = line.split()
                if len(cols) < 5:
                    continue
            first = cols[0]
            if first.startswith(("Sequence", "Name", "--------")):
                continue
            try:
                int(cols[1])
            except ValueError:
                continue  # not a data row
            n_trna += 1
            note = " ".join(cols[9:]).lower() if len(cols) > 9 else ""
            is_pseudo = "pseudo" in note
            if is_pseudo:
                n_pseudo += 1
                if mode != "legacy":
                    continue
            iso = normalise_isotype(cols[4], mode)
            if iso:
                isotypes.add(iso)
    return {"isotypes": sorted(isotypes), "n_unique": len(isotypes), "n_trna": n_trna, "n_pseudo": n_pseudo}


# --------------------------------------------------------------------------
# quality classes
# --------------------------------------------------------------------------
def rrna_presence(has5: bool, has16: bool, has23: bool, n_trna_unique: int, min_trna: int = MIN_TRNA) -> bool:
    """genome_uploader's rRNA_presence: all of 5S/16S/23S and >= 18 tRNAs."""
    return bool(has5 and has16 and has23 and n_trna_unique >= min_trna)


def mimag_quality(completeness: float, contamination: float, has5: bool, has16: bool, has23: bool,
                  n_trna_unique: int, min_trna: int = MIN_TRNA) -> Tuple[str, str]:
    """Return (class, note): high-quality / near-complete / medium-quality / low-quality."""
    rna_ok = has5 and has16 and has23 and n_trna_unique >= min_trna
    if completeness > 90 and contamination < 5 and rna_ok:
        return "high-quality", ""
    if completeness >= 90 and contamination <= 5:
        parts = []
        missing = [g for g, ok in (("5S", has5), ("16S", has16), ("23S", has23)) if not ok]
        if missing:
            parts.append("missing " + "/".join(missing) + " rRNA")
        if n_trna_unique < min_trna:
            parts.append(f"{n_trna_unique} tRNA isotypes")
        if rna_ok:
            parts.append("completeness/contamination at the MIMAG boundary (needs > 90 and < 5)")
        return "near-complete", "not MIMAG high-quality: " + ", ".join(parts)
    if completeness >= 50 and contamination < 10:
        return "medium-quality", ""
    return "low-quality", ("contamination >= 10% fails MIMAG criteria" if contamination >= 10 else "")


def genome_uploader_quality(completeness: float, contamination: float, presence: bool) -> str:
    """Replicates genome_uploader 2.5.1 compute_mag_quality (HQ/MQ only)."""
    c, k = round(float(completeness), 2), round(float(contamination), 2)
    return "HQ" if (c >= 90 and k <= 5 and presence) else "MQ"


# --------------------------------------------------------------------------
# inputs
# --------------------------------------------------------------------------
def index_dir(directory, suffixes: Sequence[str]) -> Dict[str, Path]:
    """Map MAG id -> file, stripping the given suffix(es) and a FASTA suffix."""
    idx: Dict[str, Path] = {}
    for p in sorted(Path(directory).iterdir()):
        if not p.is_file() or not p.name.endswith(tuple(suffixes)):
            continue
        mag = strip_suffixes(strip_suffixes(p.name, suffixes), SEQ_SUFFIXES)
        if mag in idx:
            raise EnaError(f"two files for MAG {mag}: {idx[mag].name}, {p.name}")
        idx[mag] = p
    return idx


def read_quality(mags_table=None, checkm2=None) -> Dict[str, Tuple[float, float]]:
    out: Dict[str, Tuple[float, float]] = {}
    if mags_table:
        for r in read_tsv(mags_table, required=("mag_id", "completeness", "contamination")):
            out[r["mag_id"]] = (float(r["completeness"]), float(r["contamination"]))
    elif checkm2:
        for r in read_tsv(checkm2, required=("Name", "Completeness", "Contamination")):
            out[r["Name"]] = (float(r["Completeness"]), float(r["Contamination"]))
    else:
        raise EnaError("give --mags-table or --checkm2")
    return out


def read_gtdbtk_domains(paths: Sequence[str]) -> Dict[str, str]:
    dom: Dict[str, str] = {}
    for p in paths or []:
        for r in read_tsv(p, required=("user_genome", "classification")):
            c = r["classification"]
            if "Archaea" in c:
                dom[r["user_genome"]] = "Archaea"
            elif "Bacteria" in c:
                dom[r["user_genome"]] = "Bacteria"
            else:
                dom[r["user_genome"]] = "unknown"
    return dom


def summarise(quality: Dict[str, Tuple[float, float]], barrnap: Dict[str, Path], trnascan: Dict[str, Path],
              domains: Optional[Dict[str, str]] = None, exclude_partial: bool = False,
              trna_mode: str = "standard20", min_trna: int = MIN_TRNA,
              allow_missing: bool = False) -> List[Dict[str, object]]:
    domains = domains or {}
    rows = []
    missing = []
    for mag in sorted(quality, key=_natural_key):
        comp, cont = quality[mag]
        notes = []
        if mag in barrnap:
            rr = parse_barrnap_gff(barrnap[mag], exclude_partial=exclude_partial)
        else:
            missing.append(f"{mag}: barrnap GFF")
            rr = {"5S": 0, "16S": 0, "23S": 0, "partial": 0, "other": 0}
            notes.append("no barrnap file")
        if mag in trnascan:
            tr = parse_trnascan(trnascan[mag], mode=trna_mode)
        else:
            missing.append(f"{mag}: tRNAscan-SE output")
            tr = {"isotypes": [], "n_unique": 0, "n_trna": 0, "n_pseudo": 0}
            notes.append("no tRNAscan-SE file")
        h5, h16, h23 = rr["5S"] > 0, rr["16S"] > 0, rr["23S"] > 0
        presence = rrna_presence(h5, h16, h23, tr["n_unique"], min_trna)
        qclass, qnote = mimag_quality(comp, cont, h5, h16, h23, tr["n_unique"], min_trna)
        gu = genome_uploader_quality(comp, cont, presence)
        if qnote:
            notes.append(qnote)
        rows.append({
            "mag_id": mag, "domain": domains.get(mag, ""),
            "rrna_5S": int(h5), "rrna_16S": int(h16), "rrna_23S": int(h23),
            "n_5S": rr["5S"], "n_16S": rr["16S"], "n_23S": rr["23S"],
            "n_rrna_partial": rr["partial"], "n_rrna_other": rr["other"],
            "trna_unique": tr["n_unique"], "trna_isotypes": ",".join(tr["isotypes"]),
            "n_trna": tr["n_trna"], "n_trna_pseudo": tr["n_pseudo"],
            "rRNA_presence": fmt_bool(presence),
            "completeness": fmt_float(comp), "contamination": fmt_float(cont),
            "mimag_quality": qclass,
            "genome_uploader_quality": gu,
            "ena_quality_conflict": fmt_bool(gu == "HQ" and qclass != "high-quality"),
            "note": "; ".join(notes),
        })
    if missing and not allow_missing:
        raise EnaError("missing RNA prediction files:\n  " + "\n  ".join(missing))
    return rows


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--barrnap-dir", required=True, help="directory with one barrnap GFF per MAG")
    ap.add_argument("--trnascan-dir", required=True, help="directory with one tRNAscan-SE table per MAG")
    q = ap.add_mutually_exclusive_group(required=True)
    q.add_argument("--mags-table", help="workflow mags.tsv (mag_id, completeness, contamination)")
    q.add_argument("--checkm2", help="CheckM2 quality_report.tsv (Name, Completeness, Contamination)")
    ap.add_argument("--gtdbtk", nargs="*", default=[], help="GTDB-Tk summary file(s), only to report the domain")
    ap.add_argument("--out", required=True, help="output TSV")
    ap.add_argument("--exclude-partial-rrna", action="store_true",
                    help="do not count barrnap hits flagged as partial")
    ap.add_argument("--trna-mode", choices=("standard20", "legacy"), default="standard20")
    ap.add_argument("--min-trna", type=int, default=MIN_TRNA)
    ap.add_argument("--allow-missing", action="store_true",
                    help="treat a missing barrnap/tRNAscan file as zero hits instead of failing")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        quality = read_quality(args.mags_table, args.checkm2)
        rows = summarise(quality, index_dir(args.barrnap_dir, GFF_SUFFIXES),
                         index_dir(args.trnascan_dir, TRNA_SUFFIXES),
                         read_gtdbtk_domains(args.gtdbtk), args.exclude_partial_rrna,
                         args.trna_mode, args.min_trna, args.allow_missing)
    except EnaError as exc:
        LOG.error("%s", exc)
        return 1
    write_tsv(args.out, rows, OUT_COLUMNS)
    n = {k: sum(1 for r in rows if r["mimag_quality"] == k)
         for k in ("high-quality", "near-complete", "medium-quality", "low-quality")}
    LOG.info("%d MAGs: rRNA_presence True for %d; %s", len(rows),
             sum(1 for r in rows if r["rRNA_presence"] == "True"), n)
    return 0


if __name__ == "__main__":
    sys.exit(main())
