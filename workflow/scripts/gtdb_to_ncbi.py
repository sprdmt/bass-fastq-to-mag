#!/usr/bin/env python3
"""GTDB-Tk classification -> NCBI lineage for ENA (genome_uploader NCBI_lineage column).

Output format
-------------
``ncbi_lineage`` is written as ``d__Bacteria;p__...;c__...;o__...;f__...;g__...;s__``
with NCBI names. genome_uploader 2.5.1 (and 3.x) only accepts name lineages
whose elements carry a rank prefix containing ``__``; it walks the lineage
from the most specific rank upwards and asks the ENA taxonomy service for the
first submittable name (``<Genus> sp.``, ``<Family> bacterium``, ``uncultured
...``). By default the lineage is truncated after the genus (``--max-rank
genus``), following ENA's guidance to identify organism-level environmental
taxa "up to genus level"; ``--max-rank species`` keeps NCBI species.

Methods
-------
majority_vote (default)
    GTDB-Tk's ``gtdb_to_ncbi_majority_vote.py`` (shipped with GTDB-Tk). For
    each genome it takes the GTDB reference genomes around the genome's
    placement in the GTDB-Tk reference tree(s) and assigns, per rank, the NCBI
    taxon carried by >= 50 % of their NCBI classifications (from the GTDB
    metadata files of the same release as the GTDB-Tk reference data). The
    script needs the GTDB-Tk package, so the workflow runs it in the GTDB-Tk
    container and this script post-processes the table (``--mv-table``).
    ``--run-mv`` calls the script directly when GTDB-Tk is installed.
mapping_table
    The route used for the earlier ENA submissions: GTDB "gtdb_vs_ncbi"
    workbooks (sheets Species/Genus/Family/Order/Class/Phylum, columns
    ``GTDB Rxxx <rank>`` and ``List of NCBI <ranks>``). For the most specific
    GTDB rank whose taxon maps to a named NCBI taxon with >= ``--min-percent``
    (default 50) of the genomes, the NCBI lineage is read from an NCBI taxdump;
    otherwise the GTDB name itself is looked up in the taxdump (direct match),
    then the next rank up, finally the domain.

consensus
    Both of the above; keeps the ranks (from the domain down) on which they
    agree and empties everything from the first disagreement. ``--cross-check``
    instead keeps the majority-vote lineage and only reports disagreements with
    the mapping-table route (column ``cross_check``) for manual review.

Post-processing (all methods): missing or "Unclassified" domains are filled
from the GTDB domain, ranks below ``--max-rank`` are emptied, and with
``--taxdump`` every remaining name is checked against NCBI names.dmp.

Usage
-----
  gtdb_to_ncbi.py --gtdbtk results/mags/gtdbtk/gtdbtk.summary.tsv \
      --mv-table results/ena/gtdbtk_majority_vote.tsv --out results/ena/ncbi_lineage.tsv
  gtdb_to_ncbi.py --gtdbtk SUMMARY --run-mv --gtdbtk-dir results/mags/gtdbtk \
      --bac120-metadata bac120_metadata_r226.tsv.gz --ar53-metadata ar53_metadata_r226.tsv.gz --out OUT
  gtdb_to_ncbi.py --gtdbtk SUMMARY --method mapping_table \
      --mapping-table gtdb_vs_ncbi_r226_bacteria.xlsx gtdb_vs_ncbi_r226_archaea.xlsx \
      --taxdump taxdmp_2025-12-01.zip --out OUT
"""
from __future__ import annotations

import argparse
import io
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
import zipfile
from array import array
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Set, Tuple
import xml.etree.ElementTree as ET

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import EnaError, read_tsv, setup_logging, write_tsv  # noqa: E402

LOG = logging.getLogger("bfm.gtdb_to_ncbi")

RANK_LETTERS = ("d", "p", "c", "o", "f", "g", "s")
RANK_NAMES = ("domain", "phylum", "class", "order", "family", "genus", "species")
NCBI_RANK_TO_LETTER = {"superkingdom": "d", "domain": "d", "phylum": "p", "class": "c", "order": "o",
                       "family": "f", "genus": "g", "species": "s"}
SHEET_RANKS = {"species": "s", "genus": "g", "family": "f", "order": "o", "class": "c", "phylum": "p"}
DOMAIN_TAXIDS = {"Bacteria": 2, "Archaea": 2157}
OUT_COLUMNS = ["mag_id", "gtdb_lineage", "ncbi_lineage", "method", "note", "ncbi_lineage_full",
               "lowest_rank", "lowest_name", "ncbi_taxid", "mapping_confidence", "cross_check"]


# --------------------------------------------------------------------------
# lineage helpers
# --------------------------------------------------------------------------
def split_lineage(lineage: str) -> List[str]:
    """Return the 7 rank names ('' when empty) of a prefixed lineage string."""
    names = [""] * 7
    for part in (lineage or "").split(";"):
        part = part.strip()
        m = re.match(r"^([dpcofgs])__(.*)$", part)
        if m:
            names[RANK_LETTERS.index(m.group(1))] = m.group(2).strip()
    return names


def join_lineage(names: Sequence[str]) -> str:
    return ";".join(f"{l}__{n}" for l, n in zip(RANK_LETTERS, names))


def truncate(names: Sequence[str], max_rank: str) -> List[str]:
    keep = RANK_NAMES.index(max_rank)
    return [n if i <= keep else "" for i, n in enumerate(names)]


def lowest(names: Sequence[str]) -> Tuple[str, str]:
    for i in range(6, -1, -1):
        if names[i]:
            return RANK_NAMES[i], names[i]
    return "", ""


def gtdb_domain(classification: str) -> str:
    if "Archaea" in classification:
        return "Archaea"
    if "Bacteria" in classification:
        return "Bacteria"
    return ""


def read_gtdbtk_summaries(paths: Sequence[str]) -> Dict[str, str]:
    """mag_id -> GTDB classification from one or more GTDB-Tk summary files.

    Accepts the per-domain files or one concatenated file with a single
    header (repeated header lines are skipped).
    """
    out: Dict[str, str] = {}
    for p in paths:
        for r in read_tsv(p, required=("user_genome", "classification")):
            if r["user_genome"] == "user_genome":
                continue
            if r["user_genome"] in out and out[r["user_genome"]] != r["classification"]:
                raise EnaError(f"{r['user_genome']} classified twice with different lineages")
            out[r["user_genome"]] = r["classification"]
    if not out:
        raise EnaError("no genomes in the GTDB-Tk summary")
    return out


def finalise(gtdb: str, ncbi_names: List[str], max_rank: str, notes: List[str]) -> Tuple[List[str], List[str]]:
    dom = gtdb_domain(gtdb)
    full = list(ncbi_names)
    if not full[0] or full[0].lower().startswith("unclassified"):
        if dom:
            full[0] = dom
            notes.append("domain taken from GTDB")
        else:
            notes.append("no domain")
    elif dom and full[0] != dom:
        notes.append(f"NCBI domain {full[0]} differs from GTDB domain {dom}")
    return truncate(full, max_rank), full


# --------------------------------------------------------------------------
# majority vote (GTDB-Tk script)
# --------------------------------------------------------------------------
def run_majority_vote(gtdbtk_dir, out_file, bac120_metadata=None, ar53_metadata=None,
                      script: str = "gtdb_to_ncbi_majority_vote.py", prefix: str = "gtdbtk") -> Path:
    exe = shutil.which(script) or (script if Path(script).exists() else None)
    if exe is None:
        raise EnaError(f"{script} not found; run it in the GTDB-Tk environment and pass --mv-table")
    if not (bac120_metadata or ar53_metadata):
        raise EnaError("--run-mv needs --bac120-metadata and/or --ar53-metadata (GTDB release metadata)")
    cmd = [exe, "--gtdbtk_output_dir", str(gtdbtk_dir), "--output_file", str(out_file), "--gtdbtk_prefix", prefix]
    if bac120_metadata:
        cmd += ["--bac120_metadata_file", str(bac120_metadata)]
    if ar53_metadata:
        cmd += ["--ar53_metadata_file", str(ar53_metadata)]
    LOG.info("running %s", " ".join(cmd))
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not Path(out_file).exists():
        raise EnaError(f"{script} failed ({res.returncode}):\n{res.stdout[-2000:]}\n{res.stderr[-2000:]}")
    return Path(out_file)


def read_mv_table(path) -> Dict[str, Tuple[str, str]]:
    rows = read_tsv(path, required=("Genome ID", "GTDB classification", "Majority vote NCBI classification"))
    return {r["Genome ID"]: (r["GTDB classification"], r["Majority vote NCBI classification"]) for r in rows}


def majority_vote_lineages(gtdb: Dict[str, str], mv: Dict[str, Tuple[str, str]], max_rank: str) -> List[Dict]:
    rows = []
    for mag in sorted(gtdb, key=_natural_key):
        notes: List[str] = []
        if mag in mv:
            raw = mv[mag][1]
            names = split_lineage(raw) if "__" in raw else [""] * 7
            if raw.startswith("Unclassified"):
                notes.append("GTDB-Tk unclassified")
            elif not any(names):
                notes.append("no NCBI majority-vote taxon")
        else:
            names = [""] * 7
            notes.append("genome missing from majority-vote table")
        names_t, full = finalise(gtdb[mag], names, max_rank, notes)
        r, n = lowest(names_t)
        rows.append({"mag_id": mag, "gtdb_lineage": gtdb[mag], "ncbi_lineage": join_lineage(names_t) if names_t[0] else "",
                     "method": "majority_vote", "note": "; ".join(notes), "ncbi_lineage_full": join_lineage(full),
                     "lowest_rank": r, "lowest_name": n, "ncbi_taxid": "", "mapping_confidence": ""})
    return rows


# --------------------------------------------------------------------------
# minimal xlsx reader (stdlib) for the GTDB comparison workbooks
# --------------------------------------------------------------------------
_NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
       "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships"}


def _col_index(ref: str) -> int:
    letters = re.match(r"[A-Z]+", ref).group(0)
    n = 0
    for ch in letters:
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_xlsx(path) -> Dict[str, List[List[str]]]:
    """Return {sheet name: rows (lists of cell strings)} using only the stdlib."""
    m = "{%s}" % _NS["m"]
    with zipfile.ZipFile(path) as z:
        shared: List[str] = []
        if "xl/sharedStrings.xml" in z.namelist():
            for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", _NS):
                shared.append("".join(t.text or "" for t in si.iter(m + "t")))
        wb = ET.fromstring(z.read("xl/workbook.xml"))
        rels = ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))
        relmap = {r.attrib["Id"]: r.attrib["Target"] for r in rels}
        sheets: Dict[str, List[List[str]]] = {}
        for s in wb.find("m:sheets", _NS):
            target = relmap[s.attrib["{%s}id" % _NS["r"]]].lstrip("/")
            if not target.startswith("xl/"):
                target = "xl/" + target
            rows: List[List[str]] = []
            for row in ET.fromstring(z.read(target)).iter(m + "row"):
                cells: Dict[int, str] = {}
                for c in row.findall("m:c", _NS):
                    t = c.attrib.get("t")
                    v = c.find("m:v", _NS)
                    if t == "s" and v is not None:
                        val = shared[int(v.text)]
                    elif t == "inlineStr":
                        val = "".join(x.text or "" for x in c.iter(m + "t"))
                    else:
                        val = v.text if (v is not None and v.text is not None) else ""
                    cells[_col_index(c.attrib.get("r", "A1"))] = val
                if cells:
                    rows.append([cells.get(i, "") for i in range(max(cells) + 1)])
            sheets[s.attrib["name"]] = rows
    return sheets


# --------------------------------------------------------------------------
# mapping-table method
# --------------------------------------------------------------------------
_ENTRY = re.compile(r"([dpcofgs])__(.*?)\s+(\d+(?:\.\d+)?)%")


def parse_ncbi_list(text: str) -> List[Tuple[str, float]]:
    """'g__Bacillus 90.1%, g__(f__Bacillaceae) 9.9%' -> [('Bacillus', 90.1)] (named taxa only)."""
    out = []
    for m in _ENTRY.finditer(text or ""):
        name = re.sub(r"\([^)]*\)", "", m.group(2)).strip()
        if name:
            out.append((name, float(m.group(3))))
    return sorted(out, key=lambda x: -x[1])


def load_mapping_tables(paths: Sequence[str], min_percent: float = 50.0) -> Dict[str, Tuple[str, float]]:
    """GTDB taxon ('g__X') -> (best NCBI name, percent) from GTDB comparison workbooks or TSVs."""
    lookup: Dict[str, Tuple[str, float]] = {}
    for p in paths:
        if str(p).endswith(".xlsx"):
            sheets = read_xlsx(p)
        else:  # TSV with the same two columns as a workbook sheet
            rows = [line.rstrip("\n").split("\t") for line in open(p)]
            sheets = {Path(p).stem: rows}
        for sname, rows in sheets.items():
            if not rows:
                continue
            header = [h.strip() for h in rows[0]]
            gcol = next((i for i, h in enumerate(header) if re.match(r"^GTDB R\d+ \w+$", h)), None)
            ncol = next((i for i, h in enumerate(header) if h.startswith("List of NCBI")), None)
            if gcol is None or ncol is None:
                LOG.warning("%s [%s]: no 'GTDB Rxxx <rank>' / 'List of NCBI ...' columns, skipped", p, sname)
                continue
            rank = header[gcol].split()[-1].lower()
            letter = SHEET_RANKS.get(rank)
            wrong = 0
            for r in rows[1:]:
                if len(r) <= max(gcol, ncol) or not r[gcol]:
                    continue
                key = r[gcol].strip()
                if letter and not key.startswith(letter + "__"):
                    wrong += 1
                    continue
                best = parse_ncbi_list(r[ncol])
                if best and best[0][1] >= min_percent:
                    lookup[key] = best[0]
            if wrong:
                LOG.warning("%s [%s]: %d rows do not carry the %s__ prefix expected for rank %s; skipped "
                            "(sheet content does not match its header)", p, sname, wrong, letter, rank)
    LOG.info("mapping lookup: %d GTDB taxa with an NCBI taxon >= %.0f%%", len(lookup), min_percent)
    return lookup


class Taxdump:
    """NCBI taxdump reader (directory, .tar.gz/.tgz or .zip) without taxonkit.

    nodes.dmp is held as compact arrays; names are loaded only for the names
    queried and the taxa on their lineages (two passes over names.dmp).
    """

    NAME_CLASSES_PRIMARY = ("scientific name",)
    NAME_CLASSES_SECONDARY = ("synonym", "equivalent name", "includes", "genbank synonym")

    def __init__(self, source):
        self.source = str(source)
        self.parent = array("i")
        self.rank = array("b")
        self._ranks: List[str] = []
        self._rank_idx: Dict[str, int] = {}
        self._load_nodes()

    # -- raw file access ------------------------------------------------
    def _open(self, member: str):
        src = self.source
        if os.path.isdir(src):
            return open(os.path.join(src, member), "r", encoding="utf-8", errors="replace")
        if src.endswith(".zip"):
            z = zipfile.ZipFile(src)
            name = next(n for n in z.namelist() if os.path.basename(n) == member)
            return io.TextIOWrapper(z.open(name), encoding="utf-8", errors="replace")
        if src.endswith((".tar.gz", ".tgz", ".tar")):
            t = tarfile.open(src)
            m = next(x for x in t.getmembers() if os.path.basename(x.name) == member)
            return io.TextIOWrapper(t.extractfile(m), encoding="utf-8", errors="replace")
        raise EnaError(f"taxdump: {src} is not a directory, .zip or .tar.gz")

    @staticmethod
    def _fields(line: str) -> List[str]:
        return [f.strip() for f in line.rstrip("\t|\n").split("\t|\t")]

    def _load_nodes(self):
        ids, parents, ranks = array("i"), array("i"), array("b")
        with self._open("nodes.dmp") as fh:
            for line in fh:
                f = self._fields(line)
                r = f[2]
                if r not in self._rank_idx:
                    self._rank_idx[r] = len(self._ranks)
                    self._ranks.append(r)
                ids.append(int(f[0]))
                parents.append(int(f[1]))
                ranks.append(self._rank_idx[r])
        n = max(ids) + 1
        self.parent = array("i", [0]) * n
        self.rank = array("b", [-1]) * n
        for i, p, r in zip(ids, parents, ranks):
            self.parent[i] = p
            self.rank[i] = r
        LOG.info("taxdump %s: %d nodes", self.source, len(ids))

    def rank_of(self, taxid: int) -> str:
        r = self.rank[taxid] if 0 < taxid < len(self.rank) else -1
        return self._ranks[r] if r >= 0 else ""

    def ancestors(self, taxid: int) -> List[int]:
        out = []
        seen = set()
        while 0 < taxid < len(self.parent) and taxid not in seen:
            out.append(taxid)
            seen.add(taxid)
            if self.parent[taxid] == taxid:
                break
            taxid = self.parent[taxid]
        return out

    def resolve(self, names: Iterable[str]) -> Tuple[Dict[str, List[Tuple[int, str]]], Dict[int, str]]:
        """Return name -> [(taxid, name class)] and taxid -> scientific name for lineages."""
        wanted = set(names)
        hits: Dict[str, List[Tuple[int, str]]] = {n: [] for n in wanted}
        with self._open("names.dmp") as fh:
            for line in fh:
                f = self._fields(line)
                if len(f) < 4:
                    continue
                if f[1] in wanted and (f[3] in self.NAME_CLASSES_PRIMARY or f[3] in self.NAME_CLASSES_SECONDARY):
                    hits[f[1]].append((int(f[0]), f[3]))
        need: Set[int] = set()
        for lst in hits.values():
            for tid, _c in lst:
                need.update(self.ancestors(tid))
        sci: Dict[int, str] = {}
        with self._open("names.dmp") as fh:
            for line in fh:
                f = self._fields(line)
                if len(f) >= 4 and f[3] == "scientific name":
                    tid = int(f[0])
                    if tid in need:
                        sci[tid] = f[1]
        return hits, sci

    def lineage_names(self, taxid: int, sci: Dict[int, str]) -> List[str]:
        names = [""] * 7
        for t in self.ancestors(taxid):
            letter = NCBI_RANK_TO_LETTER.get(self.rank_of(t))
            if letter and not names[RANK_LETTERS.index(letter)]:
                names[RANK_LETTERS.index(letter)] = sci.get(t, "")
        return names


def choose_taxid(cands: List[Tuple[int, str]], taxdump: Taxdump, domain: str, rank_letter: Optional[str]) -> Optional[int]:
    """Pick a taxid for a name: within the GTDB domain, preferring scientific names and the expected rank."""
    dom_id = DOMAIN_TAXIDS.get(domain)
    scored = []
    for tid, cls in cands:
        anc = taxdump.ancestors(tid)
        if dom_id and dom_id not in anc:
            continue
        rank_ok = rank_letter is None or NCBI_RANK_TO_LETTER.get(taxdump.rank_of(tid)) == rank_letter
        scored.append(((cls != "scientific name"), (not rank_ok), tid))
    if not scored:
        return None
    scored.sort()
    return scored[0][2]


def mapping_table_lineages(gtdb: Dict[str, str], lookup: Dict[str, Tuple[str, float]], taxdump: Taxdump,
                           max_rank: str) -> List[Dict]:
    # names to resolve: mapped NCBI names, plain GTDB names (+ 'Candidatus ' variants), domains
    queries: Set[str] = set(DOMAIN_TAXIDS)
    for cls in gtdb.values():
        names = split_lineage(cls)
        for i, n in enumerate(names):
            if not n:
                continue
            key = f"{RANK_LETTERS[i]}__{n}"
            if key in lookup:
                queries.add(lookup[key][0])
            queries.update({n, "Candidatus " + n})
    hits, sci = taxdump.resolve(queries)
    rows = []
    for mag in sorted(gtdb, key=_natural_key):
        cls = gtdb[mag]
        dom = gtdb_domain(cls)
        names = split_lineage(cls)
        notes: List[str] = []
        result = None
        for i in range(6, 0, -1):  # species .. phylum
            n = names[i]
            if not n:
                continue
            key = f"{RANK_LETTERS[i]}__{n}"
            if key in lookup:
                ncbi_name, pct = lookup[key]
                tid = choose_taxid(hits.get(ncbi_name, []), taxdump, dom, RANK_LETTERS[i])
                if tid:
                    result = (tid, f"mapped {key} -> {ncbi_name} ({pct:.1f}%)", pct)
                    break
            for cand in (n, "Candidatus " + n):
                tid = choose_taxid(hits.get(cand, []), taxdump, dom, RANK_LETTERS[i])
                if tid:
                    result = (tid, f"direct match {RANK_LETTERS[i]}__{cand}", 100.0)
                    break
            if result:
                break
        if result is None and dom:
            tid = DOMAIN_TAXIDS[dom]
            result = (tid, "fallback to domain", 100.0)
        if result is None:
            rows.append({"mag_id": mag, "gtdb_lineage": cls, "ncbi_lineage": "", "method": "mapping_table",
                         "note": "unmapped (no domain)", "ncbi_lineage_full": "", "lowest_rank": "",
                         "lowest_name": "", "ncbi_taxid": "", "mapping_confidence": ""})
            continue
        tid, how, pct = result
        notes.append(how)
        full = taxdump.lineage_names(tid, sci)
        names_t, full = finalise(cls, full, max_rank, notes)
        r, n = lowest(names_t)
        rows.append({"mag_id": mag, "gtdb_lineage": cls, "ncbi_lineage": join_lineage(names_t) if names_t[0] else "",
                     "method": "mapping_table", "note": "; ".join(notes), "ncbi_lineage_full": join_lineage(full),
                     "lowest_rank": r, "lowest_name": n, "ncbi_taxid": str(tid),
                     "mapping_confidence": f"{pct:.2f}"})
    return rows


def check_names_in_taxdump(rows: List[Dict], taxdump: Taxdump) -> None:
    """Annotate rows whose NCBI names are absent from the taxdump (outdated names)."""
    names = {n for r in rows for n in split_lineage(r["ncbi_lineage"]) if n}
    hits, _sci = taxdump.resolve(names)
    for r in rows:
        missing = [f"{l}__{n}" for l, n in zip(RANK_LETTERS, split_lineage(r["ncbi_lineage"]))
                   if n and not any(c == "scientific name" for _t, c in hits.get(n, []))]
        if missing:
            r["note"] = "; ".join(x for x in (r["note"], "not a current NCBI scientific name: " + ",".join(missing)) if x)


# --------------------------------------------------------------------------
# cross-check / consensus of the two methods
# --------------------------------------------------------------------------
def _norm(name: str) -> str:
    return re.sub(r"^Candidatus\s+", "", name or "").strip().lower()


def lineage_conflicts(a: Sequence[str], b: Sequence[str]) -> List[int]:
    """Rank indices (domain..genus) where both lineages are named but differ ('Candidatus' ignored)."""
    return [i for i in range(6) if a[i] and b[i] and _norm(a[i]) != _norm(b[i])]


def consensus_names(mv: Sequence[str], table: Sequence[str]) -> List[str]:
    out = [""] * 7
    for i in range(7):
        if mv[i] and table[i] and _norm(mv[i]) == _norm(table[i]):
            out[i] = table[i]  # current NCBI spelling (taxdump)
        else:
            break
    return out


def combine_methods(mv_rows: List[Dict], table_rows: List[Dict], max_rank: str, consensus: bool) -> List[Dict]:
    tab = {r["mag_id"]: r for r in table_rows}
    out = []
    for r in mv_rows:
        t = tab.get(r["mag_id"])
        mvf = split_lineage(r["ncbi_lineage_full"])
        tf = split_lineage(t["ncbi_lineage_full"]) if t else [""] * 7
        conf = lineage_conflicts(mvf, tf)
        cc = "; ".join(f"{RANK_NAMES[i]}: majority_vote {mvf[i]} vs mapping_table {tf[i]}" for i in conf)
        row = dict(r)
        row["cross_check"] = cc or "agree"
        if consensus:
            notes = [n for n in (r["note"],) if n]
            names = consensus_names(mvf, tf)
            names_t, full = finalise(r["gtdb_lineage"], names, max_rank, notes)
            if conf:
                notes.append(f"consensus stops above {RANK_NAMES[conf[0]]}")
            rk, nm = lowest(names_t)
            row.update({"ncbi_lineage": join_lineage(names_t) if names_t[0] else "", "method": "consensus",
                        "note": "; ".join(notes), "ncbi_lineage_full": join_lineage(full), "lowest_rank": rk,
                        "lowest_name": nm})
        out.append(row)
    return out


def _natural_key(s: str):
    return [int(t) if t.isdigit() else t for t in re.split(r"(\d+)", s)]


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtdbtk", nargs="+", required=True, help="GTDB-Tk summary TSV(s) (bac120/ar53 or concatenated)")
    ap.add_argument("--out", required=True)
    ap.add_argument("--method", choices=("majority_vote", "mapping_table", "consensus"), default="majority_vote")
    ap.add_argument("--cross-check", action="store_true",
                    help="majority_vote: also run mapping_table and report disagreements (column cross_check)")
    ap.add_argument("--max-rank", choices=RANK_NAMES[1:], default="genus",
                    help="empty all ranks below this one (default genus)")
    mv = ap.add_argument_group("majority_vote")
    mv.add_argument("--mv-table", help="output of gtdb_to_ncbi_majority_vote.py")
    mv.add_argument("--run-mv", action="store_true", help="run gtdb_to_ncbi_majority_vote.py (needs GTDB-Tk)")
    mv.add_argument("--gtdbtk-dir", help="GTDB-Tk classify_wf --out_dir (with classify/ trees)")
    mv.add_argument("--gtdbtk-prefix", default="gtdbtk")
    mv.add_argument("--bac120-metadata", help="GTDB bac120_metadata_rXXX.tsv(.gz) of the GTDB-Tk reference release")
    mv.add_argument("--ar53-metadata", help="GTDB ar53_metadata_rXXX.tsv(.gz)")
    mv.add_argument("--mv-script", default="gtdb_to_ncbi_majority_vote.py")
    mt = ap.add_argument_group("mapping_table")
    mt.add_argument("--mapping-table", nargs="+", help="GTDB gtdb_vs_ncbi workbook(s) (.xlsx) or TSV sheet exports")
    mt.add_argument("--min-percent", type=float, default=50.0)
    ap.add_argument("--taxdump", help="NCBI taxdump (dir, .zip or .tar.gz): required for mapping_table, "
                                      "optional name check for majority_vote")
    ap.add_argument("-v", "--verbose", action="count", default=0)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        gtdb = read_gtdbtk_summaries(args.gtdbtk)
        taxdump = Taxdump(args.taxdump) if args.taxdump else None
        mv_rows = table_rows = None
        if args.method in ("majority_vote", "consensus"):
            if args.run_mv:
                if not args.gtdbtk_dir:
                    raise EnaError("--run-mv needs --gtdbtk-dir")
                out_mv = Path(args.out).with_name(Path(args.out).stem + ".majority_vote_raw.tsv")
                mv_path = run_majority_vote(args.gtdbtk_dir, out_mv, args.bac120_metadata, args.ar53_metadata,
                                            args.mv_script, args.gtdbtk_prefix)
            elif args.mv_table:
                mv_path = args.mv_table
            else:
                raise EnaError(f"{args.method} needs --mv-table or --run-mv")
            mv_rows = majority_vote_lineages(gtdb, read_mv_table(mv_path), args.max_rank)
            if taxdump:
                check_names_in_taxdump(mv_rows, taxdump)
        if args.method in ("mapping_table", "consensus") or args.cross_check:
            if not (args.mapping_table and taxdump):
                raise EnaError("mapping_table, consensus and --cross-check need --mapping-table and --taxdump")
            table_rows = mapping_table_lineages(gtdb, load_mapping_tables(args.mapping_table, args.min_percent),
                                                taxdump, args.max_rank)
        if args.method == "mapping_table":
            rows = table_rows
        elif args.method == "consensus" or args.cross_check:
            rows = combine_methods(mv_rows, table_rows, args.max_rank, consensus=(args.method == "consensus"))
            bad = [r for r in rows if r["cross_check"] != "agree"]
            if bad:
                LOG.warning("%d genomes where majority_vote and mapping_table disagree (see column cross_check)", len(bad))
        else:
            rows = mv_rows
    except EnaError as exc:
        LOG.error("%s", exc)
        return 1
    write_tsv(args.out, rows, OUT_COLUMNS)
    ranks = {}
    for r in rows:
        ranks[r["lowest_rank"] or "none"] = ranks.get(r["lowest_rank"] or "none", 0) + 1
    LOG.info("%d genomes -> %s (lowest rank counts: %s)", len(rows), args.out, ranks)
    return 0 if all(r["ncbi_lineage"] for r in rows) else 2


if __name__ == "__main__":
    sys.exit(main())
