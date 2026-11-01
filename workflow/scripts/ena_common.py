#!/usr/bin/env python3
"""Shared helpers for the ENA submission scripts of bass-fastq-to-mag.

Everything here is standard library only, so the ENA scripts run in any
Python >= 3.9 environment (the workflow env, the user's ``ena_tools`` env,
or a plain interpreter for the tests).

Contents
--------
* FASTA streaming (plain / .gz / .bz2), ENA sequence rules, validated gzip copy
* ENA name and accession rules (assembly names, run / sample / study accessions)
* samples.tsv parsing (columns: sample, fq1, fq2, assembly_group,
  optional run_accession, sample_accession)
* config loading (YAML ``ena:`` block or JSON)
* Webin credentials from environment variables only, plus a log redactor
"""
from __future__ import annotations

import bz2
import csv
import gzip
import hashlib
import io
import json
import logging
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Sequence, Tuple

LOG = logging.getLogger("bfm.ena")

# --------------------------------------------------------------------------
# ENA rules and vocabularies
# --------------------------------------------------------------------------
RUN_ACCESSION_RE = re.compile(r"^[EDS]RR\d{6,}$")
ASSEMBLY_ACCESSION_RE = re.compile(r"^[EDS]RZ\d{6,}$")
SAMPLE_ACCESSION_RE = re.compile(r"^(?:[EDS]RS\d{6,}|SAM[END][A-Z]?\d+)$")
STUDY_ACCESSION_RE = re.compile(r"^(?:PRJ[EDN][A-Z]\d+|[EDS]RP\d{6,})$")
# ENA assembly name rule (webin-cli genome context)
ASSEMBLY_NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9 _#\-.]*$")
ASSEMBLY_NAME_MAX = 50
# genome_uploader (test mode) appends "_" + a 10-digit unix timestamp to each
# genome name; the result is used as sample alias and ASSEMBLYNAME.
TEST_SUFFIX_LEN = 11
# ENA accepted sequence characters (case-insensitive)
VALID_BASES = frozenset("ACGTUBDHKMNRSVWY")
MIN_SEQ_LEN = 20
TRUE_VALUES = {"true", "yes", "y", "1"}
FALSE_VALUES = {"false", "no", "n", "0"}
# INSDC missing-value term for a sample that represents a group of samples
# (co-assemblies); see ENA docs "Reporting Missing Values".
MISSING_SAMPLE_GROUP = "missing: sample group"

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fas", ".fsa")
COMPRESSION_SUFFIXES = (".gz", ".bz2")


class EnaError(Exception):
    """Base class for user-facing errors of the ENA scripts."""


class ValidationError(EnaError):
    """Raised when inputs violate ENA or genome_uploader rules."""

    def __init__(self, problems: Sequence[str], header: str = "Validation failed"):
        self.problems = list(problems)
        msg = header + ":\n  - " + "\n  - ".join(self.problems)
        super().__init__(msg)


class CredentialsError(EnaError):
    """Raised when Webin credentials are missing (message never holds secrets)."""


def setup_logging(verbose: int = 0) -> None:
    level = logging.WARNING - 10 * min(verbose, 2) if verbose else logging.INFO
    logging.basicConfig(level=level, format="%(asctime)s %(levelname)s %(name)s: %(message)s",
                        stream=sys.stderr)


# --------------------------------------------------------------------------
# small utilities
# --------------------------------------------------------------------------
def parse_bool(value, field: str = "value") -> bool:
    s = str(value).strip().lower()
    if s in TRUE_VALUES:
        return True
    if s in FALSE_VALUES or s in ("", "nan", "none"):
        return False
    raise ValueError(f"{field}: cannot interpret {value!r} as True/False")


def fmt_bool(value: bool) -> str:
    return "True" if value else "False"


def fmt_float(value: float, ndigits: int = 2) -> str:
    """Format a float with at most ``ndigits`` decimals (no trailing zeros)."""
    s = f"{float(value):.{ndigits}f}"
    if "." in s:
        s = s.rstrip("0").rstrip(".")
    return s if s not in ("-0", "") else "0"


def md5_file(path, chunk: int = 1 << 20) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(chunk), b""):
            h.update(block)
    return h.hexdigest()


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")


def atomic_write_text(path, text: str) -> None:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def strip_suffixes(name: str, suffixes: Iterable[str]) -> str:
    changed = True
    while changed:
        changed = False
        for suf in suffixes:
            if name.endswith(suf) and len(name) > len(suf):
                name = name[: -len(suf)]
                changed = True
    return name


def genome_id_from_path(path) -> str:
    """'dir/north_MAG_1.fna.gz' -> 'north_MAG_1'."""
    return strip_suffixes(Path(path).name, COMPRESSION_SUFFIXES + FASTA_SUFFIXES)


def list_fasta_files(directory, extensions: Optional[Sequence[str]] = None) -> List[Path]:
    directory = Path(directory)
    if not directory.is_dir():
        raise EnaError(f"not a directory: {directory}")
    exts = tuple(extensions) if extensions else tuple(
        s + c for s in FASTA_SUFFIXES for c in ("",) + COMPRESSION_SUFFIXES)
    files = sorted(p for p in directory.iterdir() if p.is_file() and p.name.endswith(exts))
    return files


# --------------------------------------------------------------------------
# TSV helpers (csv module; no pandas needed)
# --------------------------------------------------------------------------
def read_tsv(path, required: Sequence[str] = ()) -> List[Dict[str, str]]:
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if reader.fieldnames is None:
            raise EnaError(f"{path}: empty file")
        fields = [f.strip() if f else f for f in reader.fieldnames]
        reader.fieldnames = fields
        missing = [c for c in required if c not in fields]
        if missing:
            raise EnaError(f"{path}: missing column(s) {', '.join(missing)}; found {', '.join(fields)}")
        rows = []
        for row in reader:
            rows.append({k: (v.strip() if isinstance(v, str) else ("" if v is None else v))
                         for k, v in row.items() if k is not None})
        return rows


def write_tsv(path, rows: Sequence[Dict[str, object]], columns: Sequence[str]) -> None:
    buf = io.StringIO()
    writer = csv.writer(buf, delimiter="\t", lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
    writer.writerow(columns)
    for row in rows:
        writer.writerow(["" if row.get(c) is None else row.get(c) for c in columns])
    atomic_write_text(path, buf.getvalue())


# --------------------------------------------------------------------------
# FASTA
# --------------------------------------------------------------------------
def open_text(path):
    path = str(path)
    if path.endswith(".gz"):
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", newline="")
    if path.endswith(".bz2"):
        return io.TextIOWrapper(bz2.open(path, "rb"), encoding="utf-8", newline="")
    return open(path, "r", encoding="utf-8", newline="")


def iter_fasta(path) -> Iterator[Tuple[str, str, str]]:
    """Yield (sequence_id, full_header, sequence) for each record.

    The sequence id is the first whitespace-delimited token of the header,
    which is how webin-cli, CoverM and barrnap name sequences.
    """
    header = None
    chunks: List[str] = []
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header.split()[0] if header.split() else "", header, "".join(chunks)
                header = line[1:].strip()
                chunks = []
            else:
                if header is None:
                    raise EnaError(f"{path}: sequence data before the first header")
                chunks.append(line.strip())
        if header is not None:
            yield header.split()[0] if header.split() else "", header, "".join(chunks)


def fasta_lengths(path) -> Dict[str, int]:
    return {sid: len(seq) for sid, _h, seq in iter_fasta(path)}


def check_sequence(seq_id: str, seq: str) -> List[str]:
    """ENA sequence rules: >= 20 bp, no terminal N, IUPAC characters only."""
    problems = []
    if not seq_id:
        problems.append("empty sequence name")
    s = seq.upper()
    if len(s) < MIN_SEQ_LEN:
        problems.append(f"{seq_id}: shorter than {MIN_SEQ_LEN} bp ({len(s)} bp)")
    if s.startswith("N") or s.endswith("N"):
        problems.append(f"{seq_id}: terminal N")
    bad = set(s) - VALID_BASES
    if bad:
        problems.append(f"{seq_id}: invalid character(s) {''.join(sorted(bad))}")
    return problems


def write_validated_fasta_gz(src, dst, force: bool = False, max_problems: int = 20) -> Dict[str, object]:
    """Stream ``src`` (plain/gz/bz2) into gzipped ``dst`` while checking ENA rules.

    Returns a stats dict (n_seqs, total_bp, md5, problems, reused). Raises
    ValidationError if any record violates the ENA sequence rules, in which
    case ``dst`` is not written. If ``dst`` exists and is newer than ``src``
    (and ``force`` is False) the file is re-validated but not rewritten.
    """
    src, dst = Path(src), Path(dst)
    if not src.exists():
        raise EnaError(f"FASTA not found: {src}")
    reuse = (not force and dst.exists() and dst.stat().st_size > 0
             and dst.stat().st_mtime >= src.stat().st_mtime)
    problems: List[str] = []
    seen = set()
    n_seqs = 0
    total = 0
    dst.parent.mkdir(parents=True, exist_ok=True)
    tmp = dst.with_name(dst.name + ".tmp")
    out = None if reuse else gzip.open(tmp, "wt", compresslevel=6, encoding="utf-8")
    try:
        for sid, header, seq in iter_fasta(dst if reuse else src):
            n_seqs += 1
            total += len(seq)
            if sid in seen:
                problems.append(f"{sid}: duplicated sequence name")
            seen.add(sid)
            problems.extend(check_sequence(sid, seq))
            if out is not None:
                out.write(">" + header + "\n")
                for i in range(0, len(seq), 80):
                    out.write(seq[i:i + 80] + "\n")
        if n_seqs == 0:
            problems.append("no sequences")
    finally:
        if out is not None:
            out.close()
    if problems:
        if tmp.exists():
            tmp.unlink()
        shown = problems[:max_problems]
        if len(problems) > max_problems:
            shown.append(f"... {len(problems) - max_problems} more")
        raise ValidationError(shown, header=f"{src}: FASTA violates ENA rules")
    if not reuse:
        os.replace(tmp, dst)
    return {"n_seqs": n_seqs, "total_bp": total, "md5": md5_file(dst), "reused": reuse}


# --------------------------------------------------------------------------
# names and accessions
# --------------------------------------------------------------------------
def validate_assembly_name(name: str, reserve_test_suffix: bool = False) -> List[str]:
    problems = []
    limit = ASSEMBLY_NAME_MAX - (TEST_SUFFIX_LEN if reserve_test_suffix else 0)
    if not ASSEMBLY_NAME_RE.match(name or ""):
        problems.append(f"{name!r}: assembly names must match {ASSEMBLY_NAME_RE.pattern}")
    if len(name) > limit:
        why = " (genome_uploader appends an 11-character test suffix)" if reserve_test_suffix else ""
        problems.append(f"{name!r}: longer than {limit} characters{why}")
    return problems


def split_accessions(cell) -> List[str]:
    if cell is None:
        return []
    return [a for a in re.split(r"[,;\s]+", str(cell).strip()) if a and a.lower() != "nan"]


def dedupe(seq: Iterable[str]) -> List[str]:
    seen, out = set(), []
    for x in seq:
        if x not in seen:
            seen.add(x)
            out.append(x)
    return out


# --------------------------------------------------------------------------
# samples.tsv
# --------------------------------------------------------------------------
def load_samples(path) -> List[Dict[str, str]]:
    rows = read_tsv(path, required=("sample",))
    names = [r["sample"] for r in rows]
    dup = sorted({n for n in names if names.count(n) > 1})
    if dup:
        raise EnaError(f"{path}: duplicated sample name(s): {', '.join(dup)}")
    if not rows:
        raise EnaError(f"{path}: no samples")
    return rows


def select_group(samples: List[Dict[str, str]], assembly_group: Optional[str] = None) -> List[Dict[str, str]]:
    groups = dedupe(r.get("assembly_group", "") or "all" for r in samples)
    if assembly_group is None:
        if len(groups) > 1:
            raise EnaError("samples.tsv defines several assembly groups "
                           f"({', '.join(groups)}); pass --assembly-group")
        return list(samples)
    chosen = [r for r in samples if (r.get("assembly_group", "") or "all") == assembly_group]
    if not chosen:
        raise EnaError(f"assembly group {assembly_group!r} not found; available: {', '.join(groups)}")
    return chosen


def run_accessions(samples: List[Dict[str, str]]) -> List[str]:
    missing = [r["sample"] for r in samples if not split_accessions(r.get("run_accession"))]
    if missing:
        raise EnaError("run_accession missing in samples.tsv for: " + ", ".join(missing))
    runs = dedupe(a for r in samples for a in split_accessions(r.get("run_accession")))
    bad = [a for a in runs if not RUN_ACCESSION_RE.match(a)]
    if bad:
        raise EnaError("malformed run accession(s): " + ", ".join(bad))
    return runs


def sample_accessions(samples: List[Dict[str, str]]) -> List[str]:
    missing = [r["sample"] for r in samples if not split_accessions(r.get("sample_accession"))]
    if missing:
        raise EnaError("sample_accession missing in samples.tsv for: " + ", ".join(missing))
    accs = dedupe(a for r in samples for a in split_accessions(r.get("sample_accession")))
    bad = [a for a in accs if not SAMPLE_ACCESSION_RE.match(a)]
    if bad:
        raise EnaError("malformed sample accession(s): " + ", ".join(bad))
    return accs


# --------------------------------------------------------------------------
# config
# --------------------------------------------------------------------------
def load_config(path, section: Optional[str] = "ena") -> Dict:
    """Load a YAML (needs PyYAML) or JSON config and return ``section``."""
    if path is None:
        return {}
    text = Path(path).read_text()
    if str(path).endswith(".json"):
        data = json.loads(text)
    else:
        try:
            import yaml  # type: ignore
        except ImportError as exc:  # pragma: no cover - depends on env
            raise EnaError("PyYAML is required to read YAML configs (or pass a .json config)") from exc
        data = yaml.safe_load(text) or {}
    if section and isinstance(data, dict) and section in data:
        return data[section] or {}
    return data if isinstance(data, dict) else {}


def cfg_get(cfg: Dict, dotted: str, default=None):
    cur = cfg
    for part in dotted.split("."):
        if not isinstance(cur, dict) or part not in cur:
            return default
        cur = cur[part]
    return default if cur is None else cur


# --------------------------------------------------------------------------
# credentials (environment variables only) and redaction
# --------------------------------------------------------------------------
ENV_USER = "ENA_WEBIN"
ENV_PASSWORD = "ENA_WEBIN_PASSWORD"


def get_webin_credentials(required: bool = True) -> Tuple[Optional[str], Optional[str]]:
    """Return (user, password) from ENA_WEBIN / ENA_WEBIN_PASSWORD.

    Never reads .env files or other sources. Raises CredentialsError (without
    echoing any value) if ``required`` and a variable is missing or empty.
    """
    user = os.environ.get(ENV_USER, "").strip()
    password = os.environ.get(ENV_PASSWORD, "")
    if required and (not user or not password):
        missing = [n for n, v in ((ENV_USER, user), (ENV_PASSWORD, password)) if not v]
        raise CredentialsError("Webin credentials missing: set environment variable(s) "
                               + ", ".join(missing) + " (values are never read from files)")
    return (user or None), (password or None)


class Redactor:
    """Replace secret strings (and password-like CLI args) in text."""

    _PATTERNS = [re.compile(r"(-password(?:Env)?[=\s]+)\S+", re.I),
                 re.compile(r"(ENA_WEBIN_PASSWORD[=:]\s*)\S+")]

    def __init__(self, secrets: Iterable[Optional[str]] = ()):
        self.secrets = sorted({s for s in secrets if s and len(s) >= 3}, key=len, reverse=True)

    def __call__(self, text: str) -> str:
        if not text:
            return text
        for s in self.secrets:
            text = text.replace(s, "***")
        for pat in self._PATTERNS:
            text = pat.sub(lambda m: m.group(1) + "***" if "Env" not in m.group(0) else m.group(0), text)
        return text
