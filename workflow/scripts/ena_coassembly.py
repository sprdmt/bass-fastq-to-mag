#!/usr/bin/env python3
"""Co-assembly sample and primary-metagenome manifest for ENA.

Why
---
webin-cli accepts exactly one SAMPLE per assembly; a comma-separated list
fails with "Unknown sample ERS..,ERS.. or the sample cannot be referenced".
ENA's FAQ ("How do I register samples for co-assemblies?") asks for a new
sample that references every source sample (or run) in its ``sample derived
from`` attribute, as a comma-separated list without spaces (or an
``ERSx-ERSy`` range covering only used accessions). ENA also asks to inform
the helpdesk before submitting a primary assembly co-assembled from raw reads.

Sub-commands
------------
sample-xml     write the co-assembly sample (SAMPLE_SET), a SUBMISSION (ADD) and
               the Webin REST v2 document (WEBIN) for one assembly group
submit-sample  register it through the Webin REST API: TEST server
               (wwwdev.ebi.ac.uk) by default, production only with --live;
               credentials only from ENA_WEBIN / ENA_WEBIN_PASSWORD; the receipt
               is stored and the accession appended to an accession TSV
               (re-running skips groups already registered in that mode)
manifest       write the webin-cli manifest (-context genome, ASSEMBLY_TYPE
               'primary metagenome') with COVERAGE from the CoverM contig table,
               RUN_REF = all runs of the group, and a validated gzipped FASTA

Examples
--------
  ena_coassembly.py sample-xml --samples config/samples.tsv --group g1 \
      --params config/config.yaml --outdir results/ena/coassembly/g1
  ena_coassembly.py submit-sample --xml-dir results/ena/coassembly/g1            # TEST server
  ena_coassembly.py submit-sample --xml-dir results/ena/coassembly/g1 --live     # production
  ena_coassembly.py manifest --samples config/samples.tsv --group g1 --study PRJEB12345 \
      --accessions-tsv results/ena/coassembly/g1/coassembly_sample_accessions.tsv --mode test \
      --coverm results/mapping/g1/coverm_contig_mean.tsv --assembly-fasta results/assembly/g1/contigs.fa \
      --assembly-name g1_coassembly --program "MEGAHIT v1.2.9" --platform ILLUMINA \
      --outdir results/ena/coassembly/g1
"""
from __future__ import annotations

import argparse
import base64
import json
import logging
import os
import re
import sys
import time
import urllib.error
import urllib.request
import uuid
import xml.dom.minidom as minidom
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import (MISSING_SAMPLE_GROUP, STUDY_ACCESSION_RE, CredentialsError, EnaError,  # noqa: E402
                        ValidationError, atomic_write_text, cfg_get, fmt_float, get_webin_credentials,
                        load_config, load_samples, now_iso, read_tsv, run_accessions, sample_accessions,
                        select_group, setup_logging, validate_assembly_name, write_tsv,
                        write_validated_fasta_gz)

LOG = logging.getLogger("bfm.ena_coassembly")

ENDPOINTS = {
    ("v2", False): "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit",
    ("v2", True): "https://www.ebi.ac.uk/ena/submit/webin-v2/submit",
    ("v1", False): "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/",
    ("v1", True): "https://www.ebi.ac.uk/ena/submit/drop-box/submit/",
}
DEFAULT_CHECKLIST = "ERC000011"  # ENA default sample checklist
MANDATORY_ATTRS = ("collection date", "geographic location (country and/or sea)")
ACC_COLUMNS = ["alias_base", "alias", "accession", "biosample", "mode", "receipt", "registered_at"]
_EXISTS_RE = re.compile(r'alias: "([^"]+)"[^:]+accession: "([^"]+)"')


# --------------------------------------------------------------------------
# XML
# --------------------------------------------------------------------------
def _pretty(elem: ET.Element) -> str:
    return minidom.parseString(ET.tostring(elem, encoding="utf-8")).toprettyxml(indent="  ")


def source_accessions(group_samples: List[Dict[str, str]]) -> Tuple[List[str], str]:
    """Sample accessions of the group, else run accessions (ENA accepts either)."""
    try:
        accs = sample_accessions(group_samples)
        kind = "sample"
    except EnaError:
        accs = run_accessions(group_samples)
        kind = "run"
    prefixes = {re.match(r"^[A-Z]+", a).group(0) for a in accs}
    if len(prefixes) > 1:
        raise EnaError("mixed accession types for 'sample derived from': " + ", ".join(sorted(prefixes)))
    return accs, kind


def build_sample_xml(alias: str, title: str, description: str, taxon_id: str, scientific_name: str,
                     derived_from: Sequence[str], attributes: Dict[str, str], checklist: str = DEFAULT_CHECKLIST,
                     centre_name: Optional[str] = None) -> ET.Element:
    if not str(taxon_id).isdigit():
        raise EnaError(f"taxon id must be numeric (metagenome taxon), got {taxon_id!r}")
    sample_set = ET.Element("SAMPLE_SET")
    sample = ET.SubElement(sample_set, "SAMPLE", alias=alias)
    if centre_name:
        sample.set("center_name", centre_name)
    ET.SubElement(sample, "TITLE").text = title
    sn = ET.SubElement(sample, "SAMPLE_NAME")
    ET.SubElement(sn, "TAXON_ID").text = str(taxon_id)
    ET.SubElement(sn, "SCIENTIFIC_NAME").text = scientific_name
    ET.SubElement(sample, "DESCRIPTION").text = description
    attrs = ET.SubElement(sample, "SAMPLE_ATTRIBUTES")
    ordered = [("sample derived from", ",".join(derived_from))]
    ordered += [(k, v) for k, v in attributes.items() if k not in ("sample derived from", "ENA-CHECKLIST")]
    ordered.append(("ENA-CHECKLIST", checklist))
    for tag, value in ordered:
        a = ET.SubElement(attrs, "SAMPLE_ATTRIBUTE")
        ET.SubElement(a, "TAG").text = tag
        ET.SubElement(a, "VALUE").text = str(value)
    return sample_set


def build_submission_xml(centre_name: Optional[str] = None, hold_date: Optional[str] = None) -> ET.Element:
    sub = ET.Element("SUBMISSION")
    if centre_name:
        sub.set("center_name", centre_name)
    actions = ET.SubElement(sub, "ACTIONS")
    ET.SubElement(ET.SubElement(actions, "ACTION"), "ADD")
    if hold_date:
        ET.SubElement(ET.SubElement(actions, "ACTION"), "HOLD", HoldUntilDate=hold_date)
    return sub


def build_webin_document(sample_set: ET.Element, submission: ET.Element) -> ET.Element:
    webin = ET.Element("WEBIN")
    ss = ET.SubElement(webin, "SUBMISSION_SET")
    ss.append(submission)
    webin.append(sample_set)
    return webin


def with_alias(sample_set: ET.Element, alias: str) -> ET.Element:
    clone = ET.fromstring(ET.tostring(sample_set))
    for s in clone.iter("SAMPLE"):
        s.set("alias", alias)
    return clone


# --------------------------------------------------------------------------
# receipt and HTTP
# --------------------------------------------------------------------------
def parse_receipt(text: str) -> Dict[str, object]:
    """Parse a Webin XML receipt -> dict(success, samples{alias: (acc, biosample)}, errors, infos)."""
    try:
        root = ET.fromstring(text)
    except ET.ParseError as exc:
        raise EnaError(f"receipt is not XML: {text[:300]!r}") from exc
    if root.tag != "RECEIPT":
        root = root.find(".//RECEIPT")
        if root is None:
            raise EnaError("no RECEIPT element in response")
    samples = {}
    for s in root.findall("SAMPLE"):
        ext = s.find("EXT_ID")
        samples[s.get("alias", "")] = (s.get("accession", ""), ext.get("accession", "") if ext is not None else "")
    errors = [(e.text or "").strip() for e in root.iter("ERROR")]
    infos = [(e.text or "").strip() for e in root.iter("INFO")]
    existing = {}
    for e in errors:
        for m in _EXISTS_RE.finditer(e):
            existing[m.group(1)] = m.group(2)
    return {"success": root.get("success", "false").lower() == "true", "samples": samples,
            "errors": errors, "infos": infos, "existing": existing}


def http_post(url: str, body: bytes, headers: Dict[str, str], user: str, password: str,
              timeout: int = 300) -> Tuple[int, str]:
    """POST with HTTP basic auth (credentials never logged)."""
    req = urllib.request.Request(url, data=body, method="POST")
    for k, v in headers.items():
        req.add_header(k, v)
    token = base64.b64encode(f"{user}:{password}".encode()).decode()
    req.add_header("Authorization", "Basic " + token)
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            return resp.status, resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        return exc.code, exc.read().decode("utf-8", errors="replace")


def multipart(files: Dict[str, Tuple[str, bytes]]) -> Tuple[bytes, str]:
    boundary = uuid.uuid4().hex
    parts = []
    for field, (fname, content) in files.items():
        parts.append((f'--{boundary}\r\nContent-Disposition: form-data; name="{field}"; filename="{fname}"\r\n'
                      f"Content-Type: application/xml\r\n\r\n").encode() + content + b"\r\n")
    return b"".join(parts) + f"--{boundary}--\r\n".encode(), f"multipart/form-data; boundary={boundary}"


def submit_sample_set(sample_set: ET.Element, submission: ET.Element, live: bool = False, api: str = "v2",
                      poster=None) -> Tuple[str, Dict[str, object]]:
    """Send the sample to ENA (test server unless live=True). Returns (receipt text, parsed receipt)."""
    poster = poster or http_post
    user, password = get_webin_credentials(required=True)
    url = ENDPOINTS[(api, bool(live))]
    LOG.info("registering co-assembly sample on the %s server (%s)", "PRODUCTION" if live else "TEST", url)
    if api == "v2":
        body = ET.tostring(build_webin_document(sample_set, submission), encoding="utf-8")
        status, text = poster(url, body, {"Content-Type": "application/xml", "Accept": "application/xml"},
                              user, password)
    else:
        body, ctype = multipart({"SUBMISSION": ("submission.xml", ET.tostring(submission, encoding="utf-8")),
                                 "SAMPLE": ("sample.xml", ET.tostring(sample_set, encoding="utf-8"))})
        status, text = poster(url, body, {"Content-Type": ctype, "Accept": "application/xml"}, user, password)
    if status == 401 or status == 403:
        raise CredentialsError(f"ENA rejected the Webin credentials (HTTP {status})")
    if status >= 400 and "<RECEIPT" not in text:
        raise EnaError(f"ENA returned HTTP {status}: {text[:500]}")
    return text, parse_receipt(text)


# --------------------------------------------------------------------------
# accession log
# --------------------------------------------------------------------------
def read_accessions(path) -> List[Dict[str, str]]:
    return read_tsv(path, required=ACC_COLUMNS[:5]) if Path(path).exists() else []


def lookup_accession(path, alias_base: str, mode: str) -> Optional[Dict[str, str]]:
    hits = [r for r in read_accessions(path) if r["alias_base"] == alias_base and r["mode"] == mode]
    return hits[-1] if hits else None


def append_accession(path, row: Dict[str, str]) -> None:
    rows = read_accessions(path) + [row]
    write_tsv(path, rows, ACC_COLUMNS)


# --------------------------------------------------------------------------
# manifest
# --------------------------------------------------------------------------
def build_manifest(study: str, sample: str, assembly_name: str, coverage: float, program: str, platform: str,
                   fasta: str, run_refs: Sequence[str], description: str = "", moleculetype: str = "genomic DNA",
                   mingaplength: Optional[int] = None) -> str:
    problems = []
    if not STUDY_ACCESSION_RE.match(study or ""):
        problems.append(f"STUDY {study!r} is not a PRJEB/ERP accession")
    if not sample or "," in sample or " " in sample.strip():
        problems.append(f"SAMPLE must be one sample accession or alias (got {sample!r}); register the "
                        "co-assembly sample first (ena_coassembly.py submit-sample)")
    problems += validate_assembly_name(assembly_name)
    if not (coverage > 0):
        problems.append(f"COVERAGE must be > 0 (got {coverage})")
    if not str(fasta).endswith((".gz", ".bz2")):
        problems.append("FASTA must be gzip or bzip2 compressed")
    for k, v in (("PROGRAM", program), ("PLATFORM", platform)):
        if not v:
            problems.append(f"{k} is empty")
    if problems:
        raise ValidationError(problems, header="primary assembly manifest")
    lines = [("STUDY", study), ("SAMPLE", sample), ("ASSEMBLYNAME", assembly_name),
             ("ASSEMBLY_TYPE", "primary metagenome"), ("COVERAGE", fmt_float(coverage, 2)),
             ("PROGRAM", program), ("PLATFORM", platform), ("MOLECULETYPE", moleculetype)]
    if mingaplength:
        lines.append(("MINGAPLENGTH", str(int(mingaplength))))
    if description:
        lines.append(("DESCRIPTION", description))
    if run_refs:
        lines.append(("RUN_REF", ",".join(run_refs)))
    lines.append(("FASTA", str(fasta)))
    return "".join(f"{k}\t{v}\n" for k, v in lines)


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------
def parse_attributes(values: Sequence[str]) -> Dict[str, str]:
    out = {}
    for v in values or []:
        if "=" not in v:
            raise EnaError(f"--attribute expects TAG=VALUE, got {v!r}")
        k, val = v.split("=", 1)
        out[k.strip()] = val.strip()
    return out


def cmd_sample_xml(args) -> int:
    cfg = load_config(args.params) if args.params else {}
    cc = cfg_get(cfg, "coassembly", {}) or {}
    group_samples = select_group(load_samples(args.samples), args.group)
    accs, kind = source_accessions(group_samples)
    runs = []
    try:
        runs = run_accessions(group_samples)
    except EnaError:
        pass
    study = args.study or cfg_get(cfg, "study", "")
    group = args.group or (group_samples[0].get("assembly_group") or "all")
    alias = args.alias or f"{cc.get('sample_alias_prefix') or (study or 'bfm')}_{group}_coassembly"
    taxon_id = args.taxon_id or cc.get("taxon_id") or cfg_get(cfg, "metagenome_taxid")
    sci = args.scientific_name or cc.get("scientific_name") or cfg_get(cfg, "metagenome")
    if not taxon_id or not sci:
        raise EnaError("set the metagenome taxon: --taxon-id/--scientific-name or ena.metagenome_taxid/ena.metagenome")
    attrs: Dict[str, str] = {}
    attrs.update({str(k): str(v) for k, v in (cc.get("attributes") or {}).items()})
    attrs.update(parse_attributes(args.attribute))
    missing_value = []
    for tag in MANDATORY_ATTRS:
        if not attrs.get(tag):
            attrs[tag] = MISSING_SAMPLE_GROUP
            missing_value.append(tag)
    if missing_value:
        LOG.warning("no value for %s: using %r (INSDC term for a group of samples)",
                    ", ".join(missing_value), MISSING_SAMPLE_GROUP)
    title = args.title or cc.get("title") or f"Co-assembly of {len(accs)} metagenomic samples ({group})"
    description = args.description or cc.get("description") or (
        f"This sample represents the co-assembly of {len(accs)} metagenomic {kind}s"
        + (f" of study {study}" if study else "") + "; the source "
        + f"{kind}s are listed in 'sample derived from'.")
    centre = args.centre_name or cfg_get(cfg, "centre_name")
    sample_set = build_sample_xml(alias, title, description, str(taxon_id), sci, accs, attrs,
                                  args.checklist or cc.get("checklist") or DEFAULT_CHECKLIST, centre)
    submission = build_submission_xml(centre, args.hold_date)
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)
    atomic_write_text(out / "coassembly_sample.xml", _pretty(sample_set))
    atomic_write_text(out / "submission.xml", _pretty(submission))
    atomic_write_text(out / "webin_submission.xml", _pretty(build_webin_document(sample_set, submission)))
    info = {"alias": alias, "group": group, "study": study, "derived_from_kind": kind, "derived_from": accs,
            "runs": runs, "taxon_id": str(taxon_id), "scientific_name": sci, "created": now_iso()}
    atomic_write_text(out / "coassembly_sample.json", json.dumps(info, indent=2) + "\n")
    LOG.info("co-assembly sample %s references %d %s accessions -> %s", alias, len(accs), kind, out)
    return 0


def cmd_submit_sample(args) -> int:
    xml_dir = Path(args.xml_dir)
    info = json.loads((xml_dir / "coassembly_sample.json").read_text())
    sample_set = ET.parse(xml_dir / "coassembly_sample.xml").getroot()
    submission = ET.parse(xml_dir / "submission.xml").getroot()
    mode = "live" if args.live else "test"
    acc_tsv = Path(args.accessions_tsv or xml_dir / "coassembly_sample_accessions.tsv")
    prev = lookup_accession(acc_tsv, info["alias"], mode)
    if prev and not args.force:
        LOG.info("already registered in %s mode: %s (%s); nothing to do", mode, prev["accession"], prev["alias"])
        print(prev["accession"])
        return 0
    alias = info["alias"]
    if not args.live:
        alias = f"{alias}_{args.test_suffix or int(time.time())}"
        sample_set = with_alias(sample_set, alias)
    if args.dry_run:
        print(f"[dry-run] would POST sample {alias} to {ENDPOINTS[(args.api, bool(args.live))]}")
        return 0
    text, rec = submit_sample_set(sample_set, submission, live=args.live, api=args.api)
    receipt = xml_dir / f"receipt_{mode}_{int(time.time())}.xml"
    atomic_write_text(receipt, text)
    acc, biosample = "", ""
    if rec["success"]:
        acc, biosample = rec["samples"].get(alias, next(iter(rec["samples"].values()), ("", "")))
    elif alias in rec["existing"]:
        acc = rec["existing"][alias]
        LOG.warning("sample alias %s already exists in the account as %s; reusing it", alias, acc)
    if not acc:
        LOG.error("ENA did not register the sample (receipt %s):\n  %s", receipt, "\n  ".join(rec["errors"]) or "no message")
        return 1
    append_accession(acc_tsv, {"alias_base": info["alias"], "alias": alias, "accession": acc, "biosample": biosample,
                               "mode": mode, "receipt": str(receipt), "registered_at": now_iso()})
    for msg in rec["infos"]:
        LOG.info("ENA: %s", msg)
    print(acc)
    return 0


def cmd_manifest(args) -> int:
    from mag_coverage import assembly_coverage  # same directory
    group_samples = select_group(load_samples(args.samples), args.group)
    runs = run_accessions(group_samples)
    if args.sample_accession:
        sample = args.sample_accession
    else:
        if not args.accessions_tsv:
            raise EnaError("give --sample-accession or --accessions-tsv")
        rows = [r for r in read_accessions(args.accessions_tsv) if r["mode"] == args.mode]
        if args.alias:
            rows = [r for r in rows if r["alias_base"] == args.alias]
        if len(rows) == 0:
            raise EnaError(f"no co-assembly sample registered in {args.mode} mode in {args.accessions_tsv}")
        sample = rows[-1]["accession"]
    if args.coverage is not None:
        cov, how = float(args.coverage), "given"
    elif args.coverm:
        res = assembly_coverage(args.coverm, args.coverage_method, args.assembly_fasta)
        cov, how = float(res["coverage"]), f"{args.coverage_method} over {res['n_contigs']} contigs, {res['n_samples']} samples"
    else:
        raise EnaError("give --coverm or --coverage")
    out = Path(args.outdir)
    fasta = Path(args.assembly_fasta).resolve()
    if args.no_gzip:
        if not str(fasta).endswith((".gz", ".bz2")):
            raise EnaError("--no-gzip needs an already compressed --assembly-fasta")
        dst = fasta
    else:
        dst = Path(args.fasta_out or out / f"{args.assembly_name}.fa.gz").resolve()
        stats = write_validated_fasta_gz(fasta, dst, force=args.force)
        LOG.info("assembly FASTA %s: %d contigs, %d bp%s", dst, stats["n_seqs"], stats["total_bp"],
                 " (up to date)" if stats["reused"] else "")
    desc = args.description or f"Primary metagenome co-assembly of {len(runs)} runs ({args.program})"
    text = build_manifest(args.study, sample, args.assembly_name, cov, args.program, args.platform, str(dst), runs,
                          desc, args.moleculetype, args.mingaplength)
    mdir = out / ("manifests" if args.mode == "live" else "manifests_test")
    mdir.mkdir(parents=True, exist_ok=True)
    atomic_write_text(mdir / f"{args.assembly_name}.manifest", text)
    LOG.info("manifest %s (COVERAGE %.2f, %s)", mdir / f"{args.assembly_name}.manifest", cov, how)
    print(mdir / f"{args.assembly_name}.manifest")
    return 0


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="count", default=0)
    sub = ap.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("sample-xml", help="write the co-assembly sample XML")
    s.add_argument("--samples", required=True)
    s.add_argument("--group", help="assembly_group (default: the only group)")
    s.add_argument("--params", help="YAML/JSON with an 'ena:' block")
    s.add_argument("--study")
    s.add_argument("--alias")
    s.add_argument("--title")
    s.add_argument("--description")
    s.add_argument("--taxon-id")
    s.add_argument("--scientific-name")
    s.add_argument("--checklist", help=f"default {DEFAULT_CHECKLIST} (ENA default checklist)")
    s.add_argument("--attribute", action="append", default=[], help="TAG=VALUE, repeatable")
    s.add_argument("--centre-name")
    s.add_argument("--hold-date", help="YYYY-MM-DD (adds a HOLD action)")
    s.add_argument("--outdir", required=True)
    s.set_defaults(func=cmd_sample_xml)

    t = sub.add_parser("submit-sample", help="register the co-assembly sample (TEST server unless --live)")
    t.add_argument("--xml-dir", required=True)
    t.add_argument("--live", action="store_true", help="use the PRODUCTION server")
    t.add_argument("--api", choices=("v2", "v1"), default="v2", help="Webin REST v2 (default) or v1 drop-box")
    t.add_argument("--test-suffix", help="suffix for the test alias (default: unix time)")
    t.add_argument("--accessions-tsv")
    t.add_argument("--force", action="store_true", help="register again even if logged")
    t.add_argument("--dry-run", action="store_true")
    t.set_defaults(func=cmd_submit_sample)

    m = sub.add_parser("manifest", help="write the primary metagenome manifest")
    m.add_argument("--samples", required=True)
    m.add_argument("--group")
    m.add_argument("--study", required=True)
    m.add_argument("--mode", choices=("test", "live"), default="test")
    m.add_argument("--sample-accession")
    m.add_argument("--accessions-tsv")
    m.add_argument("--alias", help="alias_base to pick from --accessions-tsv")
    m.add_argument("--coverm", help="CoverM contig table of the group's mapping")
    m.add_argument("--coverage", type=float, help="explicit COVERAGE (overrides --coverm)")
    m.add_argument("--coverage-method", choices=("sum_mean", "length_weighted"), default="sum_mean")
    m.add_argument("--assembly-fasta", required=True)
    m.add_argument("--fasta-out")
    m.add_argument("--no-gzip", action="store_true", help="use --assembly-fasta as is (must be .gz/.bz2)")
    m.add_argument("--assembly-name", required=True)
    m.add_argument("--program", required=True, help='e.g. "MEGAHIT v1.2.9"')
    m.add_argument("--platform", default="ILLUMINA")
    m.add_argument("--description")
    m.add_argument("--moleculetype", default="genomic DNA")
    m.add_argument("--mingaplength", type=int)
    m.add_argument("--outdir", required=True)
    m.add_argument("--force", action="store_true")
    m.set_defaults(func=cmd_manifest)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        return args.func(args)
    except (EnaError, ValueError) as exc:
        LOG.error("%s", exc)
        return 2 if isinstance(exc, CredentialsError) else 1


if __name__ == "__main__":
    sys.exit(main())
