#!/usr/bin/env python3
"""Safe, resumable wrapper around genome_upload and webin_cli_handler.

Defaults are chosen so that nothing reaches the ENA production server by
accident:

* every command talks to the ENA TEST server (wwwdev.ebi.ac.uk; wiped daily)
  unless ``--live`` is given;
* ``--live`` additionally requires (a) ENA_WEBIN and ENA_WEBIN_PASSWORD in the
  environment, (b) a completed TEST run of the same input (same table MD5 /
  same manifest) recorded in the work directory, unless ``--skip-test-gate``,
  and (c) for a primary co-assembly, ``--helpdesk-notified`` (ENA asks to be
  told before co-assembled primary assemblies are submitted);
* credentials are read from the environment only (never from files, never
  printed); subprocess logs are stored with the password redacted;
* ``--dry-run`` prints the commands without contacting ENA or needing
  credentials.

Sub-commands
------------
validate-table  offline check of the genome_uploader TSV with genome_uploader's
                own validator (GenomeUpload.validate_metadata_tsv) plus extra rules
preview-taxa    taxon genome_uploader will register for each NCBI_lineage
                (read-only ENA taxonomy REST service, no credentials)
genomes         genome_upload (sample registration + manifests) -> webin-cli
                validate -> submit, for bins or MAGs
coassembly      co-assembly sample (ena_coassembly.py) -> primary metagenome
                manifest -> webin-cli validate -> submit
status          print the recorded state and accessions

State: ``<workdir>/ena_submit_state.json``; accessions:
``<workdir>/ena_accessions_<test|live>.tsv`` and the co-assembly accession TSV.
Re-running a command skips every step already recorded as done.
"""
from __future__ import annotations

import argparse
import json
import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import (CredentialsError, EnaError, Redactor, ValidationError, atomic_write_text,  # noqa: E402
                        get_webin_credentials, md5_file, now_iso, read_tsv, setup_logging,
                        validate_assembly_name, write_tsv)

LOG = logging.getLogger("bfm.ena_submit")
ERZ_RE = re.compile(r"\b([EDS]RZ\d{6,})\b")
ACCESSION_COLUMNS = ["genome_name", "alias", "sample_accession", "analysis_accession", "validated", "mode", "updated"]
_LINEAGE_RE = re.compile(r"^d__[^;]+(;[pcofgs]__[^;]*){6}$")


# --------------------------------------------------------------------------
# state
# --------------------------------------------------------------------------
class State:
    def __init__(self, workdir: Path):
        self.path = workdir / "ena_submit_state.json"
        self.data = json.loads(self.path.read_text()) if self.path.exists() else {}

    def section(self, kind: str, mode: str) -> Dict:
        return self.data.setdefault(kind, {}).setdefault(mode, {})

    def save(self) -> None:
        atomic_write_text(self.path, json.dumps(self.data, indent=2, sort_keys=True) + "\n")


# --------------------------------------------------------------------------
# subprocess helpers
# --------------------------------------------------------------------------
def _exe(name: str) -> str:
    found = shutil.which(name)
    if found:
        return found
    if Path(name).exists():
        return str(Path(name).resolve())
    raise EnaError(f"{name} not found on PATH (run inside the genome-uploader container/env)")


def run_logged(cmd: List[str], log_path: Path, redact: Redactor, cwd: Path, timeout: int = 7200,
               dry_run: bool = False) -> Tuple[int, str]:
    shown = shlex.join(redact(str(c)) for c in cmd)  # redact each argument before quoting
    if dry_run:
        print(f"[dry-run] (cd {shlex.quote(str(cwd))} && {shown})")
        return 0, ""
    LOG.info("running: %s", shown)
    cwd.mkdir(parents=True, exist_ok=True)
    try:
        res = subprocess.run(cmd, cwd=str(cwd), env=os.environ.copy(), capture_output=True, text=True,
                             timeout=timeout)
        out = redact((res.stdout or "") + (res.stderr or ""))
        code = res.returncode
    except subprocess.TimeoutExpired as exc:
        out = redact(f"timeout after {timeout}s\n{exc.stdout or ''}{exc.stderr or ''}")
        code = 124
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(f"$ {shown}\n# exit {code}  {now_iso()}\n{out}")
    return code, out


def credentials_or_exit(dry_run: bool) -> Redactor:
    if dry_run:
        return Redactor([os.environ.get("ENA_WEBIN_PASSWORD")])
    _user, password = get_webin_credentials(required=True)
    return Redactor([password])


# --------------------------------------------------------------------------
# validate-table / preview-taxa
# --------------------------------------------------------------------------
def extra_table_checks(table: Path) -> List[str]:
    rows = read_tsv(table, required=("genome_name", "genome_path", "NCBI_lineage", "accessions"))
    problems = []
    for r in rows:
        g = r["genome_name"]
        problems += [f"{g}: {p}" for p in validate_assembly_name(g, reserve_test_suffix=True)]
        if not _LINEAGE_RE.match(r["NCBI_lineage"]):
            problems.append(f"{g}: NCBI_lineage must look like d__X;p__..;c__..;o__..;f__..;g__..;s__..")
        for col in ("co-assembly", "rRNA_presence"):
            if col in r and r[col].strip().lower() not in ("true", "false", "yes", "no", "y", "n"):
                problems.append(f"{g}: {col} must be True or False (got {r[col]!r})")
        p = Path(r["genome_path"])
        if not p.is_absolute():
            problems.append(f"{g}: genome_path should be absolute ({p})")
        if not str(p).endswith((".gz", ".bz2")):
            problems.append(f"{g}: genome_path must be gzip/bzip2 compressed")
        if not p.exists():
            problems.append(f"{g}: genome_path does not exist ({p})")
    return problems


def genome_uploader_validation(table: Path, genome_type: str) -> List[str]:
    """Run genome_uploader's own TSV validation offline; returns problems ([] if valid)."""
    try:
        from genomeuploader.genome_upload import GenomeUpload  # type: ignore
    except Exception as exc:
        return [f"genome_uploader not importable ({exc}); offline validation skipped"]
    with tempfile.TemporaryDirectory() as td:
        args = {"bins": genome_type == "bins", "mags": genome_type == "mags", "live": False, "private": False,
                "tpa": False, "centre_name": "offline-validation", "force": False, "out": td,
                "upload_study": "PRJEB00000", "genome_info": str(table), "test_suffix": "0000000000",
                "chromosome_info": None}
        try:
            gu = GenomeUpload(args)
            gu.validate_metadata_tsv()
        except (ValueError, KeyError, TypeError) as exc:
            return [f"genome_uploader: {exc}"]
    return []


def cmd_validate_table(args) -> int:
    table = Path(args.genome_info)
    problems = extra_table_checks(table) + genome_uploader_validation(table, args.genome_type)
    hard = [p for p in problems if "validation skipped" not in p]
    for p in problems:
        (LOG.error if p in hard else LOG.warning)("%s", p)
    if hard:
        return 1
    print(f"{table}: valid for genome_uploader ({args.genome_type})")
    return 0


def cmd_preview_taxa(args) -> int:
    try:
        import genomeuploader.taxon_finder as tfmod  # type: ignore
    except Exception as exc:
        raise EnaError(f"genome_uploader not importable: {exc}")
    import types
    real_requests = tfmod.requests
    url_cache: Dict[str, object] = {}

    def cached_get(url, *a, **k):  # TaxonFinder re-queries the same names for many genomes
        if url not in url_cache:
            url_cache[url] = real_requests.get(url, *a, **k)
            time.sleep(args.pause)
        return url_cache[url]

    tfmod.requests = types.SimpleNamespace(get=cached_get, exceptions=real_requests.exceptions)
    try:
        rows = read_tsv(args.genome_info, required=("genome_name", "NCBI_lineage"))
        cache: Dict[str, Tuple[str, str]] = {}
        out = []
        for r in rows:
            lin = r["NCBI_lineage"]
            if lin not in cache:
                try:
                    tf = tfmod.TaxonFinder(lin)
                    cache[lin] = (str(tf.taxid), str(tf.scientific_name))
                except Exception as exc:  # network or format problems are reported per lineage
                    cache[lin] = ("", f"ERROR: {exc}")
            out.append({"genome_name": r["genome_name"], "NCBI_lineage": lin, "taxid": cache[lin][0],
                        "scientific_name": cache[lin][1]})
    finally:
        tfmod.requests = real_requests
    write_tsv(args.out, out, ["genome_name", "NCBI_lineage", "taxid", "scientific_name"])
    bad = [o for o in out if not o["taxid"]]
    LOG.info("%d genomes, %d distinct lineages, %d ENA queries, %d without a submittable taxon -> %s",
             len(out), len(cache), len(url_cache), len(bad), args.out)
    return 1 if bad else 0


# --------------------------------------------------------------------------
# webin-cli through webin_cli_handler
# --------------------------------------------------------------------------
class Webin:
    def __init__(self, exe: str, workdir: Path, mode: str, redact: Redactor, heap_max_gb: Optional[int],
                 heap_init_gb: Optional[int], dry_run: bool, timeout: int = 7200):
        self.exe = exe
        self.workdir = workdir
        self.mode = mode
        self.redact = redact
        self.dry_run = dry_run
        self.timeout = timeout
        self.heap = []
        if heap_max_gb and self._supports("--java-heap-size-max"):
            self.heap = ["--java-heap-size-initial", str(heap_init_gb or 1), "--java-heap-size-max", str(heap_max_gb)]

    def _supports(self, flag: str) -> bool:
        if self.dry_run and not shutil.which(self.exe) and not Path(self.exe).exists():
            return True
        try:
            res = subprocess.run([self.exe, "--help"], capture_output=True, text=True, timeout=120)
            return flag in (res.stdout + res.stderr)
        except Exception:
            return False

    def run(self, manifest: Path, alias: str, action: str) -> Dict[str, object]:
        safe = re.sub(r"[^A-Za-z0-9_.#-]", "_", alias)
        rundir = self.workdir / f"webin_{self.mode}" / safe
        rundir.mkdir(parents=True, exist_ok=True)
        local = rundir / f"{safe}.manifest"
        text = Path(manifest).read_text()
        fasta = re.search(r"^FASTA\s+(\S+)", text, re.M)
        if not fasta or not Path(fasta.group(1)).is_absolute() or not Path(fasta.group(1)).exists():
            raise EnaError(f"{manifest}: FASTA path must be absolute and exist")
        local.write_text(text)
        for old in [rundir / "webin-cli.report", rundir / f"{local.name}.report"]:
            if old.exists():
                old.unlink()
        cmd = [self.exe, "-m", str(local), "-c", "genome", "--mode", action]
        if self.mode == "test":
            cmd.append("--test")
        cmd += self.heap
        code, out = run_logged(cmd, self.workdir / "logs" / self.mode / f"webin_{action}_{safe}.log",
                               self.redact, rundir, self.timeout, self.dry_run)
        if self.dry_run:
            return {"ok": True, "accession": "", "dry_run": True}
        report = "\n".join(p.read_text(errors="replace") for p in rundir.glob("*.report"))
        receipts = sorted(rundir.glob(f"genome/*/{action}/receipt.xml"))
        receipt_text = receipts[-1].read_text(errors="replace") if receipts else ""
        acc = ""
        m = re.search(r'ANALYSIS[^>]*accession="([EDS]RZ\d+)"', receipt_text) or ERZ_RE.search(report) \
            or ERZ_RE.search(out)
        if m:
            acc = m.group(1)
        ok = code == 0
        if action == "validate":
            ok = ok and ("validated successfully" in (report + out) or "validation succeeded" in out.lower())
        return {"ok": ok, "accession": acc, "exit": code, "log": str(self.workdir / "logs" / self.mode)}


# --------------------------------------------------------------------------
# genomes
# --------------------------------------------------------------------------
def write_accessions(path: Path, sec: Dict, mode: str) -> None:
    rows = []
    for g, info in sorted(sec.get("genomes", {}).items()):
        rows.append({"genome_name": g, "alias": info.get("alias", ""), "sample_accession": info.get("sample", ""),
                     "analysis_accession": info.get("analysis", ""),
                     "validated": "True" if info.get("validated_md5") else "False",
                     "mode": mode, "updated": info.get("updated", "")})
    write_tsv(path, rows, ACCESSION_COLUMNS)


def cmd_genomes(args) -> int:
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    mode = "live" if args.live else "test"
    table = Path(args.genome_info).resolve()
    table_md5 = md5_file(table)
    redact = credentials_or_exit(args.dry_run)
    state = State(workdir)
    sec = state.section("genomes", mode)
    if sec.get("table_md5") and sec["table_md5"] != table_md5 and not args.restart:
        raise EnaError(f"{mode} state in {workdir} belongs to a different input table; use --restart "
                       "(or a new --workdir) to start over")
    if args.restart:
        sec.clear()
    if args.live and not args.skip_test_gate:
        t = state.section("genomes", "test")
        if not (t.get("table_md5") == table_md5 and t.get("complete")):
            raise EnaError("refusing --live: no completed TEST submission of this exact table in this workdir "
                           "(run without --live first, or pass --skip-test-gate)")
    problems = extra_table_checks(table) + genome_uploader_validation(table, args.genome_type)
    hard = [p for p in problems if "validation skipped" not in p]
    if hard:
        raise ValidationError(hard, header="genome table is not valid")
    genomes = [r["genome_name"] for r in read_tsv(table, required=("genome_name",))]
    sec["table_md5"] = table_md5
    sec.setdefault("genomes", {})
    if mode == "test":
        sec.setdefault("test_suffix", args.test_suffix or str(int(time.time())))
    if args.live:
        LOG.warning("LIVE mode: samples and assemblies will be registered on the ENA PRODUCTION server")

    # 1. sample registration + manifests (genome_upload)
    upload_root = workdir / f"genome_upload_{mode}"
    upload_dir = upload_root / ("MAG_upload" if args.genome_type == "mags" else "bin_upload")
    reg_name = ("registered_MAGs" if args.genome_type == "mags" else "registered_bins") + ("" if args.live else "_test") + ".tsv"
    manifest_dir = upload_dir / ("manifests" if args.live else "manifests_test")
    suffix = "" if args.live else "_" + sec["test_suffix"]
    registered_all = all(sec["genomes"].get(g, {}).get("sample") for g in genomes) and manifest_dir.exists()
    if not registered_all or args.force_register:
        cmd = [_exe(args.genome_upload_exe), "-u", args.study, "--genome_info", str(table),
               "--mags" if args.genome_type == "mags" else "--bins", "--out", str(upload_root),
               "--centre_name", args.centre_name]
        if args.private:
            cmd.append("--private")
        if args.live:
            cmd.append("--live")
        else:
            cmd += ["--test-suffix", sec["test_suffix"]]
        if args.force_register:
            cmd.append("--force")
        code, out = run_logged(cmd, workdir / "logs" / mode / "genome_upload.log", redact, upload_root,
                               dry_run=args.dry_run)
        if args.dry_run:
            pass
        elif code != 0:
            raise EnaError(f"genome_upload failed (exit {code}); see {workdir / 'logs' / mode / 'genome_upload.log'}")
        else:
            reg_file = upload_dir / reg_name
            reg = {}
            if reg_file.exists():
                for line in reg_file.read_text().splitlines():
                    if "\t" in line:
                        a, acc = line.split("\t", 1)
                        reg[a.strip()] = acc.strip()
            for g in genomes:
                alias = g + suffix
                if alias not in reg:
                    raise EnaError(f"{g}: alias {alias} missing from {reg_file}")
                sec["genomes"].setdefault(g, {}).update({"alias": alias, "sample": reg[alias], "updated": now_iso()})
            state.save()
    if args.dry_run:
        for g in genomes:
            print(f"[dry-run] would validate and submit {manifest_dir / (g + suffix + '.manifest')}")
        return 0
    write_accessions(workdir / f"ena_accessions_{mode}.tsv", sec, mode)
    if args.stop_after == "register":
        return 0

    # 2./3. validate and submit every manifest with webin_cli_handler
    webin = Webin(_exe(args.webin_cli_handler_exe), workdir, mode, redact, args.java_heap_max_gb,
                  args.java_heap_init_gb, args.dry_run)
    failed = []
    for g in genomes:
        info = sec["genomes"][g]
        mf = manifest_dir / f"{info['alias']}.manifest"
        if not mf.exists():
            failed.append(f"{g}: manifest {mf} missing")
            continue
        mmd5 = md5_file(mf)
        if info.get("analysis"):
            continue
        if info.get("validated_md5") != mmd5:
            res = webin.run(mf, info["alias"], "validate")
            if not res["ok"]:
                failed.append(f"{g}: validation failed (logs in {res['log']})")
                continue
            info.update({"validated_md5": mmd5, "updated": now_iso()})
            state.save()
        if args.stop_after == "validate":
            continue
        res = webin.run(mf, info["alias"], "submit")
        if res["ok"] and (res["accession"] or mode == "test"):
            # on the TEST server a re-submitted object gives a minimal report without accession
            info.update({"analysis": res["accession"] or "TEST-already-submitted", "updated": now_iso()})
            state.save()
            write_accessions(workdir / f"ena_accessions_{mode}.tsv", sec, mode)
            time.sleep(args.delay)
        else:
            failed.append(f"{g}: submission failed or no ERZ accession found (logs in {res['log']})")
    done = [g for g in genomes if sec["genomes"].get(g, {}).get("analysis")]
    sec["complete"] = (len(done) == len(genomes))
    state.save()
    write_accessions(workdir / f"ena_accessions_{mode}.tsv", sec, mode)
    LOG.info("%s mode: %d/%d genomes submitted", mode, len(done), len(genomes))
    if failed:
        LOG.error("problems:\n  %s", "\n  ".join(failed))
        return 1
    return 0


# --------------------------------------------------------------------------
# primary co-assembly
# --------------------------------------------------------------------------
def cmd_coassembly(args) -> int:
    import ena_coassembly as ec
    workdir = Path(args.workdir).resolve()
    workdir.mkdir(parents=True, exist_ok=True)
    mode = "live" if args.live else "test"
    if args.live and not args.helpdesk_notified:
        raise EnaError("refusing --live: ENA asks submitters to inform the helpdesk before submitting a primary "
                       "assembly co-assembled from raw reads; pass --helpdesk-notified once that is done")
    redact = credentials_or_exit(args.dry_run)
    state = State(workdir)
    sec = state.section("coassembly", mode)
    xml_dir = Path(args.xml_dir).resolve()
    info = json.loads((xml_dir / "coassembly_sample.json").read_text())
    key = info["alias"]
    if args.live and not args.skip_test_gate and not state.section("coassembly", "test").get(key, {}).get("analysis"):
        raise EnaError("refusing --live: the co-assembly has no completed TEST submission in this workdir "
                       "(run without --live first, or pass --skip-test-gate)")
    entry = sec.setdefault(key, {})
    acc_tsv = xml_dir / "coassembly_sample_accessions.tsv"
    # 1. co-assembly sample
    if not entry.get("sample"):
        prev = ec.lookup_accession(acc_tsv, key, mode)
        if prev:
            entry["sample"] = prev["accession"]
        elif args.dry_run:
            print(f"[dry-run] would register co-assembly sample {key} on the {mode} server")
        else:
            rc = ec.main(["submit-sample", "--xml-dir", str(xml_dir)] + (["--live"] if args.live else []))
            if rc != 0:
                raise EnaError("co-assembly sample registration failed")
            entry["sample"] = ec.lookup_accession(acc_tsv, key, mode)["accession"]
        state.save()
    # 2. manifest
    mdir = Path(args.outdir).resolve()
    manifest = mdir / ("manifests" if args.live else "manifests_test") / f"{args.assembly_name}.manifest"
    if args.dry_run and not entry.get("sample"):
        print(f"[dry-run] would write {manifest}, then webin-cli validate + submit")
        return 0
    margs = ["manifest", "--samples", args.samples, "--study", args.study, "--mode", mode,
             "--sample-accession", entry["sample"], "--assembly-fasta", args.assembly_fasta,
             "--assembly-name", args.assembly_name, "--program", args.program, "--platform", args.platform,
             "--outdir", str(mdir)]
    if args.group:
        margs += ["--group", args.group]
    if args.coverage is not None:
        margs += ["--coverage", str(args.coverage)]
    else:
        margs += ["--coverm", args.coverm, "--coverage-method", args.coverage_method]
    if args.description:
        margs += ["--description", args.description]
    if ec.main(margs) != 0:
        raise EnaError("manifest creation failed")
    # 3./4. validate + submit
    webin = Webin(_exe(args.webin_cli_handler_exe), workdir, mode, redact, args.java_heap_max_gb,
                  args.java_heap_init_gb, args.dry_run)
    mmd5 = md5_file(manifest)
    if entry.get("analysis"):
        LOG.info("primary assembly already submitted in %s mode: %s", mode, entry["analysis"])
        return 0
    if entry.get("validated_md5") != mmd5:
        res = webin.run(manifest, args.assembly_name, "validate")
        if not res["ok"]:
            raise EnaError(f"webin-cli validation failed; logs in {res['log']}")
        entry["validated_md5"] = mmd5
        state.save()
    if args.stop_after == "validate":
        return 0
    res = webin.run(manifest, args.assembly_name, "submit")
    if not (res["ok"] and (res["accession"] or mode == "test")):
        raise EnaError(f"submission failed or no ERZ accession found; logs in {res['log']}")
    entry.update({"analysis": res["accession"] or "TEST-already-submitted", "manifest": str(manifest),
                  "updated": now_iso()})
    state.save()
    print(res["accession"])
    return 0


def cmd_status(args) -> int:
    state = State(Path(args.workdir))
    print(json.dumps(state.data, indent=2, sort_keys=True))
    return 0


# --------------------------------------------------------------------------
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-v", "--verbose", action="count", default=0)
    sub = ap.add_subparsers(dest="cmd", required=True)

    def common(p):
        p.add_argument("--workdir", required=True, help="state, logs and webin-cli outputs")
        p.add_argument("--live", action="store_true", help="use the ENA PRODUCTION server (default: TEST)")
        p.add_argument("--skip-test-gate", action="store_true", help="allow --live without a recorded TEST run")
        p.add_argument("--dry-run", action="store_true", help="print commands only")
        p.add_argument("--stop-after", choices=("register", "validate", "submit"), default="submit")
        p.add_argument("--webin-cli-handler-exe", default="webin_cli_handler")
        p.add_argument("--java-heap-max-gb", type=int, default=4)
        p.add_argument("--java-heap-init-gb", type=int, default=1)

    v = sub.add_parser("validate-table", help="offline validation of the genome_uploader TSV")
    v.add_argument("--genome-info", required=True)
    v.add_argument("--genome-type", choices=("bins", "mags"), required=True)
    v.set_defaults(func=cmd_validate_table)

    t = sub.add_parser("preview-taxa", help="taxa genome_uploader will register (read-only ENA taxonomy API)")
    t.add_argument("--genome-info", required=True)
    t.add_argument("--out", required=True)
    t.add_argument("--pause", type=float, default=0.2, help="seconds between distinct lineages")
    t.set_defaults(func=cmd_preview_taxa)

    g = sub.add_parser("genomes", help="register and submit bins/MAGs")
    common(g)
    g.add_argument("--genome-info", required=True, help="ena_genome_table.py output")
    g.add_argument("--genome-type", choices=("bins", "mags"), required=True)
    g.add_argument("--study", required=True, help="PRJEB/ERP accession the genomes are submitted to")
    g.add_argument("--centre-name", required=True)
    g.add_argument("--private", dest="private", action="store_true", default=True,
                   help="raw data still private on ENA (default; genome_uploader --private)")
    g.add_argument("--public", dest="private", action="store_false", help="raw data already public")
    g.add_argument("--test-suffix", help="fixed alias suffix for TEST registrations (default: unix time, kept in state)")
    g.add_argument("--force-register", action="store_true", help="re-run genome_upload with --force")
    g.add_argument("--restart", action="store_true", help="forget the recorded state for this mode")
    g.add_argument("--delay", type=float, default=2.0, help="seconds between webin-cli submissions")
    g.add_argument("--genome-upload-exe", default="genome_upload")
    g.set_defaults(func=cmd_genomes)

    c = sub.add_parser("coassembly", help="register the co-assembly sample and submit the primary assembly")
    common(c)
    c.add_argument("--xml-dir", required=True, help="ena_coassembly.py sample-xml output directory")
    c.add_argument("--samples", required=True)
    c.add_argument("--group")
    c.add_argument("--study", required=True)
    c.add_argument("--assembly-fasta", required=True)
    c.add_argument("--assembly-name", required=True)
    c.add_argument("--program", required=True)
    c.add_argument("--platform", default="ILLUMINA")
    c.add_argument("--coverm")
    c.add_argument("--coverage", type=float)
    c.add_argument("--coverage-method", choices=("sum_mean", "length_weighted"), default="sum_mean")
    c.add_argument("--description")
    c.add_argument("--outdir", required=True, help="where manifests and the gzipped assembly go")
    c.add_argument("--helpdesk-notified", action="store_true",
                   help="confirm that the ENA helpdesk was told about this co-assembly (required for --live)")
    c.set_defaults(func=cmd_coassembly)

    s = sub.add_parser("status")
    s.add_argument("--workdir", required=True)
    s.set_defaults(func=cmd_status)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        return args.func(args)
    except CredentialsError as exc:
        LOG.error("%s", exc)
        return 2
    except (EnaError, ValueError) as exc:
        LOG.error("%s", exc)
        return 1


if __name__ == "__main__":
    sys.exit(main())
