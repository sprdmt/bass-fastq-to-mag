#!/usr/bin/env python3
"""Build the genome_uploader input table (and gzipped FASTAs) for the final MAG set.

The table has exactly the 17 columns genome_uploader 2.5.1 reads
(``MAG_MANDATORY_FIELDS + BIN_MANDATORY_FIELDS``), in the order of its README:

  genome_name genome_path accessions assembly_software binning_software
  binning_parameters stats_generation_software completeness contamination
  genome_coverage metagenome co-assembly broad_environment local_environment
  environmental_medium rRNA_presence NCBI_lineage   [chromosome_level]

``chromosome_level`` (genome_uploader >= 3.0, optional) is only written with
``--include-chromosome-level``.

Sources
-------
* --mags-table  workflow mags.tsv: mag_id, group, completeness, contamination (CheckM2)
* --mags-dir    <mag_id>.fna (plain or gzipped); copied to --fasta-outdir/<genome_name>.fa.gz
                after ENA sequence checks (>= 20 bp, no terminal N, IUPAC only, unique names)
* --rna         mimag_rna.py output (rRNA_presence, mimag_quality)
* --coverage    mag_coverage.py output (genome_coverage)
* --lineage     gtdb_to_ncbi.py output (ncbi_lineage)
* --samples     samples.tsv: run_accession of every sample of the MAG's assembly
                group -> accessions (comma-separated); co-assembly = more than one run
* --params      YAML/JSON with an ``ena:`` block (metagenome, broad_environment,
                local_environment, environmental_medium, assembly_software,
                binning_software, binning_parameters, stats_generation_software,
                optional per-group overrides under ``groups: {<group>: {...}}``);
                every key can also be given on the command line.

Validation (all problems are reported together; nothing is written on failure)
--------------------------------------------------------------------------
mandatory values present; completeness/contamination within 0-100; coverage > 0;
MAG criteria for --genome-type mags (completeness >= --min-completeness,
contamination < --max-contamination); genome names unique, ENA assembly-name
rule and <= 39 characters (genome_uploader appends an 11-character test
suffix); run accessions [EDS]RRnnnnnn; NCBI_lineage prefixed with a domain;
metagenome in genome_uploader's list (when genomeuploader is importable);
genome_uploader derives the ENA "assembly quality" text itself (HQ text when
completeness >= 90, contamination <= 5 and rRNA_presence). MIMAG high-quality
needs > 90 and < 5, so genomes exactly at 90.00 % / 5.00 % with rRNA_presence
True would receive the HQ text without being MIMAG high-quality: they stop the
build unless --allow-assembly-quality-conflict is given.

Usage
-----
  ena_genome_table.py --mags-table results/mags/mags.tsv --mags-dir results/mags/fasta \
      --rna results/mags/mimag_rna.tsv --coverage results/mags/coverage.tsv \
      --lineage results/ena/ncbi_lineage.tsv --samples config/samples.tsv \
      --params results/ena/ena_params.yaml --genome-type bins \
      --fasta-outdir results/ena/fasta --out results/ena/genome_uploader_input.tsv
"""
from __future__ import annotations

import argparse
import logging
import math
import os
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ena_common import (EnaError, ValidationError, cfg_get, fmt_float, genome_id_from_path,  # noqa: E402
                        list_fasta_files, load_config, load_samples, parse_bool, read_tsv,
                        run_accessions, select_group, setup_logging, validate_assembly_name,
                        write_tsv, write_validated_fasta_gz)

LOG = logging.getLogger("bfm.ena_genome_table")

GENOME_UPLOADER_COLUMNS = [
    "genome_name", "genome_path", "accessions", "assembly_software", "binning_software",
    "binning_parameters", "stats_generation_software", "completeness", "contamination",
    "genome_coverage", "metagenome", "co-assembly", "broad_environment", "local_environment",
    "environmental_medium", "rRNA_presence", "NCBI_lineage",
]
PARAM_KEYS = ["metagenome", "broad_environment", "local_environment", "environmental_medium",
              "assembly_software", "binning_software", "binning_parameters", "stats_generation_software"]
PROVENANCE_COLUMNS = ["genome_name", "mag_id", "group", "source_fasta", "genome_path", "n_seqs", "total_bp",
                      "md5", "completeness", "contamination", "genome_coverage", "mimag_quality",
                      "genome_uploader_quality", "rRNA_presence", "n_runs", "co-assembly",
                      "lineage_method", "lineage_note"]
_LINEAGE_RE = re.compile(r"^d__[^;]+(;[pcofgs]__[^;]*){6}$")


def genome_uploader_metagenomes() -> Optional[set]:
    try:
        from genomeuploader.constants import METAGENOMES  # type: ignore
        return set(METAGENOMES)
    except Exception:  # pragma: no cover - depends on environment
        return None


def index_by(rows: List[Dict[str, str]], key: str, what: str) -> Dict[str, Dict[str, str]]:
    out = {}
    for r in rows:
        k = r.get(key, "")
        if k in out:
            raise EnaError(f"{what}: duplicated {key} {k}")
        out[k] = r
    return out


def resolve_params(cfg: Dict, cli: Dict[str, Optional[str]], group: str) -> Dict[str, str]:
    vals = {}
    for k in PARAM_KEYS:
        v = cli.get(k)
        if v in (None, ""):
            v = cfg_get(cfg, f"groups.{group}.{k}") if group else None
        if v in (None, ""):
            v = cfg_get(cfg, k)
        vals[k] = "" if v is None else str(v).strip()
    return vals


def build_table(mags: List[Dict[str, str]], fasta_index: Dict[str, Path], rna: Dict[str, Dict[str, str]],
                coverage: Dict[str, Dict[str, str]], lineage: Dict[str, Dict[str, str]],
                samples: List[Dict[str, str]], cfg: Dict, cli: Dict[str, Optional[str]], genome_type: str,
                fasta_outdir, min_completeness: float = 50.0, max_contamination: float = 10.0,
                name_prefix: str = "", include_chromosome_level: bool = False, write_fastas: bool = True,
                force: bool = False, allow_quality_conflict: bool = False):
    problems: List[str] = []
    rows, prov = [], []
    names_seen = set()
    known_mg = genome_uploader_metagenomes()
    runs_by_group: Dict[str, List[str]] = {}
    planned = []
    for m in mags:
        mag = m["mag_id"]
        group = m.get("group", "") or ""
        gname = f"{name_prefix}{mag}"
        p = []
        if gname in names_seen:
            p.append(f"{gname}: duplicated genome name")
        names_seen.add(gname)
        p += [f"{gname}: {x}" for x in validate_assembly_name(gname, reserve_test_suffix=True)]
        src = fasta_index.get(mag)
        if src is None:
            p.append(f"{mag}: FASTA not found")
        # quality
        try:
            comp = float(m["completeness"])
            cont = float(m["contamination"])
            if not (0 <= comp <= 100 and 0 <= cont <= 100):
                p.append(f"{mag}: completeness/contamination outside 0-100 ({comp}, {cont})")
            if genome_type == "mags" and (comp < min_completeness or cont >= max_contamination):
                p.append(f"{mag}: fails MAG criteria (completeness {comp} >= {min_completeness}, "
                         f"contamination {cont} < {max_contamination})")
        except (KeyError, ValueError):
            p.append(f"{mag}: completeness/contamination missing or not numeric")
            comp = cont = math.nan
        # rna
        r = rna.get(mag)
        if r is None:
            p.append(f"{mag}: not in --rna table")
            presence = None
        else:
            try:
                presence = parse_bool(r.get("rRNA_presence", ""), "rRNA_presence")
            except ValueError as exc:
                p.append(f"{mag}: {exc}")
                presence = None
        if presence and not math.isnan(comp):
            c2, k2 = round(comp, 2), round(cont, 2)
            if (c2 >= 90 and k2 <= 5) and not (c2 > 90 and k2 < 5) and not allow_quality_conflict:
                p.append(f"{mag}: genome_uploader would write the MIMAG high-quality 'assembly quality' text "
                         f"(completeness {c2}, contamination {k2}) although MIMAG HQ needs > 90 and < 5; "
                         "use --allow-assembly-quality-conflict to accept")
        # coverage
        c = coverage.get(mag)
        cov = math.nan
        if c is None:
            p.append(f"{mag}: not in --coverage table")
        else:
            try:
                cov = float(c.get("genome_coverage", ""))
                if not (cov > 0 and math.isfinite(cov)):
                    p.append(f"{mag}: genome_coverage must be > 0 (got {cov})")
            except ValueError:
                p.append(f"{mag}: genome_coverage missing")
        # lineage
        ln = lineage.get(mag)
        lin = (ln or {}).get("ncbi_lineage", "")
        if ln is None:
            p.append(f"{mag}: not in --lineage table")
        elif not _LINEAGE_RE.match(lin):
            p.append(f"{mag}: NCBI_lineage {lin!r} is not 'd__X;p__;c__;o__;f__;g__;s__' formatted")
        # runs
        if group not in runs_by_group:
            try:
                runs_by_group[group] = run_accessions(select_group(samples, group or None))
            except EnaError as exc:
                runs_by_group[group] = []
                problems.append(f"group {group or '(all)'}: {exc}")
        runs = runs_by_group[group]
        # params
        par = resolve_params(cfg, cli, group)
        for k, v in par.items():
            if not v:
                p.append(f"{mag}: parameter {k} is empty (set ena.{k} or --{k.replace('_', '-')})")
        if known_mg is not None and par["metagenome"] and par["metagenome"] not in known_mg:
            p.append(f"{mag}: metagenome {par['metagenome']!r} is not in genome_uploader's accepted list")
        problems += p
        planned.append((m, gname, src, comp, cont, presence, cov, lin, runs, par, r, ln))
    if problems:
        raise ValidationError(problems, header=f"{len(problems)} problem(s) in the ENA genome table inputs")
    fasta_outdir = Path(fasta_outdir)
    for (m, gname, src, comp, cont, presence, cov, lin, runs, par, r, ln) in planned:
        dst = (fasta_outdir / f"{gname}.fa.gz").resolve()
        stats = {"n_seqs": "", "total_bp": "", "md5": ""}
        if write_fastas:
            stats = write_validated_fasta_gz(src, dst, force=force)
        coassembly = len(runs) > 1
        row = {"genome_name": gname, "genome_path": str(dst), "accessions": ",".join(runs),
               "completeness": fmt_float(comp), "contamination": fmt_float(cont),
               "genome_coverage": fmt_float(cov), "co-assembly": "True" if coassembly else "False",
               "rRNA_presence": "True" if presence else "False", "NCBI_lineage": lin}
        row.update({k: par[k] for k in PARAM_KEYS})
        if include_chromosome_level:
            row["chromosome_level"] = "False"
        rows.append(row)
        prov.append({"genome_name": gname, "mag_id": m["mag_id"], "group": m.get("group", ""),
                     "source_fasta": str(src), "genome_path": str(dst), "n_seqs": stats["n_seqs"],
                     "total_bp": stats["total_bp"], "md5": stats["md5"], "completeness": row["completeness"],
                     "contamination": row["contamination"], "genome_coverage": row["genome_coverage"],
                     "mimag_quality": (r or {}).get("mimag_quality", ""),
                     "genome_uploader_quality": (r or {}).get("genome_uploader_quality", ""),
                     "rRNA_presence": row["rRNA_presence"], "n_runs": len(runs), "co-assembly": row["co-assembly"],
                     "lineage_method": (ln or {}).get("method", ""), "lineage_note": (ln or {}).get("note", "")})
    return rows, prov


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--mags-table", required=True)
    ap.add_argument("--mags-dir", required=True)
    ap.add_argument("--rna", required=True)
    ap.add_argument("--coverage", required=True)
    ap.add_argument("--lineage", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--params", help="YAML/JSON with an 'ena:' block (config.yaml works)")
    ap.add_argument("--genome-type", choices=("bins", "mags"), default=None,
                    help="ENA tier: bins (binned metagenome, ERC000050) or mags (MAG, ERC000047); "
                         "default from ena.genome_type")
    ap.add_argument("--fasta-outdir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--provenance-out", help="default: <out>.provenance.tsv")
    ap.add_argument("--mag-ids", help="file with mag_ids to include (default: all rows of --mags-table)")
    ap.add_argument("--name-prefix", default="", help="prefix added to mag_id to form genome_name")
    ap.add_argument("--min-completeness", type=float, default=None)
    ap.add_argument("--max-contamination", type=float, default=None)
    ap.add_argument("--include-chromosome-level", action="store_true")
    ap.add_argument("--no-fasta", action="store_true", help="validate and write the table without gzipping FASTAs")
    ap.add_argument("--force", action="store_true", help="rewrite gzipped FASTAs even if up to date")
    ap.add_argument("--allow-assembly-quality-conflict", action="store_true",
                    help="accept genomes at exactly 90.00 %% / 5.00 %% that genome_uploader would label HQ")
    for k in PARAM_KEYS:
        ap.add_argument("--" + k.replace("_", "-"), dest=k, default=None)
    ap.add_argument("-v", "--verbose", action="count", default=0)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    setup_logging(args.verbose)
    try:
        cfg = load_config(args.params) if args.params else {}
        genome_type = args.genome_type or cfg_get(cfg, "genome_type")
        if genome_type not in ("bins", "mags"):
            raise EnaError("set --genome-type (or ena.genome_type) to 'bins' or 'mags'")
        mags = read_tsv(args.mags_table, required=("mag_id", "completeness", "contamination"))
        if args.mag_ids:
            keep = {l.strip() for l in open(args.mag_ids) if l.strip()}
            mags = [m for m in mags if m["mag_id"] in keep]
            missing = keep - {m["mag_id"] for m in mags}
            if missing:
                raise EnaError("--mag-ids not in mags table: " + ", ".join(sorted(missing)))
        fasta_index = {genome_id_from_path(p): p for p in list_fasta_files(args.mags_dir)}
        rows, prov = build_table(
            mags, fasta_index,
            index_by(read_tsv(args.rna, required=("mag_id", "rRNA_presence")), "mag_id", "rna"),
            index_by(read_tsv(args.coverage, required=("mag_id", "genome_coverage")), "mag_id", "coverage"),
            index_by(read_tsv(args.lineage, required=("mag_id", "ncbi_lineage")), "mag_id", "lineage"),
            load_samples(args.samples), cfg, {k: getattr(args, k) for k in PARAM_KEYS}, genome_type,
            args.fasta_outdir,
            min_completeness=args.min_completeness if args.min_completeness is not None
            else float(cfg_get(cfg, "min_completeness", 50)),
            max_contamination=args.max_contamination if args.max_contamination is not None
            else float(cfg_get(cfg, "max_contamination", 10)),
            name_prefix=args.name_prefix, include_chromosome_level=args.include_chromosome_level,
            write_fastas=not args.no_fasta, force=args.force,
            allow_quality_conflict=args.allow_assembly_quality_conflict)
    except (EnaError, ValueError) as exc:
        LOG.error("%s", exc)
        return 1
    cols = GENOME_UPLOADER_COLUMNS + (["chromosome_level"] if args.include_chromosome_level else [])
    write_tsv(args.out, rows, cols)
    write_tsv(args.provenance_out or str(Path(args.out).with_suffix("")) + ".provenance.tsv", prov, PROVENANCE_COLUMNS)
    LOG.info("wrote %d genomes (%s) to %s", len(rows), genome_type, args.out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
