"""Tests for the workflow helper scripts (run with: pytest tests/)."""

import csv
import gzip
import importlib.util
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[1] / "workflow" / "scripts"


def load(name):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write(path, text):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return path


# ---------------------------------------------------------------------------
def test_bins_to_table_detects_fasta_by_content(tmp_path):
    mod = load("bins_to_table")
    bins = tmp_path / "bins"
    write(bins / "bin.1.fa", ">k125_1 flag=0 multi=3.0 len=1500\nACGT\n>k125_2\nACGT\n")
    with gzip.open(bins / "bin.2.fa.gz", "wt") as fh:
        fh.write(">k125_3 desc\nACGT\n")
    write(bins / "7", ">k125_4\nACGT\n")            # GenomeFace-style name without extension
    write(bins / "bins.tsv", "0\tk125_4\n")         # excluded
    write(bins / "notes.txt", "not a fasta\n")      # ignored
    out = tmp_path / "table.tsv"
    mod.main(["--bins-dir", str(bins), "--out", str(out), "--prefix", "mb_", "--exclude", "bins.tsv"])
    rows = [line.split("\t") for line in out.read_text().splitlines()]
    assert rows == [["k125_4", "mb_7"], ["k125_1", "mb_bin.1"], ["k125_2", "mb_bin.1"], ["k125_3", "mb_bin.2"]]


def test_fix_vamb_taxonomy_truncates_and_disambiguates(tmp_path):
    mod = load("fix_vamb_taxonomy")
    contigs = write(tmp_path / "c.fa", ">c1\nA\n>c2\nA\n>c3\nA\n>c4\nA\n")
    tax = write(
        tmp_path / "tax.tsv",
        "contigs\tpredictions\n"
        "c1\tBacteria;Pseudomonadota;Gammaproteobacteria;unknown;X\n"
        "c2\tBacteria;Bacillota;Gammaproteobacteria\n"
        "c3\tBacteria;uc_Bacteria\n",
    )
    out = tmp_path / "out.tsv"
    mod.main(["--contigs", str(contigs), "--taxonomy", str(tax), "--out", str(out)])
    lines = out.read_text().splitlines()
    assert lines[0] == "contigs\tpredictions"
    assert lines[1] == "c1\tBacteria;Pseudomonadota;Gammaproteobacteria"
    assert lines[2] == "c2\tBacteria;Bacillota;Bacillota_Gammaproteobacteria"
    assert lines[3] == "c3\tBacteria"
    assert lines[4] == "c4\t"


# ---------------------------------------------------------------------------
def _binette_group(tmp_path, group, bins, layout="1.1"):
    """Binette output; layout 1.1: bin_<id>.fa + column bin_id; 1.2: binette_bin<i>.fa + column name."""
    bindir = tmp_path / group / "final_bins"
    if layout == "1.1":
        report = ["bin_id\torigin\tname\tcompleteness\tcontamination\tscore\tsize\tN50\tcontig_count"]
    else:
        report = ["name\torigin\tis_original\toriginal_name\tcompleteness\tcontamination\tscore\t"
                  "checkm2_model\tsize\tN50\tcoding_density\tcontig_count"]
    for bin_id, (comp, cont) in bins.items():
        if layout == "1.1":
            write(bindir / f"bin_{bin_id}.fa", f">{group}_c{bin_id}\n{'A' * 100}\n")
            report.append(f"{bin_id}\tmetabat2\tbin.{bin_id}.fa\t{comp}\t{cont}\t0\t100\t100\t1")
        else:
            write(bindir / f"binette_bin{bin_id}.fa", f">{group}_c{bin_id}\n{'A' * 100}\n")
            report.append(f"binette_bin{bin_id}\tmetabat2\tTrue\tbin.{bin_id}.fa\t{comp}\t{cont}\t0\t"
                          f"Gradient Boost (General Model)\t100\t100\t0.9\t1")
    rep = write(tmp_path / group / "final_bins_quality_reports.tsv", "\n".join(report) + "\n")
    return bindir, rep


def test_candidates_from_binette_12_layout(tmp_path):
    cand = load("collect_candidates")
    bins, rep = _binette_group(tmp_path, "g", {1: (95.0, 1.0), 2: (60.0, 4.0)}, layout="1.2")
    out = tmp_path / "c"
    cand.main(["--groups", "g", "--bins", str(bins), "--reports", str(rep), "--report-format", "binette",
               "--out-dir", str(out), "--genome-info", str(tmp_path / "gi.csv"), "--table", str(tmp_path / "t.tsv")])
    rows = {r["candidate"]: r for r in csv.DictReader((tmp_path / "t.tsv").open(), delimiter="\t")}
    assert sorted(rows) == ["g__binette_bin1.fa", "g__binette_bin2.fa"]
    assert rows["g__binette_bin1.fa"]["completeness"] == "95.0"
    tri = load("anvio_triage")
    o = {k: tmp_path / f"{k}.out" for k in ("triage", "all", "refine", "howto")}
    tri.main(["--bins-dir", str(bins), "--report", str(rep), "--keep-min-completeness", "50",
              "--keep-max-contamination", "10", "--triage", str(o["triage"]), "--all-table", str(o["all"]),
              "--refine-table", str(o["refine"]), "--howto", str(o["howto"])])
    assert {r["bin"] for r in csv.DictReader(o["triage"].open(), delimiter="\t")} == {"binette_bin1", "binette_bin2"}


def test_candidates_and_final_names(tmp_path):
    cand = load("collect_candidates")
    fin = load("finalize_mags")
    a_bins, a_rep = _binette_group(tmp_path, "north", {10: (95.0, 1.0), 2: (60.0, 4.0)})
    b_bins, b_rep = _binette_group(tmp_path, "baltic", {5: (40.0, 1.0)})
    genomes = tmp_path / "cands"
    info = tmp_path / "info.csv"
    table = tmp_path / "cands.tsv"
    cand.main(["--groups", "north", "baltic", "--bins", str(a_bins), str(b_bins),
               "--reports", str(a_rep), str(b_rep), "--report-format", "binette",
               "--out-dir", str(genomes), "--genome-info", str(info), "--table", str(table)])
    assert sorted(p.name for p in genomes.iterdir()) == [
        "baltic__bin_5.fa", "north__bin_10.fa", "north__bin_2.fa"]
    rows = list(csv.DictReader(info.open()))
    assert {r["genome"]: float(r["completeness"]) for r in rows}["north__bin_10.fa"] == 95.0

    # pretend dRep kept all three, each its own cluster
    write(tmp_path / "Wdb.csv", "genome,score,cluster\n" + "".join(
        f"{g},1,{i}_1\n" for i, g in enumerate(["north__bin_10.fa", "north__bin_2.fa", "baltic__bin_5.fa"], 1)))
    write(tmp_path / "Cdb.csv", "genome,secondary_cluster,primary_cluster\n" + "".join(
        f"{g},{i}_1,{i}\n" for i, g in enumerate(["north__bin_10.fa", "north__bin_2.fa", "baltic__bin_5.fa"], 1)))
    out_dir = tmp_path / "mags"
    mags = tmp_path / "mags.tsv"
    fin.main(["--winners", str(tmp_path / "Wdb.csv"), "--clusters", str(tmp_path / "Cdb.csv"),
              "--candidates", str(table), "--genomes-dir", str(genomes), "--project", "proj",
              "--out-dir", str(out_dir), "--table", str(mags)])
    result = list(csv.DictReader(mags.open(), delimiter="\t"))
    # baltic bin (40 % complete) is filtered; north bins are named in natural bin order
    assert [(r["mag_id"], r["source_bin"]) for r in result] == [
        ("proj_MAG_1", "bin_2.fa"), ("proj_MAG_2", "bin_10.fa")]
    assert [r["quality_checkm2"] for r in result] == ["medium-quality", "near-complete"]
    assert (out_dir / "proj_MAG_2.fna").read_text().startswith(">north_c10")


def test_anvio_gene_calls(tmp_path):
    mod = load("anvio_gene_calls")
    gff = write(tmp_path / "g.gff",
                "##gff-version 3\n"
                "k125_7\tProdigal_v2.6.3\tCDS\t3\t95\t10.0\t+\t0\tID=1_1;partial=00;\n"
                "k125_7\tProdigal_v2.6.3\tCDS\t100\t300\t10.0\t-\t0\tID=1_2;partial=01;\n")
    faa = write(tmp_path / "p.faa", ">k125_7_1 # 3 # 95 # 1 # ID=1_1\nMKV*\n>k125_7_2 # 100 # 300 # -1\nMAA\n")
    out = tmp_path / "calls.tsv"
    mod.main(["--gff", str(gff), "--faa", str(faa), "--out", str(out)])
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert [(r["start"], r["stop"], r["direction"], r["partial"], r["aa_sequence"]) for r in rows] == [
        ("2", "95", "f", "0", "MKV"), ("99", "300", "r", "1", "MAA")]


def test_mag_summary_joins_tables(tmp_path):
    mod = load("mag_summary")
    mags = write(tmp_path / "mags.tsv", "mag_id\tgroup\tcompleteness\tcontamination\tn_contigs\nm1\tg\t95\t1\t10\nm2\tg\t60\t2\t20\n")
    gtdb = write(tmp_path / "gtdb.tsv", "user_genome\tclassification\nm1\td__Bacteria;p__X\n")
    mimag = write(tmp_path / "mimag.tsv", "mag_id\tmimag_quality\nm1\thigh-quality\nm2\tmedium-quality\n")
    cov = write(tmp_path / "cov.tsv", "mag_id\tgroup\tn_contigs\tgenome_coverage\nm1\tg\t10\t12.5\nm2\tg\t20\t3.1\n")
    out = tmp_path / "summary.tsv"
    mod.main(["--mags", str(mags), "--gtdbtk", str(gtdb), "--mimag", str(mimag),
              "--coverage", str(cov), "--out", str(out)])
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert rows[0]["gtdb_classification"] == "d__Bacteria;p__X"
    assert rows[1]["gtdb_classification"] == ""
    assert [r["mimag_quality"] for r in rows] == ["high-quality", "medium-quality"]
    assert rows[0]["genome_coverage"] == "12.5"
    # columns present in several inputs (n_contigs in mags.tsv and coverage.tsv) appear once
    assert not [c for c in rows[0] if c.endswith(("_x", "_y"))]


def test_quality_class_boundaries():
    fin = load("finalize_mags")
    assert fin.quality_class(90.0, 5.0) == "near-complete"
    assert fin.quality_class(89.99, 1.0) == "medium-quality"
    assert fin.quality_class(95.0, 5.01) == "medium-quality"
    assert fin.quality_class(50.0, 9.99) == "medium-quality"
    assert fin.quality_class(50.0, 10.0) == "low-quality"
    assert fin.quality_class(49.9, 0.0) == "low-quality"


def _triage(tmp_path, quality, review="contaminated", refine_comp=90, refine_cont=5):
    tri = load("anvio_triage")
    bins, rep = _binette_group(tmp_path, "g", quality)
    out = {k: tmp_path / f"{k}.out" for k in ("triage", "all", "refine", "howto")}
    tri.main(["--bins-dir", str(bins), "--report", str(rep), "--keep-min-completeness", "50",
              "--keep-max-contamination", "10", "--refine-min-completeness", str(refine_comp),
              "--refine-contamination-above", str(refine_cont), "--review", review, "--group", "g",
              "--triage", str(out["triage"]), "--all-table", str(out["all"]),
              "--refine-table", str(out["refine"]), "--howto", str(out["howto"])])
    rows = {r["bin"]: r["decision"] for r in csv.DictReader(out["triage"].open(), delimiter="\t")}
    return rows, out, bins


def test_anvio_triage_decisions(tmp_path):
    quality = {1: (96.0, 2.0),    # near-complete -> keep
               2: (95.0, 7.0),    # passes 50/10 but worth cleaning below 5 % -> refine
               3: (92.0, 25.0),   # complete but contaminated -> refine
               4: (85.0, 25.0),   # contaminated, below 90 % -> drop
               5: (70.0, 4.0),    # medium quality -> keep
               6: (40.0, 1.0)}    # too incomplete -> drop
    rows, out, _ = _triage(tmp_path, quality)
    assert rows == {"bin_1": "keep", "bin_2": "refine", "bin_3": "refine", "bin_4": "drop",
                    "bin_5": "keep", "bin_6": "drop"}
    # collections: every bin in 'Binette', only the refine bins in 'Refine'
    all_bins = {line.split("\t")[1].strip() for line in out["all"].read_text().splitlines()}
    refine_bins = {line.split("\t")[1].strip() for line in out["refine"].read_text().splitlines()}
    assert all_bins == {f"bin_{i}" for i in quality}
    assert refine_bins == {"bin_2", "bin_3"}
    howto = out["howto"].read_text()
    assert "-C Refine -b bin_3" in howto and "-b bin_1`" not in howto
    # most complete bins first in the guide
    assert howto.index("bin_2 |") < howto.index("bin_3 |")


def test_anvio_triage_thresholds_and_review_all(tmp_path):
    quality = {1: (96.0, 2.0), 2: (85.0, 25.0)}
    rows, _, _ = _triage(tmp_path / "a", quality, refine_comp=80)
    assert rows == {"bin_1": "keep", "bin_2": "refine"}
    rows, out, _ = _triage(tmp_path / "b", quality, review="all")
    assert rows == {"bin_1": "review", "bin_2": "review"}
    assert len(out["refine"].read_text().splitlines()) == len(out["all"].read_text().splitlines())


def test_anvio_candidates_merges_kept_and_refined(tmp_path):
    cands = load("anvio_candidates")
    rows, out, bins = _triage(tmp_path, {1: (96.0, 2.0), 2: (92.0, 25.0), 3: (30.0, 0.0)})
    refined = tmp_path / "refined"
    write(refined / "bin_2.fa", ">g_c2\nACGT\n")
    write(refined / "Bin_2_b.fa", ">g_c2b\nACGT\n")
    report = write(tmp_path / "checkm2.tsv", "Name\tCompleteness\tContamination\n"
                   "bin_2\t91.0\t3.0\nBin_2_b\t20.0\t1.0\n")
    outdir, quality = tmp_path / "cands", tmp_path / "quality.tsv"
    cands.main(["--triage", str(out["triage"]), "--binette-bins", str(bins), "--refined-bins", str(refined),
                "--refined-checkm2", str(report), "--out-dir", str(outdir), "--quality", str(quality)])
    got = {r["bin"]: (r["completeness"], r["origin"]) for r in csv.DictReader(quality.open(), delimiter="\t")}
    assert got == {"bin_1": ("96.0", "metabat2"), "refined_bin_2": ("91.0", "anvio_refined"),
                   "refined_Bin_2_b": ("20.0", "anvio_refined")}
    assert sorted(p.name for p in outdir.iterdir()) == ["bin_1.fa", "refined_Bin_2_b.fa", "refined_bin_2.fa"]
    # the simple table feeds collect_candidates
    coll = load("collect_candidates")
    coll.main(["--groups", "g", "--bins", str(outdir), "--reports", str(quality), "--report-format", "simple",
               "--out-dir", str(tmp_path / "all"), "--genome-info", str(tmp_path / "gi.csv"),
               "--table", str(tmp_path / "c.tsv")])
    assert sorted(p.name for p in (tmp_path / "all").iterdir()) == [
        "g__bin_1.fa", "g__refined_Bin_2_b.fa", "g__refined_bin_2.fa"]


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-q"]))


# MMseqs2 18 output for a Swiss-Prot database built with an NCBI taxdump of 2025: the
# requested rank 'superkingdom' no longer exists (now 'domain'), so the --lca-ranks lineage
# (column 9) starts with 'unknown'; the full lineage (column 10) still names the domain.
MMSEQS_2025 = (
    "k125_16777\t400668\tspecies\tMarinomonas sp. MWYL1\t1\t1\t1\t1.000\t"
    "unknown;Pseudomonadota;Gammaproteobacteria;Oceanospirillales;Oceanospirillaceae;Marinomonas;"
    "Marinomonas sp. MWYL1\t-_cellular organisms;-_Bacteria;k_Pseudomonadati;p_Pseudomonadota;"
    "c_Gammaproteobacteria;o_Oceanospirillales;f_Oceanospirillaceae;g_Marinomonas;"
    "-_unclassified Marinomonas;s_Marinomonas sp. MWYL1\n"
    "k125_17977\t1236\tclass\tGammaproteobacteria\t1\t1\t1\t1.000\t"
    "unknown;Pseudomonadota;Gammaproteobacteria;uc_Gammaproteobacteria;uc_Gammaproteobacteria;"
    "uc_Gammaproteobacteria;uc_Gammaproteobacteria\t-_cellular organisms;-_Bacteria;"
    "k_Pseudomonadati;p_Pseudomonadota;c_Gammaproteobacteria\n"
    "k125_542\t2759\tdomain\tEukaryota\t1\t1\t1\t1.000\t"
    "unknown;unknown;unknown;unknown;unknown;unknown;unknown\t-_cellular organisms;-_Eukaryota\n"
    "k125_3110\t1\tno rank\troot\t1\t1\t1\t1.000\t"
    "unknown;unknown;unknown;unknown;unknown;unknown;unknown\t-_root\n"
)


def _taxconverter_table(mmseqs_text):
    rows = ["contigs\tpredictions"]
    for line in mmseqs_text.splitlines():
        parts = line.split("\t")
        rows.append(f"{parts[0]}\t{parts[8]}")
    return "\n".join(rows) + "\n"


def _read_predictions(path):
    with open(path) as fh:
        fh.readline()
        return dict(line.rstrip("\n").split("\t") for line in fh)


def test_fix_vamb_taxonomy_recovers_domain_from_full_lineage(tmp_path, capsys):
    mod = load("fix_vamb_taxonomy")
    contigs = write(tmp_path / "c.fa", "".join(f">{c}\nA\n" for c in
                                               ["k125_16777", "k125_17977", "k125_542", "k125_3110"]))
    mmseqs = write(tmp_path / "mmseqs.tsv", MMSEQS_2025)
    tax = write(tmp_path / "taxconverter.tsv", _taxconverter_table(MMSEQS_2025))
    out = tmp_path / "out.tsv"
    mod.main(["--contigs", str(contigs), "--taxonomy", str(tax), "--mmseqs", str(mmseqs),
              "--out", str(out)])
    pred = _read_predictions(out)
    assert pred["k125_16777"] == ("Bacteria;Pseudomonadota;Gammaproteobacteria;Oceanospirillales;"
                                  "Oceanospirillaceae;Marinomonas;Marinomonas sp. MWYL1")
    assert pred["k125_17977"] == "Bacteria;Pseudomonadota;Gammaproteobacteria"
    assert pred["k125_542"] == "Eukaryota"
    assert pred["k125_3110"] == ""
    assert "Domains taken from the MMseqs2 full lineage: 3" in capsys.readouterr().out


def test_fix_vamb_taxonomy_keeps_superkingdom_lineages():
    mod = load("fix_vamb_taxonomy")
    assert mod.domain_from_full_lineage("-_cellular organisms;d_Bacteria;p_Pseudomonadota") == "Bacteria"
    assert mod.domain_from_full_lineage("-_root") is None
    assert mod.domain_from_full_lineage("-_cellular organisms;-_Opisthokonta") is None
    lineage, filled = mod.fill_domain("Bacteria;Pseudomonadota", "Bacteria")
    assert (lineage, filled) == ("Bacteria;Pseudomonadota", False)


def test_fix_vamb_taxonomy_fails_without_any_lineage(tmp_path):
    mod = load("fix_vamb_taxonomy")
    contigs = write(tmp_path / "c.fa", ">k125_16777\nA\n>k125_3110\nA\n")
    tax = write(tmp_path / "taxconverter.tsv", _taxconverter_table(MMSEQS_2025))
    with pytest.raises(SystemExit) as err:
        # without --mmseqs the 2025-style lineages all start with 'unknown'
        mod.main(["--contigs", str(contigs), "--taxonomy", str(tax), "--out", str(tmp_path / "o.tsv")])
    assert "no contig has a lineage" in str(err.value)


def test_tables_have_unix_line_ends(tmp_path):
    # csv writers default to CRLF; shell steps reading the last column would see a trailing '\r'
    _, out, _ = _triage(tmp_path, {1: (96.0, 2.0), 2: (95.0, 7.0)})
    assert b"\r" not in out["triage"].read_bytes()
    scripts = Path(__file__).resolve().parents[1] / "workflow" / "scripts"
    for script in sorted(scripts.glob("*.py")):
        if script.name == "graphbin2_fast.py":  # writes CRLF like upstream GraphBin2; the rule strips it
            continue
        text = script.read_text()
        for call in ("csv.writer(", "csv.DictWriter("):
            start = 0
            while (i := text.find(call, start)) >= 0:
                depth, j = 0, i + len(call) - 1
                while True:
                    depth += {"(": 1, ")": -1}.get(text[j], 0)
                    if depth == 0:
                        break
                    j += 1
                assert "lineterminator" in text[i:j], f"{script.name}: {text[i:j][:60]}"
                start = j
