import json
import math

import pytest

import mag_coverage as mc
import mimag_rna as mr
from conftest import GFF_FULL, GFF_PARTIAL16, AA20, trna_table
from ena_common import read_tsv


# ----------------------------------------------------------------------------- mimag_rna
def test_barrnap_counts(tmp_path):
    p = tmp_path / "a.gff"
    p.write_text(GFF_FULL)
    c = mr.parse_barrnap_gff(p)
    assert (c["5S"], c["16S"], c["23S"], c["partial"]) == (1, 1, 1, 0)
    p.write_text(GFF_PARTIAL16)
    assert mr.parse_barrnap_gff(p)["16S"] == 1
    assert mr.parse_barrnap_gff(p, exclude_partial=True)["16S"] == 0


def test_barrnap_eukaryotic_models_do_not_count(tmp_path):
    p = tmp_path / "e.gff"
    p.write_text("##gff-version 3\nc1\tbarrnap:0.9\trRNA\t1\t150\t0\t+\t.\tName=5_8S_rRNA;product=5.8S ribosomal RNA\n"
                 "c1\tbarrnap:0.9\trRNA\t200\t2000\t0\t+\t.\tName=18S_rRNA;product=18S ribosomal RNA\n")
    c = mr.parse_barrnap_gff(p)
    assert (c["5S"], c["16S"], c["23S"], c["other"]) == (0, 0, 0, 2)


def test_trna_standard20_and_legacy(tmp_path):
    p = tmp_path / "t.txt"
    p.write_text(trna_table(AA20 + ["SeC", "Undet", "Leu"]))
    r = mr.parse_trnascan(p)
    assert r["n_unique"] == 20 and "Ile" in r["isotypes"] and "Met" in r["isotypes"]
    assert r["n_trna"] == 23
    p.write_text(trna_table(AA20, pseudo=("Trp", "Tyr", "Val")))
    assert mr.parse_trnascan(p)["n_unique"] == 17
    assert mr.parse_trnascan(p)["n_pseudo"] == 3
    # legacy mode (old ENA-prep Snakefile): pseudo genes and non-standard types counted
    assert mr.parse_trnascan(p, mode="legacy")["n_unique"] == 20
    p.write_text(trna_table(AA20 + ["SeC", "Sup"]))
    assert mr.parse_trnascan(p, mode="legacy")["n_unique"] == 22  # 'Fmet' instead of Met, plus Sec and Sup


@pytest.mark.parametrize("comp,cont,trna,expect_mimag,expect_gu", [
    (95.0, 1.0, 20, "high-quality", "HQ"),
    (90.0, 1.0, 20, "near-complete", "HQ"),   # MIMAG needs > 90, genome_uploader >= 90 (conflict)
    (95.0, 5.0, 20, "near-complete", "HQ"),   # MIMAG needs < 5, genome_uploader <= 5 (conflict)
    (95.0, 1.0, 17, "near-complete", "MQ"),   # tRNA criterion not met
    (95.0, 7.0, 20, "medium-quality", "MQ"),
    (70.0, 3.0, 20, "medium-quality", "MQ"),
    (49.9, 1.0, 20, "low-quality", "MQ"),
    (80.0, 10.0, 20, "low-quality", "MQ"),    # contamination >= 10 fails MIMAG
])
def test_quality_boundaries(comp, cont, trna, expect_mimag, expect_gu):
    q, _note = mr.mimag_quality(comp, cont, True, True, True, trna)
    assert q == expect_mimag
    assert mr.genome_uploader_quality(comp, cont, mr.rrna_presence(True, True, True, trna)) == expect_gu


def test_mimag_cli_contract_columns(mag_set, tmp_path):
    out = tmp_path / "mimag_rna.tsv"
    rc = mr.main(["--barrnap-dir", str(mag_set["barrnap"]), "--trnascan-dir", str(mag_set["trnascan"]),
                  "--mags-table", str(mag_set["mags"]), "--out", str(out), "--gtdbtk", str(mag_set["gtdbtk"])])
    assert rc == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    for col in ("mag_id", "rrna_5S", "rrna_16S", "rrna_23S", "trna_unique", "rRNA_presence", "mimag_quality"):
        assert col in rows["g1_MAG_1"]
    m1, m2, m3 = rows["g1_MAG_1"], rows["g1_MAG_2"], rows["g1_MAG_3"]
    assert (m1["rrna_5S"], m1["rrna_16S"], m1["rrna_23S"], m1["trna_unique"]) == ("1", "1", "1", "20")
    assert m1["rRNA_presence"] == "True" and m1["mimag_quality"] == "high-quality"
    assert m1["genome_uploader_quality"] == "HQ" and m1["ena_quality_conflict"] == "False"
    assert m2["trna_unique"] == "17" and m2["rRNA_presence"] == "False" and m2["mimag_quality"] == "near-complete"
    assert "17 tRNA isotypes" in m2["note"]
    assert m3["rrna_16S"] == "0" and m3["mimag_quality"] == "medium-quality" and m3["domain"] == "Archaea"


def test_mimag_old_layout_and_missing(tmp_path):
    b, t = tmp_path / "rrna", tmp_path / "trna"
    b.mkdir()
    t.mkdir()
    (b / "north_MAG_1.fna.gff").write_text(GFF_FULL)
    (t / "north_MAG_1.fna.txt").write_text(trna_table(AA20))
    q = tmp_path / "quality_report.tsv"
    q.write_text("Name\tCompleteness\tContamination\nnorth_MAG_1\t95.05\t0.01\nnorth_MAG_2\t60\t1\n")
    out = tmp_path / "o.tsv"
    assert mr.main(["--barrnap-dir", str(b), "--trnascan-dir", str(t), "--checkm2", str(q), "--out", str(out)]) == 1
    assert mr.main(["--barrnap-dir", str(b), "--trnascan-dir", str(t), "--checkm2", str(q), "--out", str(out),
                    "--allow-missing"]) == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    assert rows["north_MAG_1"]["rRNA_presence"] == "True"
    assert "no barrnap file" in rows["north_MAG_2"]["note"]


# ----------------------------------------------------------------------------- mag_coverage
def test_mag_coverage_contract(mag_set, tmp_path):
    out = tmp_path / "coverage.tsv"
    rc = mc.main(["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]),
                  "--coverm", f"g1={mag_set['coverm']}", "--out", str(out)])
    assert rc == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    m1 = rows["g1_MAG_1"]
    assert m1["group"] == "g1"
    assert float(m1["genome_coverage"]) == pytest.approx(20.0)          # mean of per-contig sums (16, 24)
    assert float(m1["mean_depth_S1"]) == pytest.approx(70000 / 6000, rel=1e-4)   # length-weighted
    assert float(m1["mean_depth_S2"]) == pytest.approx(25000 / 6000, rel=1e-4)
    assert float(m1["mean_depth_S3"]) == pytest.approx(1.5)
    assert float(rows["g1_MAG_2"]["genome_coverage"]) == pytest.approx(9.0)
    assert float(rows["g1_MAG_3"]["genome_coverage"]) == pytest.approx(2.0)


def test_mag_coverage_length_weighted_and_headerless(mag_set, tmp_path):
    lines = mag_set["coverm"].read_text().splitlines()[1:]
    headerless = tmp_path / "abundance.tsv"
    headerless.write_text("\n".join(lines) + "\n")
    out = tmp_path / "c.tsv"
    assert mc.main(["--mags-dir", str(mag_set["fasta"]), "--coverm", str(headerless), "--out", str(out),
                    "--method", "length_weighted"]) == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    assert float(rows["g1_MAG_1"]["genome_coverage"]) == pytest.approx((70000 + 25000 + 9000) / 6000, rel=1e-4)
    assert "mean_depth_sample_1" in rows["g1_MAG_1"]


def test_mag_coverage_missing_contigs_fail(mag_set, tmp_path):
    t = tmp_path / "partial.tsv"
    t.write_text("Contig\tS1 Mean\nk127_1\t3\n")
    assert mc.main(["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]),
                    "--coverm", f"g1={t}", "--out", str(tmp_path / "x.tsv")]) == 1


def test_mag_coverage_requires_mean_columns(tmp_path):
    t = tmp_path / "ra.tsv"
    t.write_text("Genome\tS1 Relative Abundance (%)\nMAG\t3\n")
    with pytest.raises(Exception):
        mc.read_header(t)


def test_assembly_coverage(mag_set, capsys):
    assert mc.main(["assembly", "--coverm", str(mag_set["coverm"])]) == 0
    res = json.loads(capsys.readouterr().out)
    assert res["coverage"] == pytest.approx(60.0) and res["n_contigs"] == 6 and res["n_samples"] == 3
    with pytest.raises(Exception):
        mc.assembly_coverage(mag_set["coverm"], "length_weighted", None)
    assert not math.isnan(mc.genome_depth({"a": 10}, {"a": [1.0, 2.0]}, 2)["genome_coverage"])
