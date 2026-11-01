import gzip
import json

import pytest

import ena_genome_table as egt
import gtdb_to_ncbi as g2n
import mag_coverage as mc
import mimag_rna as mr
from conftest import RUNS
from ena_common import iter_fasta, read_tsv

PARAMS = {"ena": {"genome_type": "bins", "metagenome": "seawater metagenome",
                  "broad_environment": "marine biome [ENVO:00000447]",
                  "local_environment": "sea surface microlayer [ENVO:03600079]",
                  "environmental_medium": "sea surface microlayer [ENVO:03600079]",
                  "assembly_software": "MEGAHIT_v1.2.9", "binning_software": "Binette_v1.1.2",
                  "binning_parameters": "default", "stats_generation_software": "CheckM2_v1.1.0",
                  "groups": {"g2": {"local_environment": "underlying water"}}}}
MV = ("Genome ID\tGTDB classification\tMajority vote NCBI classification\n"
      "g1_MAG_1\tx\td__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;f__Flavobacteriaceae;"
      "g__Polaribacter;s__\n"
      "g1_MAG_2\tx\td__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;f__;g__;s__\n"
      "g1_MAG_3\tx\td__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;f__Nitrosopumilaceae;"
      "g__Nitrosopumilus;s__\n")


@pytest.fixture
def upstream(mag_set, tmp_path):
    r = tmp_path / "up"
    r.mkdir()
    assert mr.main(["--barrnap-dir", str(mag_set["barrnap"]), "--trnascan-dir", str(mag_set["trnascan"]),
                    "--mags-table", str(mag_set["mags"]), "--out", str(r / "mimag_rna.tsv")]) == 0
    assert mc.main(["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]),
                    "--coverm", f"g1={mag_set['coverm']}", "--out", str(r / "coverage.tsv")]) == 0
    (r / "mv.tsv").write_text(MV)
    assert g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(r / "mv.tsv"),
                     "--out", str(r / "ncbi_lineage.tsv")]) == 0
    (r / "params.json").write_text(json.dumps(PARAMS))
    return r


def table_args(mag_set, up, out, extra=()):
    return ["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]), "--rna", str(up / "mimag_rna.tsv"),
            "--coverage", str(up / "coverage.tsv"), "--lineage", str(up / "ncbi_lineage.tsv"),
            "--samples", str(mag_set["samples"]), "--params", str(up / "params.json"),
            "--fasta-outdir", str(out.parent / "fasta"), "--out", str(out)] + list(extra)


def test_table_exact_columns_and_values(mag_set, upstream, tmp_path):
    out = tmp_path / "ena" / "genome_uploader_input.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 0
    header = out.read_text().splitlines()[0].split("\t")
    assert header == egt.GENOME_UPLOADER_COLUMNS
    assert len(header) == 17
    rows = {r["genome_name"]: r for r in read_tsv(out)}
    m1 = rows["g1_MAG_1"]
    assert m1["accessions"] == ",".join(RUNS) and m1["co-assembly"] == "True"
    assert m1["completeness"] == "95.05" and m1["contamination"] == "0.01"
    assert m1["genome_coverage"] == "20" and m1["rRNA_presence"] == "True"
    assert m1["NCBI_lineage"].startswith("d__Bacteria;p__Bacteroidota") and m1["metagenome"] == "seawater metagenome"
    assert rows["g1_MAG_2"]["rRNA_presence"] == "False"
    gz = m1["genome_path"]
    assert gz.endswith("g1_MAG_1.fa.gz")
    with gzip.open(gz, "rt") as fh:
        assert fh.read(1) == ">"
    assert [s for _i, _h, s in iter_fasta(gz)] == [s for _n, s in mag_set["contigs"]["g1_MAG_1"]]
    prov = read_tsv(str(out.with_suffix("")) + ".provenance.tsv")
    assert prov[0]["mimag_quality"] == "high-quality" and prov[0]["n_runs"] == "3"


def test_genome_uploader_accepts_table(mag_set, upstream, tmp_path):
    pytest.importorskip("genomeuploader")
    import ena_submit
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 0
    assert ena_submit.genome_uploader_validation(out, "bins") == []
    assert ena_submit.genome_uploader_validation(out, "mags") == []
    assert ena_submit.extra_table_checks(out) == []


def test_metagenome_checked_against_genome_uploader_list(mag_set, upstream, tmp_path):
    pytest.importorskip("genomeuploader")
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out, ["--metagenome", "sea surface microlayer metagenome"])) == 1
    assert not out.exists()


def test_validation_errors_are_collected(mag_set, upstream, tmp_path, caplog):
    out = tmp_path / "ena" / "t.tsv"
    p = json.loads((upstream / "params.json").read_text())
    p["ena"]["binning_software"] = ""
    (upstream / "params.json").write_text(json.dumps(p))
    rc = egt.main(table_args(mag_set, upstream, out, ["--name-prefix", "a_very_long_project_prefix_for_names_"]))
    assert rc == 1 and not out.exists()
    msg = "\n".join(caplog.messages)
    assert "binning_software is empty" in msg and "longer than 39 characters" in msg


def test_mag_tier_thresholds(mag_set, upstream, tmp_path, caplog):
    rows = mag_set["mags"].read_text().replace("g1_MAG_3\tc3\tg1\tb3\tvamb\t62.4\t3.1", "g1_MAG_3\tc3\tg1\tb3\tvamb\t62.4\t12.5")
    mag_set["mags"].write_text(rows)
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out, ["--genome-type", "mags"])) == 1
    assert "fails MAG criteria" in "\n".join(caplog.messages)
    assert egt.main(table_args(mag_set, upstream, out, ["--genome-type", "bins"])) == 0


def test_single_run_group_is_not_coassembly(mag_set, upstream, tmp_path):
    s = mag_set["samples"]
    lines = s.read_text().splitlines()
    s.write_text("\n".join([lines[0], lines[1]] + [l.replace("\tg1\t", "\tg9\t") for l in lines[2:]]) + "\n")
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 0
    r = read_tsv(out)[0]
    assert r["accessions"] == RUNS[0] and r["co-assembly"] == "False"


def test_fasta_rule_violation(mag_set, upstream, tmp_path, caplog):
    f = mag_set["fasta"] / "g1_MAG_3.fna"
    f.write_text(f.read_text().rstrip("\n") + "NNNN\n")
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 1
    assert "terminal N" in "\n".join(caplog.messages)


def test_missing_run_accession(mag_set, upstream, tmp_path, caplog):
    s = mag_set["samples"]
    s.write_text(s.read_text().replace("ERR0000002", ""))
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 1
    assert "run_accession missing" in "\n".join(caplog.messages)


def test_assembly_quality_conflict_at_boundary(mag_set, upstream, tmp_path, caplog):
    rows = mag_set["mags"].read_text().replace("\t95.05\t0.01\t", "\t90.00\t0.01\t")
    mag_set["mags"].write_text(rows)
    out = tmp_path / "ena" / "t.tsv"
    assert egt.main(table_args(mag_set, upstream, out)) == 1
    assert "MIMAG HQ needs > 90 and < 5" in "\n".join(caplog.messages)
    assert egt.main(table_args(mag_set, upstream, out, ["--allow-assembly-quality-conflict"])) == 0
