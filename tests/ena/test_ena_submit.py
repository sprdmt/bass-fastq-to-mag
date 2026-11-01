import json
import os
import stat
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import pytest

import ena_coassembly as ec
import ena_genome_table as egt
import ena_submit as es
import gtdb_to_ncbi as g2n
import mag_coverage as mc
import mimag_rna as mr
from ena_common import read_tsv

FAKE_GENOME_UPLOAD = r'''#!{py}
import argparse, csv, json, os, sys
from pathlib import Path
ap = argparse.ArgumentParser()
for o in ("-u", "--genome_info", "--out", "--centre_name", "--test-suffix"):
    ap.add_argument(o)
for o in ("--mags", "--bins", "--private", "--live", "--force"):
    ap.add_argument(o, action="store_true")
a = ap.parse_args()
Path(os.environ["FAKE_LOG"]).open("a").write(json.dumps({{"tool": "genome_upload", "argv": sys.argv[1:],
    "cwd": os.getcwd()}}) + "\n")
assert os.environ.get("ENA_WEBIN") and os.environ.get("ENA_WEBIN_PASSWORD")
up = Path(a.out) / ("MAG_upload" if a.mags else "bin_upload")
md = up / ("manifests" if a.live else "manifests_test")
md.mkdir(parents=True, exist_ok=True)
reg = []
for i, r in enumerate(csv.DictReader(open(a.genome_info), delimiter="\t")):
    alias = r["genome_name"] + ("" if a.live else "_" + a.test_suffix)
    acc = "ERS%08d" % (70000000 + i)
    reg.append(alias + "\t" + acc)
    (md / (alias + ".manifest")).write_text("STUDY\t%s\nSAMPLE\t%s\nASSEMBLYNAME\t%s\nASSEMBLY_TYPE\tbinned metagenome\n"
        "COVERAGE\t%s\nPROGRAM\tMEGAHIT\nPLATFORM\tIllumina NovaSeq 6000\nMOLECULETYPE\tgenomic DNA\nRUN_REF\t%s\n"
        "FASTA\t%s\n" % (a.u, acc, alias, r["genome_coverage"], r["accessions"], r["genome_path"]))
name = ("registered_MAGs" if a.mags else "registered_bins") + ("" if a.live else "_test") + ".tsv"
(up / name).write_text("\n".join(reg) + "\n")
print("DEBUG leaked password=" + os.environ["ENA_WEBIN_PASSWORD"])
'''

FAKE_WEBIN = r'''#!{py}
import argparse, json, os, sys, zlib
from pathlib import Path
if "--help" in sys.argv:
    print("usage: webin_cli_handler [-m M] [-c C] [--mode M] [--test] [--java-heap-size-max N]")
    sys.exit(0)
ap = argparse.ArgumentParser()
for o in ("-m", "-c", "--mode", "--java-heap-size-initial", "--java-heap-size-max"):
    ap.add_argument(o)
ap.add_argument("--test", action="store_true")
a = ap.parse_args()
Path(os.environ["FAKE_LOG"]).open("a").write(json.dumps({{"tool": "webin", "argv": sys.argv[1:]}}) + "\n")
m = Path(a.m)
name = [l.split("\t")[1].strip() for l in m.read_text().splitlines() if l.startswith("ASSEMBLYNAME")][0]
fail = os.environ.get("FAKE_WEBIN_FAIL", "")
if a.mode == fail:
    (m.parent / "webin-cli.report").write_text("ERROR: something is wrong\n")
    sys.exit(1)
if a.mode == "validate":
    (m.parent / "webin-cli.report").write_text("INFO : Submission(s) validated successfully.\n")
else:
    acc = "ERZ%08d" % (zlib.crc32(name.encode()) % 10**8)
    sub = m.parent / "genome" / name / "submit"
    sub.mkdir(parents=True, exist_ok=True)
    (sub / "receipt.xml").write_text('<RECEIPT success="true"><ANALYSIS accession="%s" alias="webin-genome-%s" '
                                     'status="PRIVATE"/></RECEIPT>' % (acc, name))
    (m.parent / "webin-cli.report").write_text(("This was a TEST submission(s).\n" if a.test else "") +
        "INFO : The submission has been completed successfully. The following analysis accession was assigned "
        "to the submission: %s\n" % acc)
'''


@pytest.fixture
def fakes(tmp_path, monkeypatch):
    b = tmp_path / "bin"
    b.mkdir()
    for name, text in (("genome_upload", FAKE_GENOME_UPLOAD), ("webin_cli_handler", FAKE_WEBIN)):
        p = b / name
        p.write_text(text.format(py=sys.executable))
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    log = tmp_path / "fake_calls.jsonl"
    monkeypatch.setenv("PATH", str(b) + os.pathsep + os.environ["PATH"])
    monkeypatch.setenv("FAKE_LOG", str(log))
    return log


def calls(log):
    return [json.loads(l) for l in log.read_text().splitlines()] if log.exists() else []


@pytest.fixture
def table(mag_set, tmp_path):
    up = tmp_path / "up"
    up.mkdir()
    assert mr.main(["--barrnap-dir", str(mag_set["barrnap"]), "--trnascan-dir", str(mag_set["trnascan"]),
                    "--mags-table", str(mag_set["mags"]), "--out", str(up / "rna.tsv")]) == 0
    assert mc.main(["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]),
                    "--coverm", f"g1={mag_set['coverm']}", "--out", str(up / "cov.tsv")]) == 0
    (up / "mv.tsv").write_text("Genome ID\tGTDB classification\tMajority vote NCBI classification\n")
    assert g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(up / "mv.tsv"),
                     "--out", str(up / "lin.tsv")]) == 0
    (up / "p.json").write_text(json.dumps({"ena": {
        "metagenome": "seawater metagenome", "broad_environment": "marine biome", "local_environment": "SML",
        "environmental_medium": "sea surface microlayer", "assembly_software": "MEGAHIT_v1.2.9",
        "binning_software": "Binette_v1.1.2", "binning_parameters": "default",
        "stats_generation_software": "CheckM2_v1.1.0"}}))
    out = tmp_path / "ena" / "genome_uploader_input.tsv"
    assert egt.main(["--mags-table", str(mag_set["mags"]), "--mags-dir", str(mag_set["fasta"]), "--rna",
                     str(up / "rna.tsv"), "--coverage", str(up / "cov.tsv"), "--lineage", str(up / "lin.tsv"),
                     "--samples", str(mag_set["samples"]), "--params", str(up / "p.json"), "--genome-type", "bins",
                     "--fasta-outdir", str(tmp_path / "ena" / "fasta"), "--out", str(out)]) == 0
    return out


def genomes_args(table, wd, *extra):
    return ["genomes", "--genome-info", str(table), "--genome-type", "bins", "--study", "PRJEB00001",
            "--centre-name", "University of Vienna", "--workdir", str(wd), "--delay", "0"] + list(extra)


def all_text(directory):
    return "\n".join(p.read_text(errors="ignore") for p in Path(directory).rglob("*") if p.is_file())


def test_test_mode_end_to_end(table, fakes, fake_creds, tmp_path):
    wd = tmp_path / "wd"
    assert es.main(genomes_args(table, wd, "--test-suffix", "1700000000")) == 0
    c = calls(fakes)
    gu = [x for x in c if x["tool"] == "genome_upload"]
    assert len(gu) == 1
    assert "--live" not in gu[0]["argv"] and "--test-suffix" in gu[0]["argv"] and "--private" in gu[0]["argv"]
    assert "--bins" in gu[0]["argv"]
    web = [x for x in c if x["tool"] == "webin"]
    assert [w["argv"][w["argv"].index("--mode") + 1] for w in web] == ["validate", "submit"] * 3
    assert all("--test" in w["argv"] and "--java-heap-size-max" in w["argv"] for w in web)
    acc = read_tsv(wd / "ena_accessions_test.tsv")
    assert len(acc) == 3 and all(r["analysis_accession"].startswith("ERZ") for r in acc)
    assert all(r["alias"].endswith("_1700000000") for r in acc)
    state = json.loads((wd / "ena_submit_state.json").read_text())
    assert state["genomes"]["test"]["complete"] is True
    text = all_text(wd)
    assert fake_creds not in text and "leaked password=***" in text


def test_resume_skips_done_work(table, fakes, fake_creds, tmp_path):
    wd = tmp_path / "wd"
    assert es.main(genomes_args(table, wd)) == 0
    n = len(calls(fakes))
    assert es.main(genomes_args(table, wd)) == 0
    assert len(calls(fakes)) == n


def test_live_gate_and_credentials(table, fakes, fake_creds, tmp_path, monkeypatch):
    wd = tmp_path / "wd"
    assert es.main(genomes_args(table, wd, "--live")) == 1          # no TEST run recorded
    assert calls(fakes) == []
    assert es.main(genomes_args(table, wd)) == 0                    # TEST run
    n = len(calls(fakes))
    monkeypatch.delenv("ENA_WEBIN_PASSWORD")
    assert es.main(genomes_args(table, wd, "--live")) == 2          # credentials missing
    assert len(calls(fakes)) == n
    monkeypatch.setenv("ENA_WEBIN_PASSWORD", fake_creds)
    assert es.main(genomes_args(table, wd, "--live")) == 0
    new = calls(fakes)[n:]
    gu = [x for x in new if x["tool"] == "genome_upload"][0]
    assert "--live" in gu["argv"] and "--test-suffix" not in gu["argv"]
    assert all("--test" not in x["argv"] for x in new if x["tool"] == "webin")
    live = read_tsv(wd / "ena_accessions_live.tsv")
    assert {r["alias"] for r in live} == {"g1_MAG_1", "g1_MAG_2", "g1_MAG_3"}


def test_changed_table_is_refused(table, fakes, fake_creds, tmp_path):
    wd = tmp_path / "wd"
    assert es.main(genomes_args(table, wd)) == 0
    table.write_text(table.read_text().replace("MEGAHIT_v1.2.9", "MEGAHIT_v1.2.8"))
    assert es.main(genomes_args(table, wd)) == 1
    assert es.main(genomes_args(table, wd, "--live")) == 1
    assert es.main(genomes_args(table, wd, "--restart")) == 0


def test_validation_failure_stops_before_submit(table, fakes, fake_creds, tmp_path, monkeypatch):
    monkeypatch.setenv("FAKE_WEBIN_FAIL", "validate")
    wd = tmp_path / "wd"
    assert es.main(genomes_args(table, wd)) == 1
    modes = [x["argv"][x["argv"].index("--mode") + 1] for x in calls(fakes) if x["tool"] == "webin"]
    assert modes == ["validate"] * 3


def test_dry_run_needs_no_credentials(table, fakes, no_creds, tmp_path, capsys):
    assert es.main(genomes_args(table, tmp_path / "wd", "--dry-run")) == 0
    out = capsys.readouterr().out
    assert "[dry-run]" in out and "genome_upload" in out and "--live" not in out
    assert calls(fakes) == []


def test_validate_table_cli(table, tmp_path):
    pytest.importorskip("genomeuploader")
    assert es.main(["validate-table", "--genome-info", str(table), "--genome-type", "mags"]) == 0
    bad = tmp_path / "bad.tsv"
    bad.write_text(table.read_text().replace("True", "maybe", 1))
    assert es.main(["validate-table", "--genome-info", str(bad), "--genome-type", "mags"]) == 1


def test_coassembly_flow(mag_set, fakes, fake_creds, tmp_path, monkeypatch):
    params = tmp_path / "p.json"
    params.write_text(json.dumps({"ena": {"study": "PRJEB00001", "metagenome": "seawater metagenome",
                                          "metagenome_taxid": 1561972}}))
    xml_dir = tmp_path / "co"
    assert ec.main(["sample-xml", "--samples", str(mag_set["samples"]), "--params", str(params),
                    "--outdir", str(xml_dir)]) == 0
    posted = []

    def fake_post(url, body, headers, user, password, timeout=300):
        posted.append(url)
        alias = ET.fromstring(body).find(".//SAMPLE").get("alias")
        return 200, (f'<RECEIPT success="true"><SAMPLE accession="ERS99999999" alias="{alias}"><EXT_ID '
                     f'accession="SAMEA1" type="biosample"/></SAMPLE></RECEIPT>')

    monkeypatch.setattr(ec, "http_post", fake_post)
    asm = tmp_path / "contigs.fa"
    asm.write_text("".join(f.read_text() for f in sorted(mag_set["fasta"].iterdir())))
    base = ["coassembly", "--xml-dir", str(xml_dir), "--samples", str(mag_set["samples"]), "--study", "PRJEB00001",
            "--assembly-fasta", str(asm), "--assembly-name", "g1_coassembly", "--program", "MEGAHIT v1.2.9",
            "--coverm", str(mag_set["coverm"]), "--outdir", str(tmp_path / "prim"), "--workdir", str(tmp_path / "wd")]
    assert es.main(base + ["--live"]) == 1 and posted == []          # helpdesk gate
    assert es.main(base + ["--live", "--helpdesk-notified"]) == 1 and posted == []   # no TEST run yet
    assert es.main(base) == 0
    assert posted == ["https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit"]
    mf = tmp_path / "prim" / "manifests_test" / "g1_coassembly.manifest"
    assert "SAMPLE\tERS99999999" in mf.read_text()
    web = [x for x in calls(fakes) if x["tool"] == "webin"]
    assert len(web) == 2 and all("--test" in w["argv"] for w in web)
    st = json.loads((tmp_path / "wd" / "ena_submit_state.json").read_text())
    assert st["coassembly"]["test"]["PRJEB00001_g1_coassembly"]["analysis"].startswith("ERZ")
    assert es.main(base + ["--live", "--helpdesk-notified"]) == 0
    assert posted[-1] == "https://www.ebi.ac.uk/ena/submit/webin-v2/submit"
    assert fake_creds not in all_text(tmp_path / "wd")
