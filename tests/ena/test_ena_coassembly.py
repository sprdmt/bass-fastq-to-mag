import json
import logging
import xml.etree.ElementTree as ET

import pytest

import ena_coassembly as ec
from conftest import RUNS, SAMPLE_ACCS
from ena_common import read_tsv

OK_RECEIPT = ('<?xml version="1.0" encoding="UTF-8"?><RECEIPT receiptDate="2026-10-31T12:00:00Z" success="true">'
              '<SAMPLE accession="ERS99999999" alias="{alias}" status="PRIVATE"><EXT_ID accession="SAMEA99999999" '
              'type="biosample"/></SAMPLE><SUBMISSION accession="ERA1" alias="sub"/><MESSAGES><INFO>This submission is '
              'a TEST submission and will be discarded within 24 hours</INFO></MESSAGES><ACTIONS>ADD</ACTIONS></RECEIPT>')
EXISTS_RECEIPT = ('<RECEIPT success="false"><MESSAGES><ERROR>In sample, alias: "{alias}". The object being added already '
                  'exists in the submission account with accession: "ERS11111111".</ERROR></MESSAGES></RECEIPT>')


@pytest.fixture
def xml_dir(mag_set, tmp_path):
    params = tmp_path / "p.json"
    params.write_text(json.dumps({"ena": {"study": "PRJEB00001", "metagenome": "seawater metagenome",
                                          "metagenome_taxid": 1561972, "centre_name": "University of Vienna",
                                          "coassembly": {"attributes": {"collection date": "2022-09",
                                                         "geographic location (country and/or sea)": "North Sea"}}}}))
    out = tmp_path / "co"
    assert ec.main(["sample-xml", "--samples", str(mag_set["samples"]), "--group", "g1", "--params", str(params),
                    "--outdir", str(out)]) == 0
    return out


def attrs_of(path):
    root = ET.parse(path).getroot()
    return {a.findtext("TAG"): a.findtext("VALUE") for a in root.iter("SAMPLE_ATTRIBUTE")}, root


def test_sample_xml_content(xml_dir):
    attrs, root = attrs_of(xml_dir / "coassembly_sample.xml")
    assert attrs["sample derived from"] == ",".join(SAMPLE_ACCS)          # comma list, no spaces
    assert attrs["ENA-CHECKLIST"] == "ERC000011"
    assert attrs["geographic location (country and/or sea)"] == "North Sea"
    assert root.find(".//TAXON_ID").text == "1561972"
    assert root.find(".//SCIENTIFIC_NAME").text == "seawater metagenome"
    assert root.find("SAMPLE").get("alias") == "PRJEB00001_g1_coassembly"
    webin = ET.parse(xml_dir / "webin_submission.xml").getroot()
    assert webin.tag == "WEBIN" and webin.find("SUBMISSION_SET/SUBMISSION/ACTIONS/ACTION/ADD") is not None
    info = json.loads((xml_dir / "coassembly_sample.json").read_text())
    assert info["runs"] == RUNS and info["derived_from_kind"] == "sample"


def test_missing_mandatory_values_use_sample_group_term(mag_set, tmp_path, caplog):
    out = tmp_path / "co2"
    with caplog.at_level(logging.WARNING):
        assert ec.main(["sample-xml", "--samples", str(mag_set["samples"]), "--taxon-id", "408172",
                        "--scientific-name", "marine metagenome", "--outdir", str(out)]) == 0
    attrs, _ = attrs_of(out / "coassembly_sample.xml")
    assert attrs["collection date"] == "missing: sample group"
    assert attrs["geographic location (country and/or sea)"] == "missing: sample group"


def test_run_accessions_used_when_no_sample_accessions(mag_set, tmp_path):
    s = mag_set["samples"]
    s.write_text("\n".join("\t".join(l.split("\t")[:5]) for l in s.read_text().splitlines()) + "\n")
    out = tmp_path / "co3"
    assert ec.main(["sample-xml", "--samples", str(s), "--taxon-id", "1561972", "--scientific-name",
                    "seawater metagenome", "--outdir", str(out)]) == 0
    attrs, _ = attrs_of(out / "coassembly_sample.xml")
    assert attrs["sample derived from"] == ",".join(RUNS)


def test_submit_requires_credentials(xml_dir, no_creds, monkeypatch):
    calls = []
    monkeypatch.setattr(ec, "http_post", lambda *a, **k: calls.append(a) or (200, ""))
    assert ec.main(["submit-sample", "--xml-dir", str(xml_dir)]) == 2
    assert calls == []


def test_submit_test_server_default_resume_and_secret(xml_dir, fake_creds, monkeypatch, caplog):
    calls = []

    def fake_post(url, body, headers, user, password, timeout=300):
        calls.append((url, body, headers, user))
        alias = ET.fromstring(body).find(".//SAMPLE").get("alias")
        return 200, OK_RECEIPT.format(alias=alias)

    monkeypatch.setattr(ec, "http_post", fake_post)
    with caplog.at_level(logging.INFO):
        assert ec.main(["submit-sample", "--xml-dir", str(xml_dir), "--test-suffix", "1700000000"]) == 0
    assert len(calls) == 1
    url, body, headers, user = calls[0]
    assert url == "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/submit"
    assert b"PRJEB00001_g1_coassembly_1700000000" in body and b"<WEBIN>" in body
    rows = read_tsv(xml_dir / "coassembly_sample_accessions.tsv")
    assert rows[0]["accession"] == "ERS99999999" and rows[0]["mode"] == "test" and rows[0]["biosample"] == "SAMEA99999999"
    # resume: nothing is posted again
    assert ec.main(["submit-sample", "--xml-dir", str(xml_dir)]) == 0
    assert len(calls) == 1
    # live goes to the production endpoint and uses the plain alias
    assert ec.main(["submit-sample", "--xml-dir", str(xml_dir), "--live"]) == 0
    assert calls[1][0] == "https://www.ebi.ac.uk/ena/submit/webin-v2/submit"
    assert b'alias="PRJEB00001_g1_coassembly"' in calls[1][1]
    for f in xml_dir.iterdir():
        assert fake_creds not in f.read_text(errors="ignore")
    assert all(fake_creds not in m for m in caplog.messages)


def test_submit_v1_multipart_and_existing_sample(xml_dir, fake_creds, monkeypatch):
    seen = {}

    def fake_post(url, body, headers, user, password, timeout=300):
        seen.update(url=url, body=body, ctype=headers["Content-Type"])
        return 200, EXISTS_RECEIPT.format(alias="PRJEB00001_g1_coassembly_42")

    monkeypatch.setattr(ec, "http_post", fake_post)
    assert ec.main(["submit-sample", "--xml-dir", str(xml_dir), "--api", "v1", "--test-suffix", "42"]) == 0
    assert seen["url"] == "https://wwwdev.ebi.ac.uk/ena/submit/drop-box/submit/"
    assert seen["ctype"].startswith("multipart/form-data") and b'name="SAMPLE"' in seen["body"] \
        and b'name="SUBMISSION"' in seen["body"]
    assert read_tsv(xml_dir / "coassembly_sample_accessions.tsv")[0]["accession"] == "ERS11111111"


def test_receipt_failure(xml_dir, fake_creds, monkeypatch):
    monkeypatch.setattr(ec, "http_post", lambda *a, **k: (400, '<RECEIPT success="false"><MESSAGES><ERROR>'
                                                                'Invalid attribute</ERROR></MESSAGES></RECEIPT>'))
    assert ec.main(["submit-sample", "--xml-dir", str(xml_dir)]) == 1
    assert not (xml_dir / "coassembly_sample_accessions.tsv").exists()


def test_manifest(mag_set, tmp_path):
    asm = tmp_path / "contigs.fa"
    asm.write_text("".join(f.read_text() for f in sorted(mag_set["fasta"].iterdir())))
    out = tmp_path / "mf"
    rc = ec.main(["manifest", "--samples", str(mag_set["samples"]), "--study", "PRJEB00001", "--mode", "test",
                  "--sample-accession", "ERS99999999", "--coverm", str(mag_set["coverm"]), "--assembly-fasta", str(asm),
                  "--assembly-name", "g1_coassembly", "--program", "MEGAHIT v1.2.9", "--outdir", str(out)])
    assert rc == 0
    text = (out / "manifests_test" / "g1_coassembly.manifest").read_text()
    fields = dict(l.split("\t", 1) for l in text.strip().splitlines())
    assert fields["SAMPLE"] == "ERS99999999" and fields["ASSEMBLY_TYPE"] == "primary metagenome"
    assert fields["COVERAGE"] == "60" and fields["RUN_REF"] == ",".join(RUNS)
    assert fields["FASTA"].endswith("g1_coassembly.fa.gz")
    # a comma-separated SAMPLE list (the failed earlier attempt) is rejected
    rc = ec.main(["manifest", "--samples", str(mag_set["samples"]), "--study", "PRJEB00001", "--sample-accession",
                  ",".join(SAMPLE_ACCS), "--coverage", "48.1", "--assembly-fasta", str(asm), "--assembly-name",
                  "g1_coassembly", "--program", "MEGAHIT v1.2.9", "--outdir", str(out)])
    assert rc == 1


def test_parse_receipt_and_endpoints():
    r = ec.parse_receipt(OK_RECEIPT.format(alias="x"))
    assert r["success"] and r["samples"]["x"] == ("ERS99999999", "SAMEA99999999")
    assert all("wwwdev" in u for (api, live), u in ec.ENDPOINTS.items() if not live)
    assert all("wwwdev" not in u for (api, live), u in ec.ENDPOINTS.items() if live)
