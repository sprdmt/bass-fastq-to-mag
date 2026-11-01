import logging

import gtdb_to_ncbi as g2n
from conftest import write_taxdump, write_xlsx
from ena_common import read_tsv

MV_HEADER = "Genome ID\tGTDB classification\tMajority vote NCBI classification\n"


def test_lineage_helpers():
    n = g2n.split_lineage("d__Bacteria;p__Bacteroidota;c__;o__X;f__;g__Polaribacter;s__Polaribacter sp. A")
    assert n == ["Bacteria", "Bacteroidota", "", "X", "", "Polaribacter", "Polaribacter sp. A"]
    assert g2n.join_lineage(g2n.truncate(n, "genus")).endswith("g__Polaribacter;s__")
    assert g2n.lowest(g2n.truncate(n, "order")) == ("order", "X")


def test_majority_vote_postprocessing(mag_set, tmp_path):
    mv = tmp_path / "mv.tsv"
    mv.write_text(MV_HEADER +
                  "g1_MAG_1\tx\td__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;"
                  "f__Flavobacteriaceae;g__Polaribacter;s__Polaribacter sp. ALD11\n"
                  "g1_MAG_2\tx\td__;p__;c__;o__;f__;g__;s__\n")
    out = tmp_path / "lin.tsv"
    rc = g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(mv), "--out", str(out)])
    assert rc == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    assert rows["g1_MAG_1"]["ncbi_lineage"] == ("d__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;"
                                                "f__Flavobacteriaceae;g__Polaribacter;s__")
    assert rows["g1_MAG_1"]["ncbi_lineage_full"].endswith("s__Polaribacter sp. ALD11")
    assert rows["g1_MAG_2"]["ncbi_lineage"] == "d__Bacteria;p__;c__;o__;f__;g__;s__"
    assert "domain taken from GTDB" in rows["g1_MAG_2"]["note"]
    assert rows["g1_MAG_3"]["ncbi_lineage"].startswith("d__Archaea;p__;")
    assert "missing from majority-vote table" in rows["g1_MAG_3"]["note"]
    # species kept on request
    g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(mv), "--out", str(out), "--max-rank", "species"])
    assert {r["mag_id"]: r for r in read_tsv(out)}["g1_MAG_1"]["lowest_rank"] == "species"


def test_unclassified_genome(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("user_genome\tclassification\nm1\tUnclassified Bacteria\n")
    mv = tmp_path / "mv.tsv"
    mv.write_text(MV_HEADER + "m1\tUnclassified Bacteria\tUnclassified Bacteria\n")
    out = tmp_path / "o.tsv"
    assert g2n.main(["--gtdbtk", str(s), "--mv-table", str(mv), "--out", str(out)]) == 0
    r = read_tsv(out)[0]
    assert r["ncbi_lineage"] == "d__Bacteria;p__;c__;o__;f__;g__;s__" and "unclassified" in r["note"]


def test_no_domain_gives_exit_2(tmp_path):
    s = tmp_path / "s.tsv"
    s.write_text("user_genome\tclassification\nm1\tUnclassified\n")
    mv = tmp_path / "mv.tsv"
    mv.write_text(MV_HEADER + "m1\tUnclassified\tUnclassified\n")
    assert g2n.main(["--gtdbtk", str(s), "--mv-table", str(mv), "--out", str(tmp_path / "o.tsv")]) == 2


def test_xlsx_reader_and_mapping_table(mag_set, tmp_path, caplog):
    wb = tmp_path / "gtdb_vs_ncbi_test.xlsx"
    write_xlsx(wb, {
        "Genus": [["GTDB R226 genus", "Number of genomes", "Number of NCBI genera", "List of NCBI genera"],
                  ["g__Polaribacter", 40, 2, "g__Polaribacter 97.5%, g__(f__Flavobacteriaceae) 2.5%"]],
        "Family": [["GTDB R226 family", "Number of genomes", "Number of NCBI families", "List of NCBI families"],
                   ["f__UBA1234", 10, 2, "f__(o__Flavobacteriales) 80.0%, f__Flavobacteriaceae 20.0%"]],
        # a sheet whose rows do not match its header (as in the archaea workbook used before)
        "Order": [["GTDB R226 order", "Number of genomes", "Number of NCBI orders", "List of NCBI orders"],
                  ["f__Nitrosopumilaceae", 5, 1, "f__Nitrosopumilaceae 100%"]],
    })
    sheets = g2n.read_xlsx(wb)
    assert sheets["Genus"][1][0] == "g__Polaribacter" and sheets["Genus"][1][1] == "40"
    tax = write_taxdump(tmp_path / "taxdump", as_zip=True)
    out = tmp_path / "lin.tsv"
    with caplog.at_level(logging.WARNING):
        rc = g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--method", "mapping_table", "--mapping-table", str(wb),
                       "--taxdump", str(tax), "--out", str(out)])
    assert rc == 0
    assert any("do not carry the o__ prefix" in m for m in caplog.messages)
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    # mapped through the workbook; the eukaryotic homonym 'Polaribacter' is ignored
    assert rows["g1_MAG_1"]["ncbi_lineage"] == ("d__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;"
                                                "f__Flavobacteriaceae;g__Polaribacter;s__")
    assert rows["g1_MAG_1"]["ncbi_taxid"] == "52959" and "mapped g__Polaribacter" in rows["g1_MAG_1"]["note"]
    # no named NCBI family >= 50 % -> direct match of the GTDB order name
    assert rows["g1_MAG_2"]["ncbi_lineage"] == "d__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;f__;g__;s__"
    assert "direct match o__Flavobacteriales" in rows["g1_MAG_2"]["note"]
    assert rows["g1_MAG_3"]["ncbi_lineage"] == ("d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;"
                                                "f__Nitrosopumilaceae;g__Nitrosopumilus;s__")


def test_taxdump_synonym_and_name_check(tmp_path):
    tax = g2n.Taxdump(write_taxdump(tmp_path / "td"))
    hits, sci = tax.resolve({"Bacteroidetes", "Flavobacteriia"})
    assert hits["Bacteroidetes"] == [(976, "synonym")]
    assert tax.lineage_names(52959, tax.resolve({"Polaribacter"})[1])[5] == "Polaribacter"
    rows = [{"ncbi_lineage": "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__;f__;g__;s__", "note": ""}]
    g2n.check_names_in_taxdump(rows, tax)
    assert "c__Bacteroidia" in rows[0]["note"]


def test_parse_ncbi_list():
    assert g2n.parse_ncbi_list("c__Thermoplasmata 87.36%, c__(p__Euryarchaeota) 12.16%, c__Methanomicrobia 0.48%") == [
        ("Thermoplasmata", 87.36), ("Methanomicrobia", 0.48)]
    assert g2n.parse_ncbi_list("c__Candidatus Nanoarchaeia(p__Nanoarchaeota) 7.9%")[0][0] == "Candidatus Nanoarchaeia"


def _workbook(path):
    write_xlsx(path, {
        "Genus": [["GTDB R226 genus", "Number of genomes", "Number of NCBI genera", "List of NCBI genera"],
                  ["g__Polaribacter", 40, 1, "g__Polaribacter 100%"]],
        "Family": [["GTDB R226 family", "Number of genomes", "Number of NCBI families", "List of NCBI families"],
                   ["f__UBA1234", 10, 1, "f__(o__Flavobacteriales) 100%"]],
    })
    return path


def test_consensus_and_cross_check(mag_set, tmp_path, caplog):
    mv = tmp_path / "mv.tsv"
    mv.write_text(MV_HEADER +
                  # family conflict with the mapping-table route (Flavobacteriaceae)
                  "g1_MAG_1\tx\td__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;"
                  "f__Weeksellaceae;g__;s__\n"
                  # majority vote shallower than the mapping table, no conflict
                  "g1_MAG_2\tx\td__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__;f__;g__;s__\n"
                  # only the 'Candidatus' spelling differs -> agreement, current NCBI spelling kept
                  "g1_MAG_3\tx\td__Archaea;p__Candidatus Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;"
                  "f__;g__;s__\n")
    wb = _workbook(tmp_path / "wb.xlsx")
    tax = write_taxdump(tmp_path / "taxdump")
    base = ["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(mv), "--mapping-table", str(wb),
            "--taxdump", str(tax)]
    out = tmp_path / "cons.tsv"
    assert g2n.main(base + ["--method", "consensus", "--out", str(out)]) == 0
    rows = {r["mag_id"]: r for r in read_tsv(out)}
    assert rows["g1_MAG_1"]["ncbi_lineage"] == ("d__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__Flavobacteriales;"
                                                "f__;g__;s__")
    assert rows["g1_MAG_1"]["method"] == "consensus"
    assert "consensus stops above family" in rows["g1_MAG_1"]["note"]
    assert rows["g1_MAG_1"]["cross_check"] == "family: majority_vote Weeksellaceae vs mapping_table Flavobacteriaceae"
    assert rows["g1_MAG_2"]["ncbi_lineage"] == "d__Bacteria;p__Bacteroidota;c__Flavobacteriia;o__;f__;g__;s__"
    assert rows["g1_MAG_2"]["cross_check"] == "agree"
    assert rows["g1_MAG_3"]["ncbi_lineage"] == ("d__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;"
                                                "f__;g__;s__")
    assert rows["g1_MAG_3"]["cross_check"] == "agree"
    # cross-check keeps the majority-vote lineage and only reports the conflict
    out2 = tmp_path / "mv_cc.tsv"
    with caplog.at_level(logging.WARNING):
        assert g2n.main(base + ["--cross-check", "--out", str(out2)]) == 0
    assert any("1 genomes where majority_vote and mapping_table disagree" in m for m in caplog.messages)
    rows2 = {r["mag_id"]: r for r in read_tsv(out2)}
    assert rows2["g1_MAG_1"]["method"] == "majority_vote"
    assert "f__Weeksellaceae" in rows2["g1_MAG_1"]["ncbi_lineage"]
    assert rows2["g1_MAG_1"]["cross_check"].startswith("family:")
    # without the mapping inputs the combined modes refuse to run
    assert g2n.main(["--gtdbtk", str(mag_set["gtdbtk"]), "--mv-table", str(mv), "--method", "consensus",
                     "--out", str(tmp_path / "x.tsv")]) == 1
