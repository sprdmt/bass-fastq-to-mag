"""Small synthetic fixtures for the ENA scripts (no network, no real data)."""
from __future__ import annotations

import random
import sys
import zipfile
from pathlib import Path

import pytest

SCRIPTS = Path(__file__).resolve().parents[2] / "workflow" / "scripts"
sys.path.insert(0, str(SCRIPTS))

RUNS = ["ERR0000003", "ERR0000002", "ERR0000001"]
SAMPLE_ACCS = ["ERS0000003", "ERS0000002", "ERS0000001"]


def rand_seq(n, seed):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


def write_fasta(path: Path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name} some description\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")


GFF_FULL = """##gff-version 3
k127_1\tbarrnap:0.9\trRNA\t10\t1500\t0\t+\t.\tName=16S_rRNA;product=16S ribosomal RNA
k127_1\tbarrnap:0.9\trRNA\t1600\t4500\t0\t+\t.\tName=23S_rRNA;product=23S ribosomal RNA
k127_2\tbarrnap:0.9\trRNA\t388\t493\t3.4e-11\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
"""
GFF_PARTIAL16 = """##gff-version 3
k127_3\tbarrnap:0.9\trRNA\t1\t700\t1e-50\t-\t.\tName=16S_rRNA;product=16S ribosomal RNA (partial);note=aligned only 45 percent of the 16S ribosomal RNA
k127_3\tbarrnap:0.9\trRNA\t800\t3800\t0\t-\t.\tName=23S_rRNA;product=23S ribosomal RNA
k127_4\tbarrnap:0.9\trRNA\t100\t210\t1e-9\t+\t.\tName=5S_rRNA;product=5S ribosomal RNA
"""
GFF_EMPTY = "##gff-version 3\n"

TRNA_HEADER = ("Sequence      \t\ttRNA\tBounds\ttRNA\tAnti\tIntron Bounds\tInf\t      \n"
               "Name          \ttRNA #\tBegin\tEnd \tType\tCodon\tBegin\tEnd\tScore\tNote\n"
               "--------      \t------\t-----\t------\t----\t-----\t-----\t----\t------\t------\n")
AA20 = ["Ala", "Arg", "Asn", "Asp", "Cys", "Gln", "Glu", "Gly", "His", "Ile2", "Leu", "Lys", "fMet",
        "Phe", "Pro", "Ser", "Thr", "Trp", "Tyr", "Val"]


def trna_table(types, pseudo=()):
    lines = [TRNA_HEADER]
    for i, t in enumerate(types, start=1):
        note = "pseudo" if t in pseudo else ""
        lines.append(f"k127_1 \t{i}\t{100 * i}\t{100 * i + 73}\t{t}\tNNN\t0\t0\t70.1\t{note}\n")
    return "".join(lines)


def write_xlsx(path: Path, sheets):
    """Minimal xlsx writer (shared strings) to test the stdlib reader without openpyxl."""
    strings, index = [], {}

    def sidx(s):
        if s not in index:
            index[s] = len(strings)
            strings.append(s)
        return index[s]

    sheet_xml = []
    for name, rows in sheets.items():
        xml_rows = []
        for r_i, row in enumerate(rows, start=1):
            cells = []
            for c_i, val in enumerate(row):
                ref = chr(65 + c_i) + str(r_i)
                if isinstance(val, (int, float)):
                    cells.append(f'<c r="{ref}"><v>{val}</v></c>')
                else:
                    cells.append(f'<c r="{ref}" t="s"><v>{sidx(str(val))}</v></c>')
            xml_rows.append(f'<row r="{r_i}">' + "".join(cells) + "</row>")
        sheet_xml.append((name, '<?xml version="1.0" encoding="UTF-8"?><worksheet xmlns="http://schemas.'
                          'openxmlformats.org/spreadsheetml/2006/main"><sheetData>' + "".join(xml_rows)
                          + "</sheetData></worksheet>"))
    ns = 'xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
    rns = 'xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships"'
    wb = (f'<?xml version="1.0" encoding="UTF-8"?><workbook {ns} {rns}><sheets>'
          + "".join(f'<sheet name="{n}" sheetId="{i}" r:id="rId{i}"/>' for i, (n, _x) in enumerate(sheet_xml, 1))
          + "</sheets></workbook>")
    rels = ('<?xml version="1.0" encoding="UTF-8"?><Relationships xmlns="http://schemas.openxmlformats.org/'
            'package/2006/relationships">' + "".join(
                f'<Relationship Id="rId{i}" Type="http://schemas.openxmlformats.org/officeDocument/2006/'
                f'relationships/worksheet" Target="worksheets/sheet{i}.xml"/>' for i in range(1, len(sheet_xml) + 1))
            + "</Relationships>")
    sst = (f'<?xml version="1.0" encoding="UTF-8"?><sst {ns} count="{len(strings)}" uniqueCount="{len(strings)}">'
           + "".join(f"<si><t>{s.replace('&', '&amp;').replace('<', '&lt;')}</t></si>" for s in strings) + "</sst>")
    with zipfile.ZipFile(path, "w") as z:
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", rels)
        z.writestr("xl/sharedStrings.xml", sst)
        for i, (_n, x) in enumerate(sheet_xml, 1):
            z.writestr(f"xl/worksheets/sheet{i}.xml", x)


# taxdump: root(1) > cellular organisms(131567) > Bacteria(2, domain) > Pseudomonadati(3379134, kingdom)
#   > Bacteroidota(976, phylum) > Flavobacteriia(117743, class) > Flavobacteriales(200644, order)
#   > Flavobacteriaceae(49546, family) > Polaribacter(52959, genus) > Polaribacter sp.(1909244, species)
# Archaea(2157) > Thermoproteota(28889) > Nitrososphaeria(1643678) > Nitrosopumilales(31932)
#   > Nitrosopumilaceae(338190) > Nitrosopumilus(338191)
# a homonym: Polaribacter (999999, genus) under Eukaryota(2759) to test domain-aware resolution
TAX_NODES = [
    (1, 1, "no rank"), (131567, 1, "cellular root"), (2, 131567, "domain"), (3379134, 2, "kingdom"),
    (976, 3379134, "phylum"), (117743, 976, "class"), (200644, 117743, "order"), (49546, 200644, "family"),
    (52959, 49546, "genus"), (1909244, 52959, "species"), (2157, 131567, "domain"), (28889, 2157, "phylum"),
    (1643678, 28889, "class"), (31932, 1643678, "order"), (338190, 31932, "family"), (338191, 338190, "genus"),
    (2759, 131567, "domain"), (999999, 2759, "genus"),
]
TAX_NAMES = [
    (1, "root", "scientific name"), (131567, "cellular organisms", "scientific name"),
    (2, "Bacteria", "scientific name"), (3379134, "Pseudomonadati", "scientific name"),
    (976, "Bacteroidota", "scientific name"), (976, "Bacteroidetes", "synonym"),
    (117743, "Flavobacteriia", "scientific name"), (200644, "Flavobacteriales", "scientific name"),
    (49546, "Flavobacteriaceae", "scientific name"), (52959, "Polaribacter", "scientific name"),
    (1909244, "Polaribacter sp.", "scientific name"), (2157, "Archaea", "scientific name"),
    (28889, "Thermoproteota", "scientific name"), (1643678, "Nitrososphaeria", "scientific name"),
    (31932, "Nitrosopumilales", "scientific name"), (338190, "Nitrosopumilaceae", "scientific name"),
    (338191, "Nitrosopumilus", "scientific name"), (2759, "Eukaryota", "scientific name"),
    (999999, "Polaribacter", "scientific name"),
]


def write_taxdump(directory: Path, as_zip: bool = False) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    nodes = "".join(f"{t}\t|\t{p}\t|\t{r}\t|\tXX\t|\t0\t|\n" for t, p, r in TAX_NODES)
    names = "".join(f"{t}\t|\t{n}\t|\t\t|\t{c}\t|\n" for t, n, c in TAX_NAMES)
    if as_zip:
        zp = directory / "taxdmp_test.zip"
        with zipfile.ZipFile(zp, "w") as z:
            z.writestr("nodes.dmp", nodes)
            z.writestr("names.dmp", names)
            z.writestr("fullnamelineage.dmp", "decoy\n")
        return zp
    (directory / "nodes.dmp").write_text(nodes)
    (directory / "names.dmp").write_text(names)
    return directory


@pytest.fixture
def mag_set(tmp_path):
    """Three MAGs of one co-assembly group with all upstream workflow outputs."""
    d = tmp_path / "results"
    fasta = d / "mags" / "fasta"
    barr = d / "mags" / "rna" / "barrnap"
    trna = d / "mags" / "rna" / "trnascan"
    for p in (fasta, barr, trna):
        p.mkdir(parents=True)
    contigs = {
        "g1_MAG_1": [("k127_1", rand_seq(5000, 1)), ("k127_2", rand_seq(1000, 2))],
        "g1_MAG_2": [("k127_3", rand_seq(4000, 3)), ("k127_4", rand_seq(2000, 4))],
        "g1_MAG_3": [("k127_5", rand_seq(3000, 5))],
    }
    for mag, recs in contigs.items():
        write_fasta(fasta / f"{mag}.fna", recs)
    (barr / "g1_MAG_1.gff").write_text(GFF_FULL)
    (barr / "g1_MAG_2.gff").write_text(GFF_PARTIAL16)
    (barr / "g1_MAG_3.gff").write_text(GFF_EMPTY)
    (trna / "g1_MAG_1.txt").write_text(trna_table(AA20 + ["SeC", "Undet", "Leu"]))
    (trna / "g1_MAG_2.txt").write_text(trna_table(AA20[:17] + ["Sup", "Pyl"]))
    (trna / "g1_MAG_3.txt").write_text(trna_table(AA20, pseudo=("Trp", "Tyr", "Val")))
    mags = d / "mags" / "mags.tsv"
    mags.write_text("mag_id\tcandidate\tgroup\tsource_bin\torigin\tcompleteness\tcontamination\tquality_from_checkm2\t"
                    "size_bp\tn_contigs\tn50\tdrep_cluster\tdrep_cluster_size\n"
                    "g1_MAG_1\tc1\tg1\tb1\tmetabat2\t95.05\t0.01\thigh\t6000\t2\t5000\t1_1\t1\n"
                    "g1_MAG_2\tc2\tg1\tb2\tsemibin2\t91.2\t4.99\thigh\t6000\t2\t4000\t2_1\t2\n"
                    "g1_MAG_3\tc3\tg1\tb3\tvamb\t62.4\t3.1\tmedium\t3000\t1\t3000\t3_1\t1\n")
    cov = d / "mapping" / "g1"
    cov.mkdir(parents=True)
    (cov / "coverm_contig_mean.tsv").write_text(
        "Contig\tS1 Mean\tS2 Mean\tS3 Mean\n"
        "k127_1\t10\t5\t1\n"
        "k127_2\t20\t0\t4\n"
        "k127_3\t2\t2\t2\n"
        "k127_4\t8\t2\t2\n"
        "k127_5\t1.5\t0.5\t0\n"
        "k127_99\t100\t100\t100\n")
    samples = tmp_path / "config" / "samples.tsv"
    samples.parent.mkdir()
    samples.write_text("sample\tfq1\tfq2\tassembly_group\trun_accession\tsample_accession\n" + "".join(
        f"S{i + 1}\tS{i + 1}_1.fq.gz\tS{i + 1}_2.fq.gz\tg1\t{RUNS[i]}\t{SAMPLE_ACCS[i]}\n" for i in range(3)))
    gtdb = d / "mags" / "gtdbtk" / "gtdbtk.summary.tsv"
    gtdb.parent.mkdir(parents=True)
    gtdb.write_text("user_genome\tclassification\tclosest_genome_reference\tnote\n"
                    "g1_MAG_1\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__Flavobacteriaceae;"
                    "g__Polaribacter;s__Polaribacter sp002163675\tGCA_1\t\n"
                    "g1_MAG_2\td__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Flavobacteriales;f__UBA1234;g__;s__\t\t\n"
                    "user_genome\tclassification\tclosest_genome_reference\tnote\n"
                    "g1_MAG_3\td__Archaea;p__Thermoproteota;c__Nitrososphaeria;o__Nitrosopumilales;"
                    "f__Nitrosopumilaceae;g__Nitrosopumilus;s__\t\t\n")
    return {"root": tmp_path, "results": d, "fasta": fasta, "barrnap": barr, "trnascan": trna, "mags": mags,
            "coverm": cov / "coverm_contig_mean.tsv", "samples": samples, "gtdbtk": gtdb, "contigs": contigs}


@pytest.fixture
def fake_creds(monkeypatch):
    monkeypatch.setenv("ENA_WEBIN", "Webin-00000")
    monkeypatch.setenv("ENA_WEBIN_PASSWORD", "S3cr3t-p@ss-word")
    return "S3cr3t-p@ss-word"


@pytest.fixture
def no_creds(monkeypatch):
    monkeypatch.delenv("ENA_WEBIN", raising=False)
    monkeypatch.delenv("ENA_WEBIN_PASSWORD", raising=False)
