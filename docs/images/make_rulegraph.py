#!/usr/bin/env python3
"""Regenerate docs/images/rulegraph.{svg,png} (run from the repository root; needs graphviz).

The rule graph of the test configuration with the optional anvi'o and ENA branches switched on,
reduced transitively (tred) and coloured by workflow stage.
"""

import re
import subprocess
import tempfile
from pathlib import Path

STAGES = {
    "Read QC": ("#4C72B0", ["cutadapt", "fastqc", "multiqc_reads"]),
    "Assembly and graph": ("#55A868", ["megahit", "contig2fastg", "fastg2gfa", "quast", "pprodigal"]),
    "Coverage": ("#8172B2", ["bowtie2_index", "bowtie2_map", "coverm_contig", "jgi_depth"]),
    "Binning": ("#C44E52", ["metabat2", "genomeface_split", "genomeface_markers", "genomeface", "semibin2",
                            "mmseqs2_taxonomy", "taxconverter", "vamb_taxonomy", "vamb", "bins_table"]),
    "Graph refinement and consensus": ("#DD8452", ["graphbin2", "binette"]),
    "Manual refinement (optional)": ("#937860", ["anvio_inputs", "anvio_contigs_db", "anvio_profile",
                                                 "anvio_merge", "anvio_triage"]),
    "MAGs": ("#2A9D8F", ["mag_candidates", "drep", "finalize_mags", "gtdbtk", "barrnap_mags", "trnascan_mags",
                         "mimag_rna", "mag_coverage", "mag_summary"]),
    "ENA (optional)": ("#6E6E6E", ["ena_params", "gtdbtk_majority_vote", "ncbi_lineage", "ena_genome_table",
                                   "ena_coassembly_xml"]),
}


def main():
    with tempfile.TemporaryDirectory() as tmp:
        variant = Path(tmp) / "variant.yaml"
        variant.write_text("anvio: {enabled: true}\nena: {enabled: true}\n")
        out = subprocess.run(["snakemake", "--directory", ".test", "--configfile", str(variant), "--rulegraph"],
                             capture_output=True, text=True, check=True).stdout
        dot = out[out.index("digraph"):]
        nodes = dict(re.findall(r'^\s*(\d+)\[label = "([^"]+)"', dot, flags=re.M))
        edges = re.findall(r"^\s*(\d+) -> (\d+)", dot, flags=re.M)
        missing = set(nodes.values()) - {r for _, rules in STAGES.values() for r in rules} - {"all"}
        if missing:
            raise SystemExit(f"assign these rules to a stage: {sorted(missing)}")
        lines = ["digraph bfm {",
                 '  graph [rankdir=TB, bgcolor=white, nodesep=0.3, ranksep=0.32, fontname="Helvetica", margin=0.15];',
                 '  node [shape=box, style="rounded,filled", fontname="Helvetica", fontsize=11, penwidth=1.4,'
                 ' height=0.28, margin="0.08,0.03"];',
                 '  edge [color="#A0A0A0", penwidth=1.0, arrowsize=0.55];']
        for stage, (color, rules) in STAGES.items():
            style = "rounded,filled,dashed" if "optional" in stage else "rounded,filled"
            for rule in rules:
                if rule in nodes.values():
                    lines.append(f'  "{rule}" [color="{color}", fillcolor="{color}26", style="{style}"];')
        for a, b in edges:
            if "all" not in (nodes[a], nodes[b]):
                lines.append(f'  "{nodes[a]}" -> "{nodes[b]}";')
        rows = "".join(f'<TR><TD BGCOLOR="{c}40" BORDER="1" COLOR="{c}" WIDTH="14" HEIGHT="10"></TD>'
                       f'<TD ALIGN="LEFT"><FONT POINT-SIZE="11">{s}</FONT></TD></TR>' for s, (c, _) in STAGES.items())
        lines.append('  legend [shape=plaintext, style="", label=<<TABLE BORDER="0" CELLBORDER="0" CELLSPACING="3">'
                     f'<TR><TD COLSPAN="2" ALIGN="LEFT"><B>Stage</B></TD></TR>{rows}'
                     '<TR><TD COLSPAN="2" ALIGN="LEFT"><FONT POINT-SIZE="10">dashed border: optional</FONT></TD></TR>'
                     "</TABLE>>];")
        lines += ["  { rank=sink; legend; }", "}"]
        styled = Path(tmp) / "styled.dot"
        styled.write_text("\n".join(lines))
        reduced = subprocess.run(["tred", str(styled)], capture_output=True, text=True, check=True).stdout
        target = Path("docs/images")
        subprocess.run(["dot", "-Tsvg", "-o", str(target / "rulegraph.svg")], input=reduced, text=True, check=True)
        subprocess.run(["dot", "-Tpng", "-Gdpi=150", "-o", str(target / "rulegraph.png")], input=reduced,
                       text=True, check=True)
    print("docs/images/rulegraph.svg and .png written")


if __name__ == "__main__":
    main()
