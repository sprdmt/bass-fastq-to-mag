"""Assembly-graph-based refinement of each binner's result with a corrected, fast GraphBin2.

workflow/scripts/graphbin2_fast.py re-implements GraphBin2 1.3.3 (MEGAHIT mode) with the
upstream algorithm, fixes for two upstream bugs and data-structure speed-ups (see
docs/graphbin2.md). It needs only the Python standard library and runs inside the stock
Bioconda GraphBin2 image.
"""


rule graphbin2:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        graph=f"{OUT}/assembly/{{group}}/assembly_graph.gfa",
        table=f"{OUT}/binning/{{group}}/{{binner}}/contig2bin.tsv",
    output:
        table=f"{OUT}/graphbin2/{{group}}/{{binner}}/contig2bin.tsv",
    log:
        f"{LOGS}/graphbin2/{{group}}/{{binner}}.log",
    benchmark:
        f"{BENCH}/graphbin2/{{group}}/{{binner}}.tsv"
    conda:
        env("graphbin2")
    container:
        image("graphbin2")
    threads: threads_for("graphbin2")
    resources:
        **resources_for("graphbin2"),
    params:
        depth=config["graphbin2"]["depth"],
        threshold=config["graphbin2"]["threshold"],
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/graphbin2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        if [ ! -s {input.table} ]; then
            echo "No initial bins: nothing to refine" >{log}
            : >{output.table}
            exit 0
        fi
        tr '\t' ',' <{input.table} >"$tmp/initial.csv"
        python $BFM_SCRIPTS/graphbin2_fast.py --contigs {input.contigs} --graph {input.graph} \
            --binned "$tmp/initial.csv" --output "$tmp/" --depth {params.depth} \
            --threshold {params.threshold} --delimiter , >{log} 2>&1
        cat "$tmp"/*graphbin2.log >>{log} 2>/dev/null || true
        # GraphBin2 writes CRLF line ends
        tr ',' '\t' <"$tmp/graphbin2_output.csv" | tr -d '\r' >{output.table}
        """
