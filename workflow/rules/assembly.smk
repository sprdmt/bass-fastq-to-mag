"""Co-assembly with MEGAHIT and construction of the assembly graph used by GraphBin2."""


rule megahit:
    input:
        r1=lambda wc: group_reads(wc.group, 1),
        r2=lambda wc: group_reads(wc.group, 2),
    output:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        log=f"{OUT}/assembly/{{group}}/megahit.log",
    log:
        f"{LOGS}/assembly/{{group}}/megahit.log",
    benchmark:
        f"{BENCH}/assembly/{{group}}/megahit.tsv"
    conda:
        env("megahit")
    container:
        image("megahit")
    threads: threads_for("megahit")
    resources:
        **resources_for("megahit"),
    params:
        r1=lambda wc, input: ",".join(input.r1),
        r2=lambda wc, input: ",".join(input.r2),
        kmers=megahit_k_params(),
        min_len=config["assembly"]["min_contig_len"],
        min_count=config["assembly"]["min_count"],
        # MEGAHIT gets an absolute byte count from the job allocation. A fraction
        # (e.g. 0.95) would refer to the whole node and can exceed the allocation.
        memory=lambda wc, resources: mem_bytes(
            resources, config["assembly"]["memory_fraction"]
        ),
        workdir=lambda wc: config["assembly"]["workdir"] or f"{OUT}/assembly/{wc.group}",
        extra=config["assembly"]["extra"],
    shell:
        r"""
        work="{params.workdir}/megahit_{wildcards.group}.tmp"
        rm -rf "$work" && mkdir -p "$(dirname "$work")"
        megahit -1 {params.r1} -2 {params.r2} -o "$work" {params.kmers} \
            --min-contig-len {params.min_len} --min-count {params.min_count} \
            --num-cpu-threads {threads} --memory {params.memory} \
            --tmp-dir "${{TMPDIR:-/tmp}}" {params.extra} >{log} 2>&1
        cp "$work/final.contigs.fa" {output.contigs}
        cp "$work/log" {output.log}
        rm -rf "$work"
        """


rule contig2fastg:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        fastg=temp(f"{OUT}/assembly/{{group}}/assembly_graph.fastg"),
        k=f"{OUT}/assembly/{{group}}/final_k.txt",
    log:
        f"{LOGS}/assembly/{{group}}/contig2fastg.log",
    benchmark:
        f"{BENCH}/assembly/{{group}}/contig2fastg.tsv"
    conda:
        env("megahit")
    container:
        image("megahit")
    threads: threads_for("contig2fastg")
    resources:
        **resources_for("contig2fastg"),
    shell:
        r"""
        # The graph must be built with the k of the final MEGAHIT iteration, which is
        # encoded in the contig names (>k<K>_<id>). A wrong k yields a graph without edges.
        k=$(head -n 1 {input.contigs} | sed -nE 's/^>k([0-9]+)_.*/\1/p')
        if [ -z "$k" ]; then
            echo "Cannot read the final k-mer size from the MEGAHIT contig names" >{log}
            exit 1
        fi
        echo "$k" >{output.k}
        megahit_toolkit contig2fastg "$k" {input.contigs} >{output.fastg} 2>{log}
        """


rule fastg2gfa:
    input:
        fastg=f"{OUT}/assembly/{{group}}/assembly_graph.fastg",
        k=f"{OUT}/assembly/{{group}}/final_k.txt",
    output:
        gfa=f"{OUT}/assembly/{{group}}/assembly_graph.gfa",
        stats=f"{OUT}/assembly/{{group}}/assembly_graph_stats.tsv",
    log:
        f"{LOGS}/assembly/{{group}}/fastg2gfa.log",
    benchmark:
        f"{BENCH}/assembly/{{group}}/fastg2gfa.tsv"
    conda:
        env("agtools")
    container:
        image("agtools")
    threads: threads_for("fastg2gfa")
    resources:
        **resources_for("fastg2gfa"),
    params:
        require_edges=config["graph"]["require_edges"],
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/fastg2gfa.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        agtools fastg2gfa -g {input.fastg} -k "$(cat {input.k})" -o "$tmp" >{log} 2>&1
        mv "$tmp/converted_graph.gfa" {output.gfa}
        awk -F'\t' '$1=="S"{{s++}} $1=="L"{{l++}} END{{printf "segments\t%d\nlinks\t%d\n", s, l}}' \
            {output.gfa} >{output.stats}
        links=$(awk '$1=="links"{{print $2}}' {output.stats})
        echo "assembly graph: $(head -n 1 {output.stats} | cut -f2) segments, $links links" >>{log}
        if [ "$links" -eq 0 ] && [ "{params.require_edges}" = "True" ]; then
            echo "ERROR: the assembly graph has no edges, so graph-based refinement would" \
                "change nothing. Check the k-mer size in final_k.txt, or set" \
                "graph.require_edges: false to continue anyway." | tee -a {log} >&2
            exit 1
        fi
        """


rule quast:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        report=f"{OUT}/assembly/{{group}}/quast/report.tsv",
    log:
        f"{LOGS}/assembly/{{group}}/quast.log",
    benchmark:
        f"{BENCH}/assembly/{{group}}/quast.tsv"
    conda:
        env("quast")
    container:
        image("quast")
    threads: threads_for("quast")
    resources:
        **resources_for("quast"),
    params:
        outdir=lambda wc: f"{OUT}/assembly/{wc.group}/quast",
        min_len=min(500, config["assembly"]["min_contig_len"]),
    shell:
        r"""
        quast.py --fast --min-contig {params.min_len} --threads {threads} \
            --labels {wildcards.group} -o {params.outdir} {input.contigs} >{log} 2>&1
        """
