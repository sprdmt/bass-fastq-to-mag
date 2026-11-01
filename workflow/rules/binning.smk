"""Four binners on each co-assembly. Every binner ends in a contig-to-bin table
(binning/{group}/{binner}/contig2bin.tsv) that GraphBin2 and Binette consume."""


# ---------------------------------------------------------------------------
# MetaBAT2
# ---------------------------------------------------------------------------
rule metabat2:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        depth=f"{OUT}/mapping/{{group}}/jgi_depth.txt",
    output:
        directory(f"{OUT}/binning/{{group}}/metabat2/run"),
    log:
        f"{LOGS}/binning/{{group}}/metabat2.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/metabat2.tsv"
    conda:
        env("metabat2")
    container:
        image("metabat2")
    threads: threads_for("metabat2")
    resources:
        **resources_for("metabat2"),
    params:
        p=config["binning"]["metabat2"],
    shell:
        r"""
        mkdir -p {output}
        metabat2 -i {input.contigs} -a {input.depth} -o {output}/bin \
            -m {params.p[min_contig]} -s {params.p[min_bin_size]} --seed {params.p[seed]} \
            -t {threads} {params.p[extra]} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# GenomeFace: marker genes are searched in parallel chunks, then binned on GPU
# ---------------------------------------------------------------------------
GF_CHUNKS = config["binning"]["genomeface"]["marker_chunks"]


rule genomeface_split:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        temp(
            expand(
                f"{OUT}/binning/{{{{group}}}}/genomeface/chunks/chunk_{{chunk}}.fa",
                chunk=range(1, GF_CHUNKS + 1),
            )
        ),
    log:
        f"{LOGS}/binning/{{group}}/genomeface_split.log",
    conda:
        env("seqkit")
    container:
        image("seqkit")
    threads: threads_for("genomeface_split")
    resources:
        **resources_for("genomeface_split"),
    params:
        n=GF_CHUNKS,
        outdir=lambda wc: f"{OUT}/binning/{wc.group}/genomeface/chunks",
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/gfsplit.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        seqkit split2 --by-part {params.n} --out-dir "$tmp" {input.contigs} >{log} 2>&1
        mkdir -p {params.outdir}
        i=1
        for f in $(ls "$tmp" | sort); do
            mv "$tmp/$f" "{params.outdir}/chunk_$i.fa"
            i=$((i + 1))
        done
        n=$((i - 1))
        if [ "$n" -ne {params.n} ]; then
            echo "Expected {params.n} chunks, got $n" >>{log}
            exit 1
        fi
        """


rule genomeface_markers:
    input:
        chunk=f"{OUT}/binning/{{group}}/genomeface/chunks/chunk_{{chunk}}.fa",
    output:
        temp(f"{OUT}/binning/{{group}}/genomeface/markers/chunk_{{chunk}}.tsv"),
    log:
        f"{LOGS}/binning/{{group}}/genomeface_markers/chunk_{{chunk}}.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/genomeface_markers/chunk_{{chunk}}.tsv"
    conda:
        env("genomeface")
    container:
        image("genomeface")
    threads: threads_for("genomeface_markers")
    resources:
        **resources_for("genomeface_markers"),
    shell:
        r"""
        # markersgf writes scratch files to the working directory
        out=$(readlink -f {output})
        inp=$(readlink -f {input.chunk})
        lg=$(readlink -f {log})
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/markersgf.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        cd "$tmp"
        markersgf -t {threads} "$inp" markers.tsv >"$lg" 2>&1
        mv markers.tsv "$out"
        """


rule genomeface:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        depth=f"{OUT}/mapping/{{group}}/jgi_depth.txt",
        markers=expand(
            f"{OUT}/binning/{{{{group}}}}/genomeface/markers/chunk_{{chunk}}.tsv",
            chunk=range(1, GF_CHUNKS + 1),
        ),
    output:
        directory(f"{OUT}/binning/{{group}}/genomeface/run"),
    log:
        f"{LOGS}/binning/{{group}}/genomeface.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/genomeface.tsv"
    conda:
        env("genomeface")
    container:
        image("genomeface")
    threads: threads_for("genomeface")
    resources:
        **resources_for("genomeface", gpu=USE_GPU),
    params:
        p=config["binning"]["genomeface"],
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/genomeface.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        cat {input.markers} >"$tmp/markers.tsv"
        genomeface -i {input.contigs} -a {input.depth} -g "$tmp/markers.tsv" \
            -o {output} -m {params.p[min_contig]} -p {params.p[min_purity]} \
            {params.p[extra]} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# SemiBin2 (self-supervised, multi-sample coverage)
# ---------------------------------------------------------------------------
rule semibin2:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        bams=lambda wc: group_bams(wc.group),
    output:
        directory(f"{OUT}/binning/{{group}}/semibin2/run"),
    log:
        f"{LOGS}/binning/{{group}}/semibin2.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/semibin2.tsv"
    conda:
        env("semibin2-gpu" if USE_GPU else "semibin2")
    container:
        image("semibin2-gpu" if USE_GPU else "semibin2")
    threads: threads_for("semibin2")
    resources:
        **resources_for("semibin2", gpu=USE_GPU),
    params:
        engine="gpu" if USE_GPU else "cpu",
        p=config["binning"]["semibin2"],
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/semibin2.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        SemiBin2 single_easy_bin -i {input.contigs} -b {input.bams} -o {output} \
            --engine {params.engine} --threads {threads} --tmpdir "$tmp" \
            --min-len {params.p[min_len]} --compression none {params.p[extra]} >{log} 2>&1
        """


# ---------------------------------------------------------------------------
# Contig-to-bin tables (one per binner)
# ---------------------------------------------------------------------------
rule bins_table:
    input:
        lambda wc: binner_rundir(wc.group, wc.binner),
    output:
        f"{OUT}/binning/{{group}}/{{binner}}/contig2bin.tsv",
    log:
        f"{LOGS}/binning/{{group}}/{{binner}}_table.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        bins=lambda wc: binner_binsdir(wc.group, wc.binner),
    shell:
        "python $BFM_SCRIPTS/bins_to_table.py --bins-dir {params.bins} --prefix {wildcards.binner}_ "
        "--exclude bins.tsv --out {output} > {log} 2>&1"
