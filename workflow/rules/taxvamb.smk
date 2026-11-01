"""TaxVAMB: MMseqs2 contig taxonomy, converted with taxconverter, then VAMB taxvamb binning."""


rule mmseqs2_taxonomy:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        f"{OUT}/binning/{{group}}/vamb/taxonomy/mmseqs2_taxonomy.tsv",
    log:
        f"{LOGS}/binning/{{group}}/mmseqs2_taxonomy.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/mmseqs2_taxonomy.tsv"
    conda:
        env("mmseqs2")
    container:
        image("mmseqs2")
    threads: threads_for("mmseqs2_taxonomy")
    resources:
        **resources_for(
            "mmseqs2_taxonomy", gpu=USE_GPU and config["binning"]["vamb"]["mmseqs_gpu"]
        ),
    params:
        db=config["databases"]["mmseqs2_taxonomy"],
        gpu="--gpu 1" if (USE_GPU and config["binning"]["vamb"]["mmseqs_gpu"]) else "",
        extra=config["binning"]["vamb"]["mmseqs_extra"],
    shell:
        r"""
        if [ -z "{params.db}" ]; then
            echo "Set databases.mmseqs2_taxonomy in the config" >{log}
            exit 1
        fi
        work=$(mktemp -d "${{TMPDIR:-/tmp}}/mmseqs.XXXXXX")
        trap 'rm -rf "$work"' EXIT
        mmseqs createdb {input.contigs} "$work/contigDB" >{log} 2>&1
        mmseqs taxonomy "$work/contigDB" {params.db} "$work/result" "$work/tmp" \
            --threads {threads} {params.gpu} {params.extra} >>{log} 2>&1
        mmseqs createtsv "$work/contigDB" "$work/result" {output} >>{log} 2>&1
        """


rule taxconverter:
    input:
        f"{OUT}/binning/{{group}}/vamb/taxonomy/mmseqs2_taxonomy.tsv",
    output:
        f"{OUT}/binning/{{group}}/vamb/taxonomy/taxconverter.tsv",
    log:
        f"{LOGS}/binning/{{group}}/taxconverter.log",
    conda:
        env("taxconverter")
    container:
        image("taxconverter")
    threads: threads_for("taxconverter")
    resources:
        **resources_for("taxconverter"),
    shell:
        "taxconverter mmseqs2 -i {input} -o {output} > {log} 2>&1"


rule vamb_taxonomy:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        taxonomy=f"{OUT}/binning/{{group}}/vamb/taxonomy/taxconverter.tsv",
        mmseqs=f"{OUT}/binning/{{group}}/vamb/taxonomy/mmseqs2_taxonomy.tsv",
    output:
        f"{OUT}/binning/{{group}}/vamb/taxonomy/taxonomy_vamb.tsv",
    log:
        f"{LOGS}/binning/{{group}}/vamb_taxonomy.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    shell:
        "python $BFM_SCRIPTS/fix_vamb_taxonomy.py --contigs {input.contigs} --taxonomy {input.taxonomy} "
        "--mmseqs {input.mmseqs} --out {output} > {log} 2>&1"


rule vamb:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        taxonomy=f"{OUT}/binning/{{group}}/vamb/taxonomy/taxonomy_vamb.tsv",
        bams=lambda wc: group_bams(wc.group),
    output:
        directory(f"{OUT}/binning/{{group}}/vamb/run"),
    log:
        f"{LOGS}/binning/{{group}}/vamb.log",
    benchmark:
        f"{BENCH}/binning/{{group}}/vamb.tsv"
    conda:
        env("vamb-gpu" if USE_GPU else "vamb")
    container:
        image("vamb-gpu" if USE_GPU else "vamb")
    threads: threads_for("vamb")
    resources:
        **resources_for("vamb", gpu=USE_GPU),
    params:
        bamdir=lambda wc: f"{OUT}/mapping/{wc.group}/bam",
        cuda="--cuda" if USE_GPU else "",
        p=config["binning"]["vamb"],
    shell:
        r"""
        rm -rf {output} # VAMB refuses an existing output directory
        vamb bin taxvamb --outdir {output} --fasta {input.contigs} --bamdir {params.bamdir} \
            --taxonomy {input.taxonomy} -p {threads} {params.cuda} -m {params.p[min_contig]} \
            -o "" --minfasta {params.p[min_fasta]} {params.p[extra]} >{log} 2>&1
        """
