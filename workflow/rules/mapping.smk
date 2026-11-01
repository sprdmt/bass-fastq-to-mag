"""Read mapping to the co-assembly and per-contig coverage/depth tables."""

BT2_EXT = (".1.bt2l", ".2.bt2l", ".3.bt2l", ".4.bt2l", ".rev.1.bt2l", ".rev.2.bt2l")


rule bowtie2_index:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        multiext(f"{OUT}/mapping/{{group}}/index/contigs", *BT2_EXT),
    log:
        f"{LOGS}/mapping/{{group}}/bowtie2_index.log",
    benchmark:
        f"{BENCH}/mapping/{{group}}/bowtie2_index.tsv"
    conda:
        env("mapping")
    container:
        image("mapping")
    threads: threads_for("bowtie2_index")
    resources:
        **resources_for("bowtie2_index"),
    params:
        prefix=lambda wc: f"{OUT}/mapping/{wc.group}/index/contigs",
    shell:
        "bowtie2-build --large-index --threads {threads} {input.contigs} {params.prefix} > {log} 2>&1"


rule bowtie2_map:
    input:
        r1=lambda wc: group_member_reads(wc)[0],
        r2=lambda wc: group_member_reads(wc)[1],
        index=multiext(f"{OUT}/mapping/{{group}}/index/contigs", *BT2_EXT),
    output:
        bam=f"{OUT}/mapping/{{group}}/bam/{{sample}}.bam",
        bai=f"{OUT}/mapping/{{group}}/bam/{{sample}}.bam.bai",
    log:
        f"{LOGS}/mapping/{{group}}/{{sample}}.bowtie2.log",
    benchmark:
        f"{BENCH}/mapping/{{group}}/{{sample}}.bowtie2.tsv"
    conda:
        env("mapping")
    container:
        image("mapping")
    threads: threads_for("bowtie2_map")
    resources:
        **resources_for("bowtie2_map"),
    params:
        prefix=lambda wc: f"{OUT}/mapping/{wc.group}/index/contigs",
        extra=config["mapping"]["bowtie2_extra"],
        sort_threads=lambda wc, threads: max(1, threads // 2),
        sort_mem=lambda wc, resources, threads: per_thread_mb(
            resources, max(1, threads // 2)
        ),
    shell:
        r"""
        set -o pipefail
        bowtie2 -x {params.prefix} -1 {input.r1} -2 {input.r2} --threads {threads} \
            {params.extra} 2>{log} \
            | samtools sort -@ {params.sort_threads} -m {params.sort_mem}M \
                -T "${{TMPDIR:-/tmp}}/{wildcards.group}.{wildcards.sample}.sort" \
                -o {output.bam} - 2>>{log}
        samtools index -@ {threads} {output.bam} {output.bai} 2>>{log}
        """


rule coverm_contig:
    input:
        bams=lambda wc: group_bams(wc.group),
    output:
        f"{OUT}/mapping/{{group}}/coverm_contig_mean.tsv",
    log:
        f"{LOGS}/mapping/{{group}}/coverm_contig.log",
    benchmark:
        f"{BENCH}/mapping/{{group}}/coverm_contig.tsv"
    conda:
        env("coverm")
    container:
        image("coverm")
    threads: threads_for("coverm_contig")
    resources:
        **resources_for("coverm_contig"),
    shell:
        "coverm contig --bam-files {input.bams} --methods mean --threads {threads} "
        "--output-file {output} > {log} 2>&1"


rule jgi_depth:
    input:
        bams=lambda wc: group_bams(wc.group),
    output:
        f"{OUT}/mapping/{{group}}/jgi_depth.txt",
    log:
        f"{LOGS}/mapping/{{group}}/jgi_depth.log",
    benchmark:
        f"{BENCH}/mapping/{{group}}/jgi_depth.tsv"
    conda:
        env("metabat2")
    container:
        image("metabat2")
    threads: threads_for("jgi_depth")
    resources:
        **resources_for("jgi_depth"),
    shell:
        r"""
        export OMP_NUM_THREADS={threads}
        jgi_summarize_bam_contig_depths --outputDepth {output} {input.bams} >{log} 2>&1
        """
