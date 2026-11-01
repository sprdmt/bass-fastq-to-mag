"""Read QC: adapter/quality trimming with cutadapt and FastQC reports on raw and trimmed reads."""


rule cutadapt:
    input:
        r1=lambda wc: raw_reads(wc.sample)[0],
        r2=lambda wc: raw_reads(wc.sample)[1],
    output:
        r1=f"{OUT}/qc/trimmed/{{sample}}_R1.fastq.gz",
        r2=f"{OUT}/qc/trimmed/{{sample}}_R2.fastq.gz",
        json=f"{OUT}/qc/cutadapt/{{sample}}.cutadapt.json",
        report=f"{OUT}/qc/cutadapt/{{sample}}.cutadapt.txt",
    log:
        f"{LOGS}/qc/cutadapt/{{sample}}.log",
    benchmark:
        f"{BENCH}/qc/cutadapt/{{sample}}.tsv"
    conda:
        env("cutadapt")
    container:
        image("cutadapt")
    threads: threads_for("cutadapt")
    resources:
        **resources_for("cutadapt"),
    params:
        args=cutadapt_args(),
    shell:
        "cutadapt {params.args} -j {threads} --json {output.json} "
        "-o {output.r1} -p {output.r2} {input.r1} {input.r2} "
        "> {output.report} 2> {log}"


rule fastqc:
    input:
        r1=lambda wc: fastqc_reads(wc)[0],
        r2=lambda wc: fastqc_reads(wc)[1],
    output:
        zip1=f"{OUT}/qc/fastqc/{{stage}}/{{sample}}_R1_fastqc.zip",
        zip2=f"{OUT}/qc/fastqc/{{stage}}/{{sample}}_R2_fastqc.zip",
        html1=f"{OUT}/qc/fastqc/{{stage}}/{{sample}}_R1_fastqc.html",
        html2=f"{OUT}/qc/fastqc/{{stage}}/{{sample}}_R2_fastqc.html",
    log:
        f"{LOGS}/qc/fastqc/{{stage}}/{{sample}}.log",
    benchmark:
        f"{BENCH}/qc/fastqc/{{stage}}/{{sample}}.tsv"
    wildcard_constraints:
        stage="raw|trimmed",
    conda:
        env("fastqc")
    container:
        image("fastqc")
    threads: threads_for("fastqc")
    resources:
        **resources_for("fastqc"),
    shell:
        r"""
        # Link the reads under sample-based names so the report names are predictable.
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/fastqc.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        ext() {{ basename "$1" | sed -E 's/.*(\.(fastq|fq)(\.gz)?)$/\1/'; }}
        ln -s "$(readlink -f {input.r1})" "$tmp/{wildcards.sample}_R1$(ext {input.r1})"
        ln -s "$(readlink -f {input.r2})" "$tmp/{wildcards.sample}_R2$(ext {input.r2})"
        fastqc --quiet --threads {threads} --outdir "$(dirname {output.zip1})" "$tmp"/* >{log} 2>&1
        """


rule multiqc_reads:
    input:
        multiqc_inputs,
    output:
        f"{OUT}/reports/multiqc_reads.html",
    log:
        f"{LOGS}/qc/multiqc_reads.log",
    conda:
        env("multiqc")
    container:
        image("multiqc")
    threads: threads_for("multiqc")
    resources:
        **resources_for("multiqc"),
    shell:
        r"""
        multiqc --force --no-data-dir --filename $(basename {output}) \
            --outdir $(dirname {output}) {input} >{log} 2>&1
        """
