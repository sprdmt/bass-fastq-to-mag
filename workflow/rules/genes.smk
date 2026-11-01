"""Gene prediction on the co-assembly (proteins are reused by Binette and anvi'o)."""


rule pprodigal:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
    output:
        faa=f"{OUT}/genes/{{group}}/proteins.faa",
        fna=f"{OUT}/genes/{{group}}/genes.fna",
        gff=f"{OUT}/genes/{{group}}/genes.gff",
    log:
        f"{LOGS}/genes/{{group}}/pprodigal.log",
    benchmark:
        f"{BENCH}/genes/{{group}}/pprodigal.tsv"
    conda:
        env("pprodigal")
    container:
        image("pprodigal")
    threads: threads_for("pprodigal")
    resources:
        **resources_for("pprodigal"),
    params:
        procedure=config["genes"]["procedure"],
        chunksize=config["genes"]["chunksize"],
    shell:
        r"""
        pprodigal --tasks {threads} --chunksize {params.chunksize} -p {params.procedure} \
            -i {input.contigs} -a {output.faa} -d {output.fna} -f gff -o {output.gff} >{log} 2>&1
        """
