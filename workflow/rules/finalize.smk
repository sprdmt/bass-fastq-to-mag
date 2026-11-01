"""Final MAG set: naming, taxonomy (GTDB-Tk), rRNA/tRNA genes for MIMAG, coverage and summary."""


rule finalize_mags:
    input:
        winners=f"{OUT}/mags/drep/data_tables/Wdb.csv",
        clusters=f"{OUT}/mags/drep/data_tables/Cdb.csv",
        genomes=f"{OUT}/mags/drep/dereplicated_genomes",
        candidates=f"{OUT}/mags/candidates/candidates.tsv",
    output:
        fasta=directory(f"{OUT}/mags/fasta"),
        table=f"{OUT}/mags/mags.tsv",
    log:
        f"{LOGS}/mags/finalize.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        project=PROJECT,
        min_completeness=config["mags"]["min_completeness"],
        max_contamination=config["mags"]["max_contamination"],
    shell:
        "python $BFM_SCRIPTS/finalize_mags.py --winners {input.winners} --clusters {input.clusters} "
        "--candidates {input.candidates} --genomes-dir {input.genomes} --project {params.project} "
        "--min-completeness {params.min_completeness} --max-contamination {params.max_contamination} "
        "--out-dir {output.fasta} --table {output.table} > {log} 2>&1"


rule gtdbtk:
    input:
        mags=f"{OUT}/mags/fasta",
    output:
        summary=f"{OUT}/mags/gtdbtk/gtdbtk.summary.tsv",
    log:
        f"{LOGS}/mags/gtdbtk.log",
    benchmark:
        f"{BENCH}/mags/gtdbtk.tsv"
    conda:
        env("gtdbtk")
    container:
        image("gtdbtk")
    threads: threads_for("gtdbtk")
    resources:
        **resources_for("gtdbtk"),
    params:
        outdir=lambda wc, output: os.path.dirname(output.summary),
        db=config["databases"]["gtdbtk"],
        extra=config["gtdbtk"]["extra"],
    shell:
        r"""
        if [ -z "{params.db}" ]; then
            echo "Set databases.gtdbtk in the config" >{log}
            exit 1
        fi
        export GTDBTK_DATA_PATH={params.db}
        if ! ls {input.mags}/*.fna >/dev/null 2>&1; then
            echo "No MAGs to classify" >{log}
            printf "user_genome\tclassification\n" >{output.summary}
            exit 0
        fi
        gtdbtk classify_wf --genome_dir {input.mags} --extension fna --out_dir {params.outdir} \
            --cpus {threads} --tmpdir "${{TMPDIR:-/tmp}}" {params.extra} >{log} 2>&1
        # one table for bacteria and archaea
        first=1
        for f in {params.outdir}/gtdbtk.bac120.summary.tsv {params.outdir}/gtdbtk.ar53.summary.tsv; do
            [ -s "$f" ] || continue
            if [ $first -eq 1 ]; then
                cat "$f"
                first=0
            else tail -n +2 "$f"; fi
        done >{output.summary}
        """


# Shell snippet: prints "arc" for MAGs classified as archaea by GTDB-Tk, else "bac".
_DOMAIN_SH = r"""
domain() {{
    if [ -n "{params.gtdbtk}" ] && awk -F'\t' -v m="$1" '$1==m && $2 ~ /^d__Archaea/ {{f=1}} END{{exit !f}}' {params.gtdbtk}; then
        echo arc; else echo bac; fi
}}
"""


rule barrnap_mags:
    input:
        mags=f"{OUT}/mags/fasta",
        gtdbtk=mag_domains_input,
    output:
        directory(f"{OUT}/mags/rna/barrnap"),
    log:
        f"{LOGS}/mags/barrnap.log",
    benchmark:
        f"{BENCH}/mags/barrnap.tsv"
    conda:
        env("barrnap")
    container:
        image("barrnap")
    threads: threads_for("mag_rna")
    resources:
        **resources_for("mag_rna"),
    params:
        gtdbtk=lambda wc, input: input.gtdbtk if input.gtdbtk else "",
    shell:
        _DOMAIN_SH + r"""
        mkdir -p {output}
        : > {log}
        for f in {input.mags}/*.fna; do
            [ -e "$f" ] || continue
            m=$(basename "$f" .fna); k=$(domain "$m")
            barrnap --kingdom "$k" --threads {threads} --quiet "$f" > {output}/"$m".gff 2>> {log}
        done
        """


rule trnascan_mags:
    input:
        mags=f"{OUT}/mags/fasta",
        gtdbtk=mag_domains_input,
    output:
        directory(f"{OUT}/mags/rna/trnascan"),
    log:
        f"{LOGS}/mags/trnascan.log",
    benchmark:
        f"{BENCH}/mags/trnascan.tsv"
    conda:
        env("trnascan")
    container:
        image("trnascan")
    threads: threads_for("mag_rna")
    resources:
        **resources_for("mag_rna"),
    params:
        gtdbtk=lambda wc, input: input.gtdbtk if input.gtdbtk else "",
    shell:
        _DOMAIN_SH + r"""
        mkdir -p {output}
        : > {log}
        for f in {input.mags}/*.fna; do
            [ -e "$f" ] || continue
            m=$(basename "$f" .fna)
            if [ "$(domain "$m")" = arc ]; then mode=-A; else mode=-B; fi
            tRNAscan-SE $mode --quiet --thread {threads} -o {output}/"$m".txt "$f" >> {log} 2>&1
        done
        """


rule mimag_rna:
    input:
        mags=f"{OUT}/mags/mags.tsv",
        barrnap=f"{OUT}/mags/rna/barrnap",
        trnascan=f"{OUT}/mags/rna/trnascan",
        gtdbtk=mag_domains_input,
    output:
        f"{OUT}/mags/mimag_rna.tsv",
    log:
        f"{LOGS}/mags/mimag_rna.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        gtdbtk=lambda wc, input: f"--gtdbtk {input.gtdbtk}" if input.gtdbtk else "",
        rna=mimag_rna_args(),
    shell:
        "python $BFM_SCRIPTS/mimag_rna.py --barrnap-dir {input.barrnap} --trnascan-dir {input.trnascan} "
        "--mags-table {input.mags} {params.gtdbtk} {params.rna} --out {output} > {log} 2>&1"


rule mag_coverage:
    input:
        mags=f"{OUT}/mags/mags.tsv",
        fasta=f"{OUT}/mags/fasta",
        coverm=expand(f"{OUT}/mapping/{{group}}/coverm_contig_mean.tsv", group=GROUPS),
    output:
        f"{OUT}/mags/coverage.tsv",
    log:
        f"{LOGS}/mags/coverage.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        coverm=" ".join(f"{g}={OUT}/mapping/{g}/coverm_contig_mean.tsv" for g in GROUPS),
        method=config["mags"]["coverage_method"],
    shell:
        "python $BFM_SCRIPTS/mag_coverage.py --mags-table {input.mags} --mags-dir {input.fasta} "
        "--coverm {params.coverm} --method {params.method} --out {output} > {log} 2>&1"


rule mag_summary:
    input:
        mags=f"{OUT}/mags/mags.tsv",
        gtdbtk=mag_domains_input,
        mimag=f"{OUT}/mags/mimag_rna.tsv",
        coverage=f"{OUT}/mags/coverage.tsv",
    output:
        f"{OUT}/mags/mag_summary.tsv",
    log:
        f"{LOGS}/mags/summary.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        gtdbtk=lambda wc, input: input.gtdbtk if input.gtdbtk else "",
    shell:
        "python $BFM_SCRIPTS/mag_summary.py --mags {input.mags} --gtdbtk '{params.gtdbtk}' "
        "--mimag {input.mimag} --coverage {input.coverage} --out {output} > {log} 2>&1"
