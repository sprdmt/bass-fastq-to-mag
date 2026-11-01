"""Optional anvi'o workspace for manual refinement of the Binette bins (one per assembly group).

`snakemake anvio` builds the workspace and triages the bins (anvio_triage.py): bins that pass
the MAG thresholds are kept as they are, nearly complete but contaminated bins go into the
collection 'Refine' (anvio.review: all = every bin), the rest is dropped. All Binette bins
stay available in the collection 'Binette'. After refining with anvi-refine, set
`mags.source: anvio`: the final candidates are the kept bins plus the bins of 'Refine'
(re-checked with CheckM2), dereplicated across groups. See docs/manual_refinement.md.
"""

# Shell snippet: exit status 0 if the profile database holds a collection 'Refine'.
_HAS_REFINE_SH = r"""
has_refine() {{
    python - "$1" 2>/dev/null <<'PY'
import sqlite3, sys
db = sqlite3.connect(sys.argv[1])
row = db.execute("select 1 from collections_info where collection_name = 'Refine'").fetchone()
sys.exit(0 if row else 1)
PY
}}
"""


rule anvio_inputs:
    input:
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        gff=f"{OUT}/genes/{{group}}/genes.gff",
        faa=f"{OUT}/genes/{{group}}/proteins.faa",
    output:
        fasta=f"{OUT}/anvio/{{group}}/contigs.fa",
        calls=f"{OUT}/anvio/{{group}}/external_gene_calls.tsv",
    log:
        f"{LOGS}/anvio/{{group}}/inputs.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        min_len=config["anvio"]["min_contig_len"],
    shell:
        r"""
        # anvi'o needs plain contig names that match the BAM files (first header token)
        cut -d ' ' -f 1 {input.contigs} >{output.fasta}
        python $BFM_SCRIPTS/anvio_gene_calls.py --gff {input.gff} --faa {input.faa} --out {output.calls} >{log} 2>&1
        """


rule anvio_contigs_db:
    input:
        fasta=f"{OUT}/anvio/{{group}}/contigs.fa",
        calls=f"{OUT}/anvio/{{group}}/external_gene_calls.tsv",
    output:
        db=f"{OUT}/anvio/{{group}}/contigs.db",
        done=f"{OUT}/anvio/{{group}}/contigs_db.done",
    log:
        f"{LOGS}/anvio/{{group}}/contigs_db.log",
    benchmark:
        f"{BENCH}/anvio/{{group}}/contigs_db.tsv"
    conda:
        env("anvio")
    container:
        image("anvio")
    threads: threads_for("anvio_contigs_db")
    resources:
        **resources_for("anvio_contigs_db"),
    params:
        name=lambda wc: anvio_name(f"{PROJECT}_{wc.group}"),
        scg=config["anvio"]["scg_taxonomy"],
        scg_dir=config["databases"]["anvio_scg"],
    shell:
        r"""
        anvi-gen-contigs-database -f {input.fasta} -o {output.db} -n {params.name} \
            --external-gene-calls {input.calls} --skip-predict-frame \
            --ignore-internal-stop-codons -T {threads} >{log} 2>&1
        anvi-run-hmms -c {output.db} -T {threads} --just-do-it >>{log} 2>&1
        if [ "{params.scg}" = "True" ]; then
            scg_dir=""
            [ -n "{params.scg_dir}" ] && scg_dir="--scgs-taxonomy-data-dir {params.scg_dir}"
            anvi-run-scg-taxonomy -c {output.db} -T {threads} $scg_dir >>{log} 2>&1
        fi
        touch {output.done}
        """


rule anvio_profile:
    input:
        bam=f"{OUT}/mapping/{{group}}/bam/{{sample}}.bam",
        db=f"{OUT}/anvio/{{group}}/contigs.db",
        done=f"{OUT}/anvio/{{group}}/contigs_db.done",
    output:
        f"{OUT}/anvio/{{group}}/profiles/{{sample}}/PROFILE.db",
    log:
        f"{LOGS}/anvio/{{group}}/profile_{{sample}}.log",
    benchmark:
        f"{BENCH}/anvio/{{group}}/profile_{{sample}}.tsv"
    conda:
        env("anvio")
    container:
        image("anvio")
    threads: threads_for("anvio_profile")
    resources:
        **resources_for("anvio_profile"),
    params:
        outdir=lambda wc: f"{OUT}/anvio/{wc.group}/profiles/{wc.sample}",
        name=lambda wc: anvio_name(wc.sample),
        min_len=config["anvio"]["min_contig_len"],
    shell:
        r"""
        anvi-profile -i {input.bam} -c {input.db} -o {params.outdir} -S {params.name} \
            --min-contig-length {params.min_len} -W -T {threads} >{log} 2>&1
        """


rule anvio_merge:
    input:
        db=f"{OUT}/anvio/{{group}}/contigs.db",
        profiles=lambda wc: expand(
            f"{OUT}/anvio/{wc.group}/profiles/{{sample}}/PROFILE.db",
            sample=group_samples(wc.group),
        ),
    output:
        f"{OUT}/anvio/{{group}}/merged/PROFILE.db",
    log:
        f"{LOGS}/anvio/{{group}}/merge.log",
    benchmark:
        f"{BENCH}/anvio/{{group}}/merge.tsv"
    conda:
        env("anvio")
    container:
        image("anvio")
    threads: threads_for("anvio_merge")
    resources:
        **resources_for("anvio_merge"),
    params:
        outdir=lambda wc: f"{OUT}/anvio/{wc.group}/merged",
    shell:
        r"""
        n=$(echo {input.profiles} | wc -w)
        if [ "$n" -gt 1 ]; then
            anvi-merge {input.profiles} -c {input.db} -o {params.outdir} -W >{log} 2>&1
        else
            # anvi-merge needs two profiles; a single profile is used directly
            rm -rf {params.outdir}
            cp -r $(dirname {input.profiles}) {params.outdir}
            echo "single sample: profile copied" >{log}
        fi
        """


rule anvio_triage:
    input:
        db=f"{OUT}/anvio/{{group}}/contigs.db",
        db_done=f"{OUT}/anvio/{{group}}/contigs_db.done",
        profile=f"{OUT}/anvio/{{group}}/merged/PROFILE.db",
        bins=f"{OUT}/binette/{{group}}/final_bins",
        report=f"{OUT}/binette/{{group}}/final_bins_quality_reports.tsv",
    output:
        triage=f"{OUT}/anvio/{{group}}/triage.tsv",
        all_table=f"{OUT}/anvio/{{group}}/collection_Binette.tsv",
        refine_table=f"{OUT}/anvio/{{group}}/collection_Refine.tsv",
        howto=f"{OUT}/anvio/{{group}}/HOW_TO_REFINE.md",
        done=f"{OUT}/anvio/{{group}}/collections.done",
    log:
        f"{LOGS}/anvio/{{group}}/triage.log",
    conda:
        env("anvio")
    container:
        image("anvio")
    threads: threads_for("anvio_collection")
    resources:
        **resources_for("anvio_collection"),
    params:
        keep_comp=config["mags"]["min_completeness"],
        keep_cont=config["mags"]["max_contamination"],
        refine_comp=config["anvio"]["refine_min_completeness"],
        refine_cont=config["anvio"]["refine_contamination_above"],
        review=config["anvio"]["review"],
        image=image("anvio") or "",
    shell:
        _HAS_REFINE_SH + r"""
        python $BFM_SCRIPTS/anvio_triage.py --bins-dir {input.bins} --report {input.report} \
            --keep-min-completeness {params.keep_comp} --keep-max-contamination {params.keep_cont} \
            --refine-min-completeness {params.refine_comp} \
            --refine-contamination-above {params.refine_cont} --review {params.review} \
            --group {wildcards.group} --profile {input.profile} --contigs-db {input.db} \
            --image "{params.image}" --triage {output.triage} --all-table {output.all_table} \
            --refine-table {output.refine_table} --howto {output.howto} >{log} 2>&1
        anvi-import-collection {output.all_table} -p {input.profile} -c {input.db} -C Binette \
            --contigs-mode >>{log} 2>&1
        # never replace a 'Refine' collection that may already hold manual work
        if has_refine {input.profile}; then
            echo "Collection 'Refine' exists already in {input.profile}; not replaced" >>{log}
        elif [ -s {output.refine_table} ]; then
            anvi-import-collection {output.refine_table} -p {input.profile} -c {input.db} \
                -C Refine --contigs-mode >>{log} 2>&1
        fi
        touch {output.done}
        """


rule anvio_refined_bins:
    input:
        db=f"{OUT}/anvio/{{group}}/contigs.db",
        profile=f"{OUT}/anvio/{{group}}/merged/PROFILE.db",
        done=f"{OUT}/anvio/{{group}}/collections.done",
    output:
        bins=directory(f"{OUT}/anvio/{{group}}/refined_bins"),
        export=f"{OUT}/anvio/{{group}}/refined_collection.txt",
    log:
        f"{LOGS}/anvio/{{group}}/summarize.log",
    conda:
        env("anvio")
    container:
        image("anvio")
    threads: threads_for("anvio_collection")
    resources:
        **resources_for("anvio_collection"),
    params:
        summary=lambda wc: f"{OUT}/anvio/{wc.group}/summary_Refine",
        prefix=lambda wc: f"{OUT}/anvio/{wc.group}/refined_collection",
    shell:
        _HAS_REFINE_SH + r"""
        mkdir -p {output.bins}
        rm -rf {params.summary}
        if ! has_refine {input.profile}; then
            echo "No collection 'Refine' (nothing to refine)" >{log}
            : >{output.export}
            exit 0
        fi
        # text copy of the refined collection (backup of the manual work)
        anvi-export-collection -p {input.profile} -C Refine -O {params.prefix} >{log} 2>&1
        anvi-summarize -c {input.db} -p {input.profile} -C Refine -o {params.summary} >>{log} 2>&1
        for f in {params.summary}/bin_by_bin/*/*-contigs.fa; do
            [ -e "$f" ] || continue
            b=$(basename "$f" -contigs.fa)
            cp "$f" {output.bins}/"$b".fa
        done
        """


rule checkm2_refined:
    input:
        bins=f"{OUT}/anvio/{{group}}/refined_bins",
    output:
        f"{OUT}/anvio/{{group}}/refined_checkm2/quality_report.tsv",
    log:
        f"{LOGS}/anvio/{{group}}/checkm2.log",
    benchmark:
        f"{BENCH}/anvio/{{group}}/checkm2.tsv"
    conda:
        env("checkm2")
    container:
        image("checkm2")
    threads: threads_for("checkm2")
    resources:
        **resources_for("checkm2"),
    params:
        outdir=lambda wc: f"{OUT}/anvio/{wc.group}/refined_checkm2",
        db=config["databases"]["checkm2"],
    shell:
        r"""
        if ! ls {input.bins}/*.fa >/dev/null 2>&1; then
            mkdir -p {params.outdir}
            printf "Name\tCompleteness\tContamination\n" >{output}
            echo "No refined bins" >{log}
            exit 0
        fi
        checkm2 predict --input {input.bins} -x fa --output-directory {params.outdir} \
            --threads {threads} --database_path {params.db} --tmpdir "${{TMPDIR:-/tmp}}" \
            --force >{log} 2>&1
        """


rule anvio_candidates:
    input:
        triage=f"{OUT}/anvio/{{group}}/triage.tsv",
        binette=f"{OUT}/binette/{{group}}/final_bins",
        refined=f"{OUT}/anvio/{{group}}/refined_bins",
        checkm2=f"{OUT}/anvio/{{group}}/refined_checkm2/quality_report.tsv",
    output:
        bins=directory(f"{OUT}/anvio/{{group}}/candidates"),
        quality=f"{OUT}/anvio/{{group}}/candidates_quality.tsv",
    log:
        f"{LOGS}/anvio/{{group}}/candidates.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    shell:
        "python $BFM_SCRIPTS/anvio_candidates.py --triage {input.triage} --binette-bins {input.binette} "
        "--refined-bins {input.refined} --refined-checkm2 {input.checkm2} "
        "--out-dir {output.bins} --quality {output.quality} > {log} 2>&1"
