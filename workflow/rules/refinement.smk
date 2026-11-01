"""Consensus bins with Binette (per assembly group), then candidate collection and
dereplication with dRep across all groups."""


rule binette:
    input:
        tables=lambda wc: [t for _, t in binette_tables(wc.group)],
        contigs=f"{OUT}/assembly/{{group}}/contigs.fa",
        proteins=f"{OUT}/genes/{{group}}/proteins.faa",
    output:
        bins=directory(f"{OUT}/binette/{{group}}/final_bins"),
        report=f"{OUT}/binette/{{group}}/final_bins_quality_reports.tsv",
    log:
        f"{LOGS}/binette/{{group}}.log",
    benchmark:
        f"{BENCH}/binette/{{group}}.tsv"
    conda:
        env("binette")
    container:
        image("binette")
    threads: threads_for("binette")
    resources:
        **resources_for("binette"),
    params:
        outdir=lambda wc: f"{OUT}/binette/{wc.group}",
        db=config["databases"]["checkm2"],
        min_completeness=config["binette"]["min_completeness"],
        extra=config["binette"]["extra"],
        pairs=lambda wc: " ".join(
            f"{label}={t}" for label, t in binette_tables(wc.group)
        ),
    shell:
        r"""
        if [ -z "{params.db}" ]; then
            echo "Set databases.checkm2 in the config" >{log}
            exit 1
        fi
        # one named copy per non-empty table: the names become Binette's bin set names (origin)
        sets=$(mktemp -d "${{TMPDIR:-/tmp}}/binette_sets.XXXXXX")
        trap 'rm -rf "$sets"' EXIT
        tables=""
        for pair in {params.pairs}; do
            label=${{pair%%=*}}
            t=${{pair#*=}}
            if [ -s "$t" ]; then
                cp "$t" "$sets/$label.tsv"
                tables="$tables $sets/$label.tsv"
            fi
        done
        if [ -z "$tables" ]; then
            echo "All binners returned empty results" >{log}
            exit 1
        fi
        binette --contig2bin_tables $tables --contigs {input.contigs} --proteins {input.proteins} \
            --outdir {params.outdir} --threads {threads} \
            --min_completeness {params.min_completeness} --checkm2_db {params.db} \
            {params.extra} >{log} 2>&1
        """


rule mag_candidates:
    input:
        unpack(candidate_inputs),
    output:
        genomes=directory(f"{OUT}/mags/candidates/genomes"),
        info=f"{OUT}/mags/candidates/genome_info.csv",
        table=f"{OUT}/mags/candidates/candidates.tsv",
    log:
        f"{LOGS}/mags/candidates.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        groups=GROUPS,
        fmt="simple" if config["mags"]["source"] == "anvio" else "binette",
    shell:
        "python $BFM_SCRIPTS/collect_candidates.py --groups {params.groups} --bins {input.bins} "
        "--reports {input.reports} --report-format {params.fmt} --out-dir {output.genomes} "
        "--genome-info {output.info} --table {output.table} > {log} 2>&1"


rule drep:
    input:
        genomes=f"{OUT}/mags/candidates/genomes",
        info=f"{OUT}/mags/candidates/genome_info.csv",
    output:
        genomes=directory(f"{OUT}/mags/drep/dereplicated_genomes"),
        winners=f"{OUT}/mags/drep/data_tables/Wdb.csv",
        clusters=f"{OUT}/mags/drep/data_tables/Cdb.csv",
    log:
        f"{LOGS}/mags/drep.log",
    benchmark:
        f"{BENCH}/mags/drep.tsv"
    conda:
        env("drep")
    container:
        image("drep")
    threads: threads_for("drep")
    resources:
        **resources_for("drep"),
    params:
        workdir=lambda wc, output: os.path.dirname(output.genomes),
        p=config["drep"],
    shell:
        r"""
        ls {input.genomes}/*.fa >"{params.workdir}.genomes.txt" 2>/dev/null || true
        n=$(wc -l <"{params.workdir}.genomes.txt")
        if [ "$n" -lt 2 ]; then
            # dRep needs at least two genomes; a single candidate is its own representative.
            mkdir -p {output.genomes} $(dirname {output.winners})
            echo "genome,score,cluster" >{output.winners}
            echo "genome,secondary_cluster,primary_cluster" >{output.clusters}
            for g in $(cat "{params.workdir}.genomes.txt"); do
                cp "$g" {output.genomes}/
                echo "$(basename $g),0,1_1" >>{output.winners}
                echo "$(basename $g),1_1,1" >>{output.clusters}
            done
            echo "Only $n candidate genome(s); dRep skipped" >{log}
            exit 0
        fi
        dRep dereplicate {params.workdir} -g "{params.workdir}.genomes.txt" \
            --genomeInfo {input.info} -comp {params.p[min_completeness]} \
            -con {params.p[max_contamination]} -sa {params.p[ani]} -p {threads} \
            {params.p[extra]} >{log} 2>&1
        """
