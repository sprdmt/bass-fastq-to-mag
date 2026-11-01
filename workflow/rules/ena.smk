"""ENA submission files (docs/ena_submission.md).

`snakemake ena_prepare` builds the genome_uploader table with validated gzipped FASTAs, the NCBI
lineages and, for every assembly group with at least two samples, the co-assembly sample XML.
`snakemake ena_validate` checks the table offline and previews the ENA taxa (needs internet).
`snakemake ena_submit` submits to the ENA TEST server only (credentials from ENA_WEBIN and
ENA_WEBIN_PASSWORD). Live submission is deliberately not a workflow target: run
workflow/scripts/ena_submit.py with --live by hand after a successful test submission.
"""

import json
from pathlib import Path


localrules:
    ena_params,
    ena_validate,
    ena_submit_test,


rule ena_params:
    output:
        f"{OUT}/ena/ena_params.json",
    log:
        f"{LOGS}/ena/params.log",
    params:
        block=json.dumps({"ena": ena_effective_config()}, indent=1, sort_keys=True),
    run:
        require_ena_params()  # this job runs first: stop early on an incomplete ena: block
        Path(output[0]).write_text(f"{params.block}\n")
        Path(log[0]).write_text(f"{len(params.block)} characters written\n")


rule gtdbtk_majority_vote:
    input:
        summary=f"{OUT}/mags/gtdbtk/gtdbtk.summary.tsv",
        metadata=gtdb_metadata_inputs,
    output:
        f"{OUT}/ena/gtdbtk_majority_vote.tsv",
    log:
        f"{LOGS}/ena/gtdbtk_majority_vote.log",
    conda:
        env("gtdbtk")
    container:
        image("gtdbtk")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        gtdbtk_dir=lambda wc, input: os.path.dirname(input.summary),
        bac120=config["ena"]["gtdb_metadata"]["bac120"],
        ar53=config["ena"]["gtdb_metadata"]["ar53"],
    shell:
        r"""
        if [ -z "{params.bac120}" ] || [ -z "{params.ar53}" ]; then
            echo "Set ena.gtdb_metadata.bac120 and .ar53 (snakemake download_databases fetches them)" >{log}
            exit 1
        fi
        gtdb_to_ncbi_majority_vote.py --gtdbtk_output_dir {params.gtdbtk_dir} \
            --bac120_metadata_file {params.bac120} --ar53_metadata_file {params.ar53} \
            --output_file {output} >{log} 2>&1
        """


rule ncbi_lineage:
    input:
        unpack(ena_lineage_inputs),
    output:
        f"{OUT}/ena/ncbi_lineage.tsv",
    log:
        f"{LOGS}/ena/ncbi_lineage.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        args=ena_lineage_args,
    shell:
        "python $BFM_SCRIPTS/gtdb_to_ncbi.py --gtdbtk {input.summary} {params.args} --out {output} > {log} 2>&1"


rule ena_genome_table:
    input:
        mags=f"{OUT}/mags/mags.tsv",
        fasta=f"{OUT}/mags/fasta",
        rna=f"{OUT}/mags/mimag_rna.tsv",
        coverage=f"{OUT}/mags/coverage.tsv",
        lineage=f"{OUT}/ena/ncbi_lineage.tsv",
        samples=config["samples"],
        params=f"{OUT}/ena/ena_params.json",
    output:
        table=f"{OUT}/ena/genome_uploader_input.tsv",
        provenance=f"{OUT}/ena/genome_uploader_input.provenance.tsv",
        fasta=directory(f"{OUT}/ena/fasta"),
    log:
        f"{LOGS}/ena/genome_table.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    params:
        genome_type=config["ena"]["genome_type"],
    shell:
        "python $BFM_SCRIPTS/ena_genome_table.py --mags-table {input.mags} --mags-dir {input.fasta} "
        "--rna {input.rna} --coverage {input.coverage} --lineage {input.lineage} "
        "--samples {input.samples} --params {input.params} --genome-type {params.genome_type} "
        "--fasta-outdir {output.fasta} --out {output.table} "
        "--provenance-out {output.provenance} --force > {log} 2>&1"


rule ena_coassembly_xml:
    input:
        samples=config["samples"],
        params=f"{OUT}/ena/ena_params.json",
    output:
        directory(f"{OUT}/ena/coassembly/{{group}}/sample"),
    log:
        f"{LOGS}/ena/coassembly_xml_{{group}}.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    threads: threads_for("small")
    resources:
        **resources_for("small"),
    shell:
        "python $BFM_SCRIPTS/ena_coassembly.py sample-xml --samples {input.samples} --group {wildcards.group} "
        "--params {input.params} --outdir {output} > {log} 2>&1"


rule ena_validate:
    input:
        table=f"{OUT}/ena/genome_uploader_input.tsv",
    output:
        report=f"{OUT}/ena/validation/validate_table.txt",
        taxa=f"{OUT}/ena/validation/taxa_preview.tsv",
    log:
        f"{LOGS}/ena/validate.log",
    conda:
        env("ena")
    container:
        image("ena")
    params:
        genome_type=config["ena"]["genome_type"],
    shell:
        r"""
        python $BFM_SCRIPTS/ena_submit.py validate-table --genome-info {input.table} \
            --genome-type {params.genome_type} >{output.report} 2>&1
        # asks the ENA taxonomy service which taxon each genome would be registered as
        python $BFM_SCRIPTS/ena_submit.py preview-taxa --genome-info {input.table} --out {output.taxa} >{log} 2>&1
        """


rule ena_submit_test:
    input:
        unpack(ena_submit_inputs),
    output:
        f"{OUT}/ena/submit/test_submission.done",
    log:
        f"{LOGS}/ena/submit_test.log",
    conda:
        env("ena")
    container:
        image("ena")
    params:
        workdir=lambda wc, output: os.path.dirname(output[0]),
        genome_type=config["ena"]["genome_type"],
        study=config["ena"]["study"],
        centre=config["ena"]["centre_name"],
        private="--private" if config["ena"]["private"] else "--public",
        heap=config["ena"]["webin"]["java_heap_max_gb"],
        delay=config["ena"]["webin"]["delay_s"],
        coassembly=ena_coassembly_commands,
    shell:
        r"""
        if [ -z "${{ENA_WEBIN:-}}" ] || [ -z "${{ENA_WEBIN_PASSWORD:-}}" ]; then
            echo "Export ENA_WEBIN and ENA_WEBIN_PASSWORD before running ena_submit" | tee {log}
            exit 1
        fi
        if [ -z "{params.study}" ] || [ -z "{params.centre}" ]; then
            echo "Set ena.study and ena.centre_name in the config" | tee {log}
            exit 1
        fi
        # TEST server only; live submission: see docs/ena_submission.md
        python $BFM_SCRIPTS/ena_submit.py genomes --genome-info {input.table} --genome-type {params.genome_type} \
            --study {params.study} --centre-name "{params.centre}" {params.private} \
            --workdir {params.workdir} --java-heap-max-gb {params.heap} --delay {params.delay} >{log} 2>&1
        {params.coassembly}
        python $BFM_SCRIPTS/ena_submit.py status --workdir {params.workdir} >>{log} 2>&1
        touch {output}
        """
