"""Optional download of the reference databases (`snakemake download_databases`).

Downloads run on the machine that runs Snakemake (they need internet access) into
`databases.download_dir`. Afterwards point `databases.checkm2`, `databases.gtdbtk`,
`databases.mmseqs2_taxonomy`, `databases.anvio_scg` and `ena.gtdb_metadata` to the printed
paths (database_paths.yaml).
"""

DB_DIR = config["databases"].get("download_dir", "databases")
MMSEQS_DB = config["databases"].get("mmseqs2_download", "GTDB")
# directory name without "/" (e.g. UniProtKB/Swiss-Prot -> UniProtKB_Swiss-Prot)
MMSEQS_DIR = re.sub(r"[^A-Za-z0-9._-]+", "_", MMSEQS_DB)


localrules:
    download_checkm2,
    download_gtdbtk,
    download_gtdb_metadata,
    download_mmseqs2,
    download_anvio_scg,
    download_databases,


rule download_checkm2:
    output:
        f"{DB_DIR}/checkm2/CheckM2_database/uniref100.KO.1.dmnd",
    log:
        f"{LOGS}/databases/checkm2.log",
    conda:
        env("checkm2")
    container:
        image("checkm2")
    params:
        outdir=lambda wc, output: os.path.dirname(os.path.dirname(output[0])),
    shell:
        "checkm2 database --download --path {params.outdir} > {log} 2>&1"


rule download_gtdbtk:
    output:
        touch(f"{DB_DIR}/gtdbtk/download.done"),
    log:
        f"{LOGS}/databases/gtdbtk.log",
    conda:
        env("gtdbtk")
    container:
        image("gtdbtk")
    params:
        url=config["databases"].get(
            "gtdbtk_url",
            "https://data.gtdb.ecogenomic.org/releases/release226/226.0/auxillary_files/"
            "gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz",
        ),
        outdir=lambda wc, output: os.path.dirname(output[0]),
    shell:
        r"""
        mkdir -p {params.outdir}
        wget -q -O - {params.url} | tar -xz -C {params.outdir} >{log} 2>&1
        echo "GTDB-Tk data: $(ls -d {params.outdir}/*/ | head -n 1)" >>{log}
        """


rule download_gtdb_metadata:
    output:
        bac120=f"{DB_DIR}/gtdb_metadata/bac120_metadata_r226.tsv.gz",
        ar53=f"{DB_DIR}/gtdb_metadata/ar53_metadata_r226.tsv.gz",
    log:
        f"{LOGS}/databases/gtdb_metadata.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    params:
        base=config["databases"]["gtdb_metadata_url"],
    shell:
        r"""
        for f in {output.bac120} {output.ar53}; do
            python -c "import sys, urllib.request; urllib.request.urlretrieve(sys.argv[1], sys.argv[2])" \
                {params.base}/$(basename $f) $f >>{log} 2>&1
        done
        """


rule download_mmseqs2:
    output:
        touch(f"{DB_DIR}/mmseqs2/{MMSEQS_DIR}/download.done"),
    log:
        f"{LOGS}/databases/mmseqs2.log",
    conda:
        env("mmseqs2")
    container:
        image("mmseqs2")
    threads: 16
    params:
        prefix=lambda wc, output: os.path.join(os.path.dirname(output[0]), "db"),
        name=MMSEQS_DB,
        gpu=USE_GPU and config["binning"]["vamb"]["mmseqs_gpu"],
    shell:
        r"""
        tmp=$(mktemp -d "${{TMPDIR:-/tmp}}/mmseqsdb.XXXXXX")
        trap 'rm -rf "$tmp"' EXIT
        mmseqs databases {params.name} {params.prefix} "$tmp" --threads {threads} >{log} 2>&1
        if [ "{params.gpu}" = "True" ]; then
            # GPU searches need a padded copy of the database; makepaddedseqdb does not carry
            # over the taxonomy files, so link them
            mmseqs makepaddedseqdb {params.prefix} {params.prefix}_padded >>{log} 2>&1
            for ext in _mapping _taxonomy _nodes.dmp _names.dmp _merged.dmp; do
                if [ -e "{params.prefix}$ext" ]; then
                    ln -sf "$(basename {params.prefix})$ext" "{params.prefix}_padded$ext"
                fi
            done
            echo "use databases.mmseqs2_taxonomy: {params.prefix}_padded" >>{log}
        fi
        """


rule download_anvio_scg:
    output:
        touch(f"{DB_DIR}/anvio_scg/download.done"),
    log:
        f"{LOGS}/databases/anvio_scg.log",
    conda:
        env("anvio")
    container:
        image("anvio")
    params:
        outdir=lambda wc, output: os.path.join(os.path.dirname(output[0]), "data"),
    shell:
        "anvi-setup-scg-taxonomy --scgs-taxonomy-data-dir {params.outdir} -T 4 > {log} 2>&1"


rule download_databases:
    input:
        checkm2=rules.download_checkm2.output[0],
        gtdbtk=rules.download_gtdbtk.output[0],
        bac120=rules.download_gtdb_metadata.output.bac120,
        ar53=rules.download_gtdb_metadata.output.ar53,
        mmseqs2=rules.download_mmseqs2.output[0],
        anvio=rules.download_anvio_scg.output[0],
    output:
        f"{DB_DIR}/database_paths.yaml",
    log:
        f"{LOGS}/databases/paths.log",
    conda:
        env("helpers")
    container:
        image("helpers")
    params:
        db_dir=lambda wc, output: os.path.dirname(output[0]),
        mmseqs=f"{DB_DIR}/mmseqs2/{MMSEQS_DIR}/db",
        padded=(
            "_padded" if (USE_GPU and config["binning"]["vamb"]["mmseqs_gpu"]) else ""
        ),
    shell:
        r"""
        gtdb=$(ls -d {params.db_dir}/gtdbtk/*/ | head -n 1)
        cat >{output} <<EOF
# Copy into config/config.yaml
databases:
  checkm2: $(readlink -f {input.checkm2})
  gtdbtk: $(readlink -f "$gtdb")
  mmseqs2_taxonomy: $(readlink -f {params.mmseqs}){params.padded}
  anvio_scg: $(readlink -f {params.db_dir}/anvio_scg/data)
ena:
  gtdb_metadata:
    bac120: $(readlink -f {input.bac120})
    ar53: $(readlink -f {input.ar53})
EOF
        cat {output} >{log}
        """
