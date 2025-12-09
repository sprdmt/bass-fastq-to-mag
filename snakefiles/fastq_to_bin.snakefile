import os
import glob

configfile: "snakemake-config/project_config.yaml"
PROJECT = config["project"]
FASTQ_DIR = config.get("fastq_dir", "fastq")

SAMPLES = sorted([
    os.path.basename(f).replace("_1.fastq.gz", "")
    for f in glob.glob(os.path.join(FASTQ_DIR, "*_1.fastq.gz"))
])

print("Detected samples:", SAMPLES)

#######################################################################
# Final targets
#######################################################################
rule all:
    input:
        # Assembly and graph
        coassembly = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        gfa        = f"assembly/coassembly/{PROJECT}_coassembly_29_127.gfa",

        # Gene prediction
        prodigal_faa = f"genes/prodigal/{PROJECT}_genes.faa",
        prodigal_fna = f"genes/prodigal/{PROJECT}_genes.fna",
        prodigal_gff = f"genes/prodigal/{PROJECT}_genes.gff",

        # Bowtie2 mappings on contigs
        coassembly_idx = [
            f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.{x}.bt2l"
            for x in [1, 2, 3, 4, 'rev.1', 'rev.2']
        ],
        coassembly_bam = expand("mapping/bowtie2/coassembly/bam/{sample}.bam", sample=SAMPLES),
        coassembly_abundance = expand("mapping/bowtie2/coassembly/bam/{sample}_abundance.txt", sample=SAMPLES),

        # GraphBin2 prerequisites and results (markers only)
        graphbin2_abundance = "binning/graphbin2/abundance.tsv",
        graphbin2_metabat2_done   = "binning/graphbin2/metabat2/.done",
        graphbin2_genomeface_done = "binning/graphbin2/genomeface/.done",
        graphbin2_semibin2_done   = "binning/graphbin2/semibin2/.done",
        graphbin2_vamb_done       = "binning/graphbin2/vamb/.done",

        # Binette → dRep → anvio
        binette_final_done = "binning/binette/final_bins/.done",
        drep_done          = "binning/dRep/binette/.done",
        anvio_contigs_db      = "anvio/contigs/contigs.db",
        anvio_profiles        = expand("anvio/profile/single/{sample}/PROFILE.db", sample=SAMPLES),
        anvio_profile_merged  = "anvio/profile/merged/PROFILE.db",
        anvio_bins_original   = "anvio/bins/non_refined/.done",
        anvio_bins_imported   = "anvio/profile/merged/collection_OriginalBins.done",
        anvio_hmms_done       = "anvio/contigs/hmms.done",
        anvio_taxonomy_done   = "anvio/contigs/scg_taxonomy.done"


#######################################################################
# Assembly and graph
#######################################################################
rule coassemble_reads:
    input:
        r1 = expand(os.path.join(FASTQ_DIR, "{sample}_1.fastq.gz"), sample=SAMPLES),
        r2 = expand(os.path.join(FASTQ_DIR, "{sample}_2.fastq.gz"), sample=SAMPLES)
    output:
        assembly = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    log:
        f"logs/coassemble_reads/{PROJECT}.log"
    params:
        tmpdir = lambda wc: os.environ.get("TMPDIR", "/tmp"),
        r1_csv = lambda wildcards, input: ",".join(input.r1),
        r2_csv = lambda wildcards, input: ",".join(input.r2)
    shell:
        r"""
        mkdir -p logs/coassemble_reads assembly/coassembly
        ml MEGAHIT
        megahit \
            -1 {params.r1_csv} \
            -2 {params.r2_csv} \
            -o {params.tmpdir}/megahit_out \
            --k-min 25 --k-max 125 --k-step 10 \
            --min-contig-len 1000 \
            --num-cpu-threads {resources.cpus_per_task} \
            --memory 0.95 \
            --min-count 2 \
            --force \
            --tmp-dir $TMPDIR 
            
        cp {params.tmpdir}/megahit_out/final.contigs.fa {output.assembly}
        
        """

rule megahit_contig2fastg:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        fastg = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fastg"
    log:
        f"logs/graphbuild/{PROJECT}.contig2fastg.log"
    shell:
        r"""
        mkdir -p logs/graphbuild assembly/coassembly
        ml MEGAHIT
        megahit_toolkit contig2fastg 127 {input.contigs} > {output.fastg} 2> {log}
        """

rule agtools_fastg2gfa:
    input:
        fastg = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fastg"
    output:
        gfa = f"assembly/coassembly/{PROJECT}_coassembly_29_127.gfa"
    log:
        f"logs/graphbuild/{PROJECT}.fastg2gfa.log"
    shell:
        r"""
        ml Conda && conda activate agtools
        agtools fastg2gfa -g {input.fastg} -k 127 -o assembly/coassembly/ &> {log}
        if [ -s assembly/coassembly/converted_graph.gfa ]; then
            mv -f assembly/coassembly/converted_graph.gfa {output.gfa}
        else
            echo "agtools did not produce converted_graph.gfa" >> {log}
            exit 2
        fi
        """

#######################################################################
# Gene prediction
#######################################################################
rule predict_genes_pprodigal:
    input:
        assembly = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        faa = f"genes/prodigal/{PROJECT}_genes.faa",
        fna = f"genes/prodigal/{PROJECT}_genes.fna",
        gff = f"genes/prodigal/{PROJECT}_genes.gff"
    log:
        f"logs/predict_genes_pprodigal/{PROJECT}.log"
    params:
        tasks = 40,
        chunksize = 20000
    shell:
        """
        mkdir -p genes/prodigal logs/predict_genes_pprodigal
        ml Conda
        conda activate pprodigal
        pprodigal \
            --tasks {params.tasks} \
            --chunksize {params.chunksize} \
            -i {input.assembly} \
            -a {output.faa} \
            -d {output.fna} \
            -f gff \
            -o {output.gff} \
            > {log} 2>&1
        echo "pprodigal finished on {input.assembly}" >> {log}
        """

#######################################################################
# Mapping and coverage
#######################################################################
rule build_bowtie2_coassembly_index:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        idx1 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.1.bt2l",
        idx2 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.2.bt2l",
        idx3 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.3.bt2l",
        idx4 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.4.bt2l",
        idxr1 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.rev.1.bt2l",
        idxr2 = f"mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127.rev.2.bt2l"
    log:
        f"logs/build_index/{PROJECT}_coassembly.log"
    shell:
        r"""
        mkdir -p mapping/bowtie2/coassembly/index logs/build_index
        ml Bowtie2
        bowtie2-build {input.contigs} mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127 --threads {resources.cpus_per_task} --large-index &> {log}
        """

rule map_reads_coassembly:
    input:
        idx1 = rules.build_bowtie2_coassembly_index.output.idx1,
        r1 = lambda wc: os.path.join(FASTQ_DIR, f"{wc.sample}_1.fastq.gz"),
        r2 = lambda wc: os.path.join(FASTQ_DIR, f"{wc.sample}_2.fastq.gz")
    output:
        bam = "mapping/bowtie2/coassembly/bam/{sample}.bam",
        bai = "mapping/bowtie2/coassembly/bam/{sample}.bam.bai",
        abundance = "mapping/bowtie2/coassembly/bam/{sample}_abundance.txt"
    log:
        "logs/bowtie2_coassembly/{sample}.log"
    shell:
        r"""
        mkdir -p mapping/bowtie2/coassembly/bam logs/bowtie2_coassembly
        ml Bowtie2
        ml SAMtools

        bowtie2 -x mapping/bowtie2/coassembly/index/{PROJECT}_coassembly_29_127 \
          -1 {input.r1} -2 {input.r2} \
          -S ${{TMPDIR}}/{wildcards.sample}.sam \
          --no-unal --threads {resources.cpus_per_task} &> {log}

        samtools view -@ {resources.cpus_per_task} -b ${{TMPDIR}}/{wildcards.sample}.sam \
          > ${{TMPDIR}}/{wildcards.sample}-raw.bam

        samtools sort -@ {resources.cpus_per_task} \
          -o {output.bam} ${{TMPDIR}}/{wildcards.sample}-raw.bam

        samtools index {output.bam} {output.bai}
        samtools idxstats {output.bam} > {output.abundance}

        rm -f ${{TMPDIR}}/{wildcards.sample}.sam ${{TMPDIR}}/{wildcards.sample}-raw.bam
        """


rule coverm_contig_coverage:
    input:
        bams = expand("mapping/bowtie2/coassembly/bam/{sample}.bam", sample=SAMPLES)
    output:
        tsv = "binning/graphbin2/abundance.tsv"
    log:
        f"logs/graphbin2/{PROJECT}.coverm.log"
    shell:
        r"""
        mkdir -p binning/graphbin2 logs/graphbin2
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
        
        
        ml Conda
        conda activate CoverM-0.7.0
        coverm contig \
          -b mapping/bowtie2/coassembly/bam/*.bam \
          -o {output.tsv} \
          -t {resources.cpus_per_task} \
          &> {log}
          
          
        sed -i '1d' {output.tsv}
        """

#######################################################################
# Binning backends
#######################################################################
rule jgi_depth:
    input:
        assembly = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        bams = expand("mapping/bowtie2/coassembly/bam/{sample}.bam", sample=SAMPLES)
    output:
        depth = f"binning/metabat2/{PROJECT}_depth.txt"
    log:
        f"logs/binning/rgi_depth/{PROJECT}.log"
    shell:
        """
        ml MetaBAT
        mkdir -p binning/metabat2
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            --referenceFasta {input.assembly} \
            {input.bams}
        """

rule metabat2_binning:
    input:
        assembly = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        depth = f"binning/metabat2/{PROJECT}_depth.txt"
    output:
        done = f"binning/metabat2/{PROJECT}/bins/.done"
    log:
        f"logs/binning/metabat2_binning/{PROJECT}.log"
    threads: 24
    shell:
        """
        set -eo pipefail
        ml MetaBAT
        
        NCORES={resources.cpus_per_task}
        
        mkdir -p binning/metabat2/{PROJECT}/bins
        metabat2 \
            -i {input.assembly} \
            -a {input.depth} \
            -o binning/metabat2/{PROJECT}/bins/bin \
            -m 1500 -s 100000 -t {threads}
        for f in binning/metabat2/{PROJECT}/bins/*.fa; do
          awk '/^>/ {{split($0, a, /[\\t ]/); print a[1]; next}} {{print}}' "$f" > "${{f}}.fixed"
          mv "${{f}}.fixed" "${{f}}"
        done
        rm binning/metabat2/{PROJECT}/bins/*BinInfo.txt
        touch {output.done}
        """

rule rename_metabat2_bins:
    input:
        marker = f"binning/metabat2/{PROJECT}/bins/.done"
    output:
        renamed = f"binning/metabat2/{PROJECT}/bins/.renamed"
    run:
        import glob, shutil, os
        bins_dir = f"binning/metabat2/{PROJECT}/bins"
        files = sorted(glob.glob(f"{bins_dir}/*.fa"))
        for i, f in enumerate(files, start=1):
            target = f"{bins_dir}/bin.{i}.fa"
            if os.path.abspath(f) != os.path.abspath(target):
                shutil.move(f, target)
        open(output.renamed, "w").close()

# Configuration for parallel marker search
MARKER_CHUNKS = 24
THREADS_PER_CHUNK = 4

# Split contigs into chunks
rule split_contigs_for_markers:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        chunks = expand(f"binning/genomeface/{PROJECT}/chunks/chunk_{{chunk}}.fasta", 
                       chunk=range(1, MARKER_CHUNKS + 1)),
        done = f"binning/genomeface/{PROJECT}/chunks/.split_done"
    log:
        f"logs/binning/gf_split/{PROJECT}.log"
    params:
        n_chunks = MARKER_CHUNKS,
        chunk_dir = f"binning/genomeface/{PROJECT}/chunks"
    shell:
        r"""
        set -euo pipefail
        
        # Load seqkit using conda stack
        ml SeqKit
        
        mkdir -p {params.chunk_dir}
        
        # Split into chunks using seqkit in TMPDIR
        seqkit split -p {params.n_chunks} -O $TMPDIR/seqkit_split {input.contigs} 2>&1 | tee {log}
        
        # Rename and copy chunks to final location
        i=1
        for f in $TMPDIR/seqkit_split/*.fasta; do
            if [ -f "$f" ]; then
                cp "$f" {params.chunk_dir}/chunk_$i.fasta
                i=$((i+1))
            fi
        done
        
        # Clean up temp files
        rm -rf $TMPDIR/seqkit_split
        
        touch {output.done}
        echo "Split into {params.n_chunks} chunks" >> {log}
        """

# Run markersgf on each chunk in parallel
rule genomeface_markers_chunk:
    input:
        chunk = f"binning/genomeface/{PROJECT}/chunks/chunk_{{chunk}}.fasta"
    output:
        markers = f"binning/genomeface/{PROJECT}/markers_chunks/chunk_{{chunk}}.tsv"
    threads: THREADS_PER_CHUNK
    log:
        f"logs/binning/gf_markers/{PROJECT}_chunk_{{chunk}}.log"
    shell:
        r"""
        set -euo pipefail

        # 1) Remove path pollution from loaded modules
        module purge || true

        # 2) Sanitize Python import behavior
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1

        # 3) Activate the intended env
        ml Conda
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate genomeface

        # 4) Preflight check
        python - <<'PY'
import sys, numpy, Bio
print("PY:", sys.executable)
print("NP FILE:", getattr(numpy, "__file__", None))
print("NP VER:", numpy.__version__)
print("BIO:", Bio.__version__)
assert "/app/snakemake/" not in (numpy.__file__ or ""), "Shadowed by Snakemake module"
PY

        # 5) Create output directories
        mkdir -p binning/genomeface/{PROJECT}/markers_chunks
        mkdir -p logs/binning/gf_markers
        
        # Store absolute paths before changing directory
        WORKDIR=$(pwd)
        INPUT_ABS="$WORKDIR/{input.chunk}"
        OUTPUT_ABS="$WORKDIR/{output.markers}"
        LOG_ABS="$WORKDIR/{log}"
        
        # Copy input to TMPDIR and run there
        cp "$INPUT_ABS" $TMPDIR/chunk_{wildcards.chunk}.fasta
        
        # Run markersgf in TMPDIR
        cd $TMPDIR
        env -u PYTHONPATH PYTHONNOUSERSITE=1 \
          markersgf -t {threads} chunk_{wildcards.chunk}.fasta chunk_{wildcards.chunk}.tsv 2>&1 | tee "$LOG_ABS"
        
        # Copy result back
        cp $TMPDIR/chunk_{wildcards.chunk}.tsv "$OUTPUT_ABS"
        
        # Clean up
        rm -f $TMPDIR/chunk_{wildcards.chunk}.fasta $TMPDIR/chunk_{wildcards.chunk}.tsv
        """

# Merge all marker chunks into single TSV
rule merge_markers:
    input:
        chunks = expand(f"binning/genomeface/{PROJECT}/markers_chunks/chunk_{{chunk}}.tsv",
                       chunk=range(1, MARKER_CHUNKS + 1))
    output:
        markers = f"binning/genomeface/{PROJECT}/markers.tsv"
    log:
        f"logs/binning/gf_merge/{PROJECT}.log"
    shell:
        """
        set -euo pipefail
        
        # Store absolute paths before changing directory
        WORKDIR=$(pwd)
        OUTPUT_ABS="$WORKDIR/{output.markers}"
        LOG_ABS="$WORKDIR/{log}"
        
        # Convert all input paths to absolute
        INPUT_ABS=""
        for chunk in {input.chunks}; do
            INPUT_ABS="$INPUT_ABS $WORKDIR/$chunk"
        done
        
        # Create output directory
        mkdir -p binning/genomeface/{PROJECT}
        mkdir -p logs/binning/gf_merge
        
        # Work in TMPDIR for merging
        cd $TMPDIR
        
        # Simple concatenation using absolute paths
        cat $INPUT_ABS > merged_markers.tsv
        
        # Copy final result back
        cp merged_markers.tsv "$OUTPUT_ABS"
        
        # Log stats
        echo "Total chunks: {MARKER_CHUNKS}" > "$LOG_ABS"
        echo "Merged lines: $(wc -l < merged_markers.tsv)" >> "$LOG_ABS"
        
        # Clean up
        rm -f merged_markers.tsv
        """



rule run_genomeface:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        depth   = f"binning/metabat2/{PROJECT}_depth.txt",
        markers = f"binning/genomeface/{PROJECT}/markers.tsv"
    output:
        done = f"binning/genomeface/{PROJECT}/.done"
    log:
        f"logs/binning/genomeface/{PROJECT}.log"
    shell:
        """
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
        
        ml Conda
        conda activate genomeface
        genomeface -i {input.contigs} -a {input.depth} -g {input.markers} -o binning/genomeface/{PROJECT}/bins
        mv binning/genomeface/{PROJECT}/bins/bins.tsv binning/genomeface/{PROJECT}/
        touch {output.done}
        """

rule rename_genomeface_bins:
    input:
        marker = f"binning/genomeface/{PROJECT}/.done"
    output:
        renamed = f"binning/genomeface/{PROJECT}/.renamed"
    run:
        import glob, shutil, os
        bins_dir = f"binning/genomeface/{PROJECT}/bins"
        files = sorted(
            f for f in glob.glob(f"{bins_dir}/*") if os.path.isfile(f) and '.' not in os.path.basename(f)
        )
        for i, f in enumerate(files, start=1):
            target = f"{bins_dir}/bin.{i}.fa"
            if os.path.abspath(f) != os.path.abspath(target):
                shutil.move(f, target)
        open(output.renamed, "w").close()

rule semibin2_single_easy_bin:
    input:
        contigs   = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        bam_files = expand("mapping/bowtie2/coassembly/bam/{sample}.bam", sample=SAMPLES)
    output:
        done = f"binning/semibin2/{PROJECT}/single_easy/.done"
    log:
        f"logs/binning/semibin2_single_easy_bin/{PROJECT}.log"
    shell:
        """
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
        
        
        ml Conda
        conda activate SemiBin-2.2.0
        mkdir -p binning/semibin2/{PROJECT}/single_easy
        SemiBin2 single_easy_bin \
          -i {input.contigs} \
          --engine gpu \
          -o binning/semibin2/{PROJECT}/single_easy \
          -p {resources.cpus_per_task} \
          --tmpdir $TMPDIR/semibin2 \
          -m 2500 \
          -b {input.bam_files}
        touch {output.done}
        """

rule rename_semibin2_bins:
    input:
        marker = f"binning/semibin2/{PROJECT}/single_easy/.done"
    output:
        renamed = f"binning/semibin2/{PROJECT}/single_easy/.renamed"
    run:
        import glob, shutil, os, subprocess
        bins_dir = f"binning/semibin2/{PROJECT}/single_easy/output_bins"
        gz_files = sorted(glob.glob(os.path.join(bins_dir, "*.fa.gz")))
        for gz_file in gz_files:
            fa_file = gz_file[:-3]
            subprocess.run(["gunzip", "-f", gz_file], check=False)
            assert os.path.exists(fa_file), f"Decompression failed for {gz_file}"
        files = sorted(glob.glob(os.path.join(bins_dir, "*.fa")))
        for i, f in enumerate(files, start=1):
            target = os.path.join(bins_dir, f"bin.{i}.fa")
            if os.path.abspath(f) != os.path.abspath(target):
                shutil.move(f, target)
        open(output.renamed, "w").close()

#######################################################################
# TaxVAMB
#######################################################################
rule mmseqs2_taxonomy:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        tsv         = f"binning/vamb/{PROJECT}/taxonomy.tsv",
        tsv_taxconv = f"binning/vamb/{PROJECT}/taxonomy_taxconv.tsv"
    params:
        mmseqs_uniref_db = "/lisc/data/scratch/oceanography/spriahailo/db/mmseqs2/uniref90/db"
    log:
        f"logs/binning/mmseqs2_taxonomy/{PROJECT}.log"
    shell:
        """
        
        module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
        
        ml MMseqs2
        #source /lisc/app/mmseqs2/17-b804f/util/bash-completion.sh
        mkdir -p binning/vamb/{PROJECT}/mmseqs
        mmseqs createdb {input.contigs} binning/vamb/{PROJECT}/mmseqs/contigDB
        mmseqs taxonomy binning/vamb/{PROJECT}/mmseqs/contigDB {params.mmseqs_uniref_db} \
            binning/vamb/{PROJECT}/mmseqs/result $TMPDIR \
            --threads {resources.cpus_per_task} \
            --gpu $CUDA_VISIBLE_DEVICES \
            --tax-lineage 1 \
            --lca-ranks superkingdom,phylum,class,order,family,genus,species \
            --max-seq-len 1500000
        mmseqs createtsv binning/vamb/{PROJECT}/mmseqs/contigDB binning/vamb/{PROJECT}/mmseqs/result {output.tsv}
        
        module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
        ml Conda
        conda activate taxconv
        taxconverter mmseqs2 -i {output.tsv} -o {output.tsv_taxconv}
        """

rule fix_taxonomy_for_vamb:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        taxconv = f"binning/vamb/{PROJECT}/taxonomy_taxconv.tsv",
    output:
        fixed = f"binning/vamb/{PROJECT}/taxonomy_fixed.tsv",
    log:
        f"logs/binning/fix_taxonomy_for_vamb/{PROJECT}.log",
    shell:
        r"""
        ml Python
        python - <<'PY' "{input.contigs}" "{input.taxconv}" "{output.fixed}" "{log}"
import sys

fasta, taxconv, outp, logp = sys.argv[1:5]

# 1) Read contig IDs from FASTA (first token of header)
contigs = []
with open(fasta, 'r') as f:
    for ln in f:
        if ln.startswith('>'):
            contigs.append(ln[1:].strip().split()[0])
n_fa = len(contigs)

# 2) Load taxonomy and minimally disambiguate tokens that occur under multiple parents
first_parent = dict()   # token -> parent_token of its FIRST occurrence
taxmap = dict()
n_rows = 0
n_ambig = 0

def clean_lineage(raw):
    global n_ambig
    if not raw:
        return ""
    toks = [t.strip() for t in raw.split(';') if t.strip()]
    if not toks:
        return ""
    # stop at unknown / uc_*
    kept = []
    for t in toks:
        if t == "unknown" or t.startswith("uc_"):
            break
        kept.append(t)
    if not kept:
        return ""
    out = []
    for i in range(len(kept)):
        tok = kept[i]
        base_parent = kept[i-1] if i > 0 else "root"
        prev = first_parent.get(tok)
        if prev is None:
            first_parent[tok] = base_parent
            out.append(tok)
        elif prev == base_parent:
            out.append(tok)
        else:
            out.append(base_parent + "_" + tok)
            n_ambig += 1
    return ";".join(out)

# taxconv header usually: "contigs\tpredictions"
with open(taxconv, 'r') as f:
    _header = f.readline()
    for ln in f:
        if not ln.strip():
            continue
        n_rows += 1
        parts = ln.rstrip('\n').split('\t')
        name = parts[0].strip()
        lineage = parts[1].strip() if len(parts) > 1 else ""
        taxmap[name] = clean_lineage(lineage)

# 3) Write fixed file for ALL FASTA contigs (preserve order)
missing = 0
with open(outp, 'w') as out:
    out.write("contigs\tpredictions\n")
    for name in contigs:
        if name in taxmap:
            out.write(name + "\t" + taxmap[name] + "\n")
        else:
            missing += 1
            out.write(name + "\t\n")

with open(logp, 'w') as lg:
    lg.write("FASTA contigs: " + str(n_fa) + "\n")
    lg.write("Taxonomy rows read: " + str(n_rows) + "\n")
    lg.write("Missing contigs filled: " + str(missing) + "\n")
    lg.write("Ambiguous tokens disambiguated: " + str(n_ambig) + "\n")
PY
        """


rule run_taxvamb:
    input:
        contigs  = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        taxonomy = f"binning/vamb/{PROJECT}/taxonomy_fixed.tsv",
        bams     = expand("mapping/bowtie2/coassembly/bam/{sample}.bam", sample=SAMPLES)
    output:
        done = f"binning/vamb/{PROJECT}/.done"
    log:
        f"logs/binning/taxvamb/{PROJECT}.log"
    shell:
        """
        set -euo pipefail
        
        # Store absolute paths before changing directory
        WORKDIR=$(pwd)
        CONTIGS_ABS="$WORKDIR/{input.contigs}"
        TAXONOMY_ABS="$WORKDIR/{input.taxonomy}"
        OUTPUT_DIR_ABS="$WORKDIR/binning/vamb/{PROJECT}"
        LOG_ABS="$WORKDIR/{log}"
        
        # Create output directories in working directory
        mkdir -p "$OUTPUT_DIR_ABS" "$(dirname $LOG_ABS)"
        
        # Remove path pollution from loaded modules
        module purge || true
        
        # Sanitize Python import behavior
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        
        ml vamb
        
        # Create unique temporary directory in $TMPDIR
        TMPWORK=$(mktemp -d "${{TMPDIR}}/taxvamb.XXXXXXXXXX")
        trap "rm -rf '$TMPWORK'" EXIT
        
        echo "Working in temporary directory: $TMPWORK" > "$LOG_ABS"
        echo "Start time: $(date)" >> "$LOG_ABS"
        
        # Copy input files to temporary directory
        echo "Copying input files to $TMPDIR..." >> "$LOG_ABS"
        cp "$CONTIGS_ABS" "$TMPWORK/contigs.fasta"
        cp "$TAXONOMY_ABS" "$TMPWORK/taxonomy.tsv"
        
        # Copy BAM files
        mkdir -p "$TMPWORK/bams"
        for bam in {input.bams}; do
            BAM_ABS="$WORKDIR/$bam"
            cp "$BAM_ABS" "$TMPWORK/bams/"
            echo "  Copied $(basename $bam)" >> "$LOG_ABS"
        done
        
        echo "All input files copied successfully" >> "$LOG_ABS"
        
        # Change to temporary directory
        cd "$TMPWORK"
        
        # Run TaxVamb with local paths
        echo "Starting TaxVamb training..." >> "$LOG_ABS"
        vamb bin taxvamb \
            --fasta contigs.fasta \
            --bamdir bams \
            --taxonomy taxonomy.tsv \
            --outdir taxvamb_bins \
            --cuda \
            -p {resources.cpus_per_task} \
            -m 1000 \
            -o "" \
            --minfasta 200000 \
            2>&1 | tee -a "$LOG_ABS"
        
        echo "TaxVamb completed successfully" >> "$LOG_ABS"
        
        # Copy results back to working directory
        echo "Copying results back to $OUTPUT_DIR_ABS..." >> "$LOG_ABS"
        
        # Copy the entire output directory
        if [[ -d taxvamb_bins ]]; then
            cp -r taxvamb_bins "$OUTPUT_DIR_ABS/"
            
            # Count bins
            BIN_COUNT=$(find taxvamb_bins -name "*.fna" -o -name "*.fa" -o -name "*.fasta" | wc -l)
            echo "Copied $BIN_COUNT bin files" >> "$LOG_ABS"
        else
            echo "Warning: taxvamb_bins directory not found" >> "$LOG_ABS"
        fi
        
        # Copy any other output files (logs, CSV, TSV, etc.)
        for ext in csv tsv txt log; do
            for f in *.$ext; do
                [[ -f "$f" ]] && cp "$f" "$OUTPUT_DIR_ABS/"
            done
        done
        
        # Create done marker
        touch "$OUTPUT_DIR_ABS/.done"
        
        echo "End time: $(date)" >> "$LOG_ABS"
        echo "TaxVamb workflow completed successfully" >> "$LOG_ABS"
        
        # Clean up happens automatically via trap
        """


rule rename_vamb_bins:
    input:
        marker = f"binning/vamb/{PROJECT}/.done"
    output:
        renamed = f"binning/vamb/{PROJECT}/taxvamb_bins/.renamed"
    run:
        import glob, shutil, os
        bins_dir = f"binning/vamb/{PROJECT}/taxvamb_bins/bins"
        files = sorted(glob.glob(f"{bins_dir}/*.f*a"))
        for i, f in enumerate(files, start=1):
            target = f"{bins_dir}/bin.{i}.fa"
            if os.path.abspath(f) != os.path.abspath(target):
                shutil.move(f, target)
        open(output.renamed, "w").close()

#######################################################################
# GraphBin2 (generic rules)
#######################################################################
def _gb_bins_dir(binner):
    if binner == "metabat2":
        return f"binning/metabat2/{PROJECT}/bins"
    if binner == "genomeface":
        return f"binning/genomeface/{PROJECT}/bins"
    if binner == "semibin2":
        return f"binning/semibin2/{PROJECT}/single_easy/output_bins"
    if binner == "vamb":
        return f"binning/vamb/{PROJECT}/taxvamb_bins/bins"
    raise ValueError(binner)

def _gb_marker(binner):
    if binner == "metabat2":
        return f"binning/metabat2/{PROJECT}/bins/.renamed"
    if binner == "genomeface":
        return f"binning/genomeface/{PROJECT}/.renamed"
    if binner == "semibin2":
        return f"binning/semibin2/{PROJECT}/single_easy/.renamed"
    if binner == "vamb":
        return f"binning/vamb/{PROJECT}/taxvamb_bins/.renamed"
    raise ValueError(binner)

#######################################################################
# GraphBin2: prepResult -> {binner}_initial_contig_bins.csv
#######################################################################
rule graphbin2_prep:
    input:
        marker = lambda wc: _gb_marker(wc.binner)
    output:
        csv = "binning/graphbin2/{binner}_initial_contig_bins.csv"
    log:
        "logs/graphbin2/" + PROJECT + ".{binner}.prep.log"
    params:
        outdir = "binning/graphbin2",
        prefix = "{binner}",
        bins_dir = lambda wc: _gb_bins_dir(wc.binner)
    shell:
        r"""
        mkdir -p binning/graphbin2 logs/graphbin2
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda && conda activate graphbin2-patched
        prepResult \
          --binned {params.bins_dir} \
          --output {params.outdir} \
          --prefix {params.prefix} \
          &> {log}
        test -s {output.csv}
        """

#######################################################################
# GraphBin2: run refinement per binner
#######################################################################
rule graphbin2_run:
    input:
        contigs = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta",
        graph   = f"assembly/coassembly/{PROJECT}_coassembly_29_127.gfa",
        binned  = "binning/graphbin2/{binner}_initial_contig_bins.csv",
        abund   = "binning/graphbin2/abundance.tsv"
    output:
        csv  = "binning/graphbin2/{binner}/graphbin2_output.csv",
        done = "binning/graphbin2/{binner}/.done"
    log:
        "logs/graphbin2/" + PROJECT + ".{binner}.run.log"
    params:
        outdir = "binning/graphbin2/{binner}"
    threads: 12
    shell:
        r"""
        set -euo pipefail
        
        # Store absolute paths before changing directory
        WORKDIR=$(pwd)
        CONTIGS_ABS="$WORKDIR/{input.contigs}"
        GRAPH_ABS="$WORKDIR/{input.graph}"
        BINNED_ABS="$WORKDIR/{input.binned}"
        ABUND_ABS="$WORKDIR/{input.abund}"
        OUTPUT_DIR_ABS="$WORKDIR/{params.outdir}"
        LOG_ABS="$WORKDIR/{log}"
        
        # Create output directories in working directory
        mkdir -p "$OUTPUT_DIR_ABS" "$(dirname $LOG_ABS)"
        
        # Remove path pollution from loaded modules
        module purge || true
        
        # Sanitize Python import behavior
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        
        ml Conda
        conda activate graphbin2-patched
        
        # Create unique temporary directory in $TMPDIR
        TMPWORK=$(mktemp -d "${{TMPDIR}}/graphbin2.XXXXXXXXXX")
        trap "rm -rf '$TMPWORK'" EXIT
        
        echo "Working in temporary directory: $TMPWORK" > "$LOG_ABS"
        
        # Copy input files to temporary directory
        echo "Copying input files to $TMPDIR..." >> "$LOG_ABS"
        cp "$CONTIGS_ABS" "$TMPWORK/contigs.fasta"
        cp "$GRAPH_ABS" "$TMPWORK/graph.gfa"
        cp "$BINNED_ABS" "$TMPWORK/binned.csv"
        cp "$ABUND_ABS" "$TMPWORK/abundance.tsv"
        
        # Change to temporary directory
        cd "$TMPWORK"
        
        # 1) Trim + dedupe (keep bin column intact; detect delimiter != space)
        BINNED_TRIMMED="binned.trimmed.csv"
        CHANGED_COUNT="$(
        python - <<'PY' "binned.csv" "$BINNED_TRIMMED"
import sys,csv
if len(sys.argv)!=3:
    sys.stderr.write('usage: IN OUT\n'); sys.exit(2)
inp=sys.argv[1]; outp=sys.argv[2]
def first_nonempty_line(p):
    f=open(p,'r',newline='')
    try:
        for ln in f:
            if ln.strip():
                return ln.rstrip('\n')
    finally:
        f.close()
    return ''
probe=first_nonempty_line(inp)
delim=','
for d in (',','\t',';','|'):
    if d in probe:
        delim=d; break
seen=set(); total=0; changed=0; dup=0
fin=open(inp,'r',newline=''); fout=open(outp,'w',newline='')
try:
    r=csv.reader(fin,delimiter=delim)
    w=csv.writer(fout,delimiter=',',quotechar='"',lineterminator='\n')
    for row in r:
        if not row: continue
        total+=1
        if len(row)<2: row=(row+[''])[:2]
        contig_raw=row[0]; bin_raw=row[1]
        contig_trim=contig_raw.split()[0] if contig_raw else contig_raw
        if contig_trim!=contig_raw: changed+=1
        key=(contig_trim,bin_raw)
        if key in seen: dup+=1; continue
        seen.add(key)
        w.writerow([contig_trim,bin_raw])
finally:
    fin.close(); fout.close()
sys.stderr.write('[trim_bins] Input: '+inp+'\n')
sys.stderr.write('[trim_bins] Output: '+outp+'\n')
sys.stderr.write('[trim_bins] Rows read: '+str(total)+', trimmed: '+str(changed)+', deduplicated: '+str(dup)+'\n')
print(changed)
PY
        )"
        
        # 2) Safe switch: compute contig_count from FASTA and drop problematic rows
        BINNED_SAFE="binned.safe.csv"
        python - <<'PY' "contigs.fasta" "$BINNED_TRIMMED" "$BINNED_SAFE"
import sys,re,csv
fa=sys.argv[1]; inp=sys.argv[2]; outp=sys.argv[3]
contig_count=0
f=open(fa,'r')
try:
    for ln in f:
        if ln.startswith('>'):
            contig_count+=1
finally:
    f.close()
dropped=0; kept=0
fin=open(inp,'r',newline=''); fout=open(outp,'w',newline='')
try:
    r=csv.reader(fin,delimiter=',')
    w=csv.writer(fout,delimiter=',',quotechar='"',lineterminator='\n')
    for row in r:
        if not row: continue
        name=row[0]; binname=row[1] if len(row)>1 else ''
        # accept "210535" or "k127_210535"
        m=re.match(r'(?:k\d+_)?(\d+)$', name)
        num=int(m.group(1)) if m else None
        if num is not None and num==contig_count:
            dropped+=1
            continue
        w.writerow([name,binname]); kept+=1
finally:
    fin.close(); fout.close()
sys.stderr.write('[safe_drop] contig_count='+str(contig_count)+' dropped='+str(dropped)+' kept='+str(kept)+'\n')
PY
        
        BINNED_RUNIN="$BINNED_SAFE"
        
        # 3) Import-time patch to avoid NameError in GB2 logging
        PYFIX_DIR="_pyfix"; mkdir -p "$PYFIX_DIR"
        cat > "$PYFIX_DIR/sitecustomize.py" <<'PY'
from collections import defaultdict
try:
    import graphbin2.graphbin2_MEGAHIT as gb
    if not hasattr(gb, 'contig_names'):
        gb.contig_names = defaultdict(lambda: '<contig>')
except Exception:
    pass
PY
        export PYTHONPATH="$PYFIX_DIR:${{PYTHONPATH:-}}"
        export OMP_NUM_THREADS="{threads}"
        
        # 4) Run GraphBin2 in local temporary directory
        echo "Running GraphBin2..." >> "$LOG_ABS"
        graphbin2 \
          --assembler megahit \
          --contigs   "contigs.fasta" \
          --graph     "graph.gfa" \
          --binned    "$BINNED_RUNIN" \
          --output    "./" \
          --abundance "abundance.tsv" \
          --nthreads  {resources.cpus_per_task} \
          &> graphbin2.log
        
        # 5) If we trimmed IDs, normalize headers in per-bin FASTAs
        if [[ "$CHANGED_COUNT" =~ ^[0-9]+$ ]] && [[ "$CHANGED_COUNT" -gt 0 ]]; then
          BINS_DIR="bins"
          if [[ -d "$BINS_DIR" ]]; then
            python - <<'PY' "$BINS_DIR"
import sys,os
from pathlib import Path
bins_dir=Path(sys.argv[1])
def clean(path):
    tmp=Path(str(path)+'.tmp'); total=0; changed=0
    fin=open(path,'r'); fout=open(tmp,'w')
    try:
        for ln in fin:
            total+=1
            if ln.startswith('>'):
                hdr=ln[1:].rstrip('\n')
                trimmed=hdr.split()[0] if hdr else hdr
                if trimmed!=hdr: changed+=1
                fout.write('>'+trimmed+'\n')
            else:
                fout.write(ln)
    finally:
        fin.close(); fout.close()
    os.replace(str(tmp),str(path))
    return total,changed
files=list(bins_dir.glob('*.fasta'))+list(bins_dir.glob('*.fa'))
for f in files:
    t,c=clean(f)
    if c>0:
        sys.stderr.write('[clean_bins] '+f.name+' headers_changed='+str(c)+'\n')
sys.stderr.write('[clean_bins] files='+str(len(files))+'\n')
PY
          fi
        fi
        
        # 6) Prune empty FASTAs
        BINS_DIR="bins"
        if [[ -d "$BINS_DIR" ]]; then
          python - <<'PY' "$BINS_DIR"
import sys
from pathlib import Path
bins_dir=Path(sys.argv[1]); checked=0; deleted=0
files=list(bins_dir.glob('*.fasta'))+list(bins_dir.glob('*.fa'))
for f in files:
    checked+=1
    try:
        if f.stat().st_size==0:
            f.unlink(); deleted+=1
            sys.stderr.write('[prune_bins] zero-byte '+f.name+'\n'); continue
        has_header=False
        fin=open(f,'r')
        try:
            for ln in fin:
                if ln.startswith('>'):
                    has_header=True; break
        finally:
            fin.close()
        if not has_header:
            f.unlink(); deleted+=1
            sys.stderr.write('[prune_bins] no headers '+f.name+'\n')
    except Exception as e:
        sys.stderr.write('[prune_bins] warn '+f.name+': '+str(e)+'\n')
sys.stderr.write('[prune_bins] checked='+str(checked)+' deleted='+str(deleted)+'\n')
PY
        fi
        
        # 7) Copy all results back to working directory
        echo "Copying results back to $OUTPUT_DIR_ABS..." >> "$LOG_ABS"
        
        # Copy all output files
        if [[ -d bins ]]; then
            mkdir -p "$OUTPUT_DIR_ABS/bins"
            cp -r bins/* "$OUTPUT_DIR_ABS/bins/"
        fi
        
        # Copy main output files
        for f in *.csv *.tsv *.txt; do
            [[ -f "$f" ]] && cp "$f" "$OUTPUT_DIR_ABS/"
        done
        
        # Append local log to final log
        if [[ -f graphbin2.log ]]; then
            cat graphbin2.log >> "$LOG_ABS"
        fi
        
        # Create done marker
        touch "$OUTPUT_DIR_ABS/.done"
        
        echo "GraphBin2 completed successfully" >> "$LOG_ABS"
        
        # Clean up happens automatically via trap
        """



#######################################################################
# Binette (consume GB2 bins), finalize, dRep, CheckM2
#######################################################################
rule binette_run:
    input:
        contigs  = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fixed.fa",
        proteins = f"genes/prodigal/{PROJECT}_genes.faa",
        gb_metabat2_done   = "binning/graphbin2/metabat2/.done",
        gb_genomeface_done = "binning/graphbin2/genomeface/.done",
        gb_semibin2_done   = "binning/graphbin2/semibin2/.done",
        gb_vamb_done       = "binning/graphbin2/vamb/.done"
    output:
        done = "binning/binette/.done"
    log:
        f"logs/binette/{PROJECT}.log"
    params:
        checkm2_db = "/lisc/scratch/mirror/checkm2/1.1.0/CheckM2_database/uniref100.KO.1.dmnd",
        outdir = "binning/binette",
        gb_metabat2 = "binning/graphbin2/metabat2/bins",
        gb_genomeface = "binning/graphbin2/genomeface/bins",
        gb_semibin2 = "binning/graphbin2/semibin2/bins",
        gb_vamb = "binning/graphbin2/vamb/bins"
    shell:
        r"""
        set -euo pipefail
        
        # Create local log directory
        mkdir -p logs/binette
        
        # Set up temporary working directory
        WORKDIR="${{TMPDIR}}/binette_$$_$RANDOM"
        mkdir -p "$WORKDIR"
        
        # Ensure cleanup on exit
        trap 'rm -rf "$WORKDIR"' EXIT
        
        # Copy input files to TMPDIR
        echo "Copying input files to $WORKDIR..." > {log}
        mkdir -p "$WORKDIR/input"
        cp {input.contigs} "$WORKDIR/input/" 2>> {log}
        cp {input.proteins} "$WORKDIR/input/" 2>> {log}
        
        # Copy bin directories to TMPDIR
        mkdir -p "$WORKDIR/bins"
        cp -r {params.gb_metabat2} "$WORKDIR/bins/metabat2_bins" 2>> {log}
        cp -r {params.gb_genomeface} "$WORKDIR/bins/genomeface_bins" 2>> {log}
        cp -r {params.gb_semibin2} "$WORKDIR/bins/semibin2_bins" 2>> {log}
        cp -r {params.gb_vamb} "$WORKDIR/bins/vamb_bins" 2>> {log}
        
        # Create output directory in TMPDIR
        mkdir -p "$WORKDIR/output"
        
        echo "Starting binette run..." >> {log}
        
        # Remove path pollution from loaded modules
        module purge || true
        
        # Sanitize Python import behavior
        unset PYTHONPATH
        export PYTHONNOUSERSITE=1
        
        # Load environment
        ml Conda
        conda activate binette
        
        # Run binette in TMPDIR
        binette \
          -d "$WORKDIR/bins/metabat2_bins" \
             "$WORKDIR/bins/genomeface_bins" \
             "$WORKDIR/bins/semibin2_bins" \
             "$WORKDIR/bins/vamb_bins" \
          -c "$WORKDIR/input/$(basename {input.contigs})" \
          -p "$WORKDIR/input/$(basename {input.proteins})" \
          -m 70 \
          -t {resources.cpus_per_task} \
          -o "$WORKDIR/output" \
          -e fasta \
          --checkm2_db {params.checkm2_db} \
          -v \
          &>> {log}
        
        echo "Copying results back..." >> {log}
        
        # Copy results back to target directory
        mkdir -p {params.outdir}
        cp -r "$WORKDIR/output"/* {params.outdir}/ 2>> {log}
        
        # Create done file
        touch {output.done}
        
        echo "Binette run completed successfully" >> {log}
        """


# Marker to assert binette/final_bins/ contains FASTAs (no directory outputs created)
rule binette_final_marker:
    input:
        "binning/binette/.done"
    output:
        flag = "binning/binette/final_bins/.done"
    log:
        f"logs/binette/{PROJECT}.final_marker.log"
    shell:
        r"""
        set -euo pipefail
        dir="binning/binette/final_bins"
        mkdir -p "$dir"
        shopt -s nullglob
        files=( "$dir"/*.fa "$dir"/*.fna "$dir"/*.fasta )
        if [ "${{#files[@]}}" -eq 0 ]; then
            echo "No FASTA bins present in $dir" | tee {log}
            exit 2
        fi
        touch {output.flag}
        """

DREP_SA = "0.995"

rule drep_binette:
    input:
        "binning/binette/final_bins/.done"
    output:
        done = "binning/dRep/binette/.done"
    params:
        outdir = "binning/dRep/binette",
        sa     = DREP_SA,
        extra  = ""
    log:
        f"logs/drep/{PROJECT}_binette.log"
    shell:
        r"""
        set -euo pipefail
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml dRep
        mkdir -p {params.outdir} logs/drep

        BASE="binning/binette/final_bins"
        if ! compgen -G "$BASE/*.fa" >/dev/null && \
           ! compgen -G "$BASE/*.fna" >/dev/null && \
           ! compgen -G "$BASE/*.fasta" >/dev/null; then
            echo "No FASTA files found in $BASE" >&2
            exit 2
        fi

        mapfile -t GENOMES < <(find "$BASE" -maxdepth 1 -type f \( -name '*.fa' -o -name '*.fna' -o -name '*.fasta' \) -print)
        dRep dereplicate {params.outdir} -g "${{GENOMES[@]}}" -sa {params.sa} -p {resources.cpus_per_task} {params.extra} &> {log}

        touch {output.done}
        """

rule checkm2_predict:
    input:
        drep_done = "binning/dRep/binette/.done",
        db = "/lisc/data/scratch/mirror/checkm2/1.1.0/CheckM2_database/uniref100.KO.1.dmnd"
    output:
        predicted_quality = f"binning/checkm2/{PROJECT}/quality_report.tsv"
    log:
        f"logs/binning/checkm2_predict/{PROJECT}.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml CheckM2
        
        BINS_DIR="binning/dRep/binette/dereplicated_genomes"
        if [ ! -d "$BINS_DIR" ]; then
            echo "Expected dRep output folder not found: $BINS_DIR" >&2
            exit 3
        fi
        checkm2 predict \
            -x fa \
            --input "$BINS_DIR" \
            --output-directory binning/checkm2/{PROJECT} \
            --threads {resources.cpus_per_task} \
            --tmpdir $TMPDIR \
            --force \
            --database_path {input.db}
        ls -lh {output.predicted_quality}
        """

#######################################################################
# Anvio
#######################################################################
rule fix_assembly_headers_copy:
    input:
        src = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fasta"
    output:
        dst = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fixed.fa"
    log:
        f"logs/anvio/{PROJECT}.clean_headers.log"
    resources:
        mem_mb = 6000,
        cpus_per_task = 2,
        runtime = 20
    run:
        import os, re, gzip
        src = input.src
        dst = output.dst
        def clean_header(h):
            h = h[1:]
            m = re.match(r'^([A-Za-z0-9_\-]+)', h)
            if m:
                return '>' + m.group(1)
            return '>' + re.split(r'[ \t|]', h)[0]
        def process_fasta(infile, outfile):
            op_in = gzip.open if infile.endswith('.gz') else open
            op_out = gzip.open if outfile.endswith('.gz') else open
            n = 0
            os.makedirs(os.path.dirname(outfile), exist_ok=True)
            with op_in(infile, 'rt') as fin, op_out(outfile, 'wt') as fout:
                for line in fin:
                    if line.startswith('>'):
                        fout.write(clean_header(line.strip()) + '\n')
                        n += 1
                    else:
                        fout.write(line)
            return n
        n = process_fasta(src, dst)
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        with open(log[0], "w") as lg:
            lg.write(f"Cleaned {n} headers into: {dst}\n")

rule build_external_gene_calls:
    input:
        gff = f"genes/prodigal/{PROJECT}_genes.gff",
        faa = f"genes/prodigal/{PROJECT}_genes.faa"
    output:
        tsv = f"anvio/contigs/{PROJECT}_external_gene_calls.tsv"
    params:
        version = "v2.6.3"
    log:
        f"logs/anvio/{PROJECT}.build_external_gene_calls.log"
    resources:
        mem_mb = 16000,
        cpus_per_task = 1,
        runtime = 60
    run:
        import os, re
        gff_file = input.gff
        faa_file = input.faa
        out_tsv  = output.tsv
        version  = params.version
        def simplify(s):
            return re.split(r'[ \t|#]', s)[0]
        seqs = {}
        with open(faa_file) as faa:
            gene_id = None
            seq_lines = []
            for line in faa:
                if line.startswith('>'):
                    if gene_id is not None:
                        seqs[gene_id] = ''.join(seq_lines)
                    gene_id = simplify(line.strip()[1:])
                    seq_lines = []
                else:
                    seq_lines.append(line.strip())
            if gene_id is not None:
                seqs[gene_id] = ''.join(seq_lines)
        os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
        with open(out_tsv, "w") as out_f, open(log[0], "w") as lg:
            out_f.write("gene_callers_id\tcontig\tstart\tstop\tdirection\tpartial\tcall_type\tsource\tversion\taa_sequence\n")
            gene_counter = 1
            kept = 0
            with open(gff_file) as gff:
                for ln in gff:
                    if ln.startswith('#'):
                        continue
                    parts = ln.rstrip('\n').split('\t')
                    if len(parts) != 9:
                        continue
                    contig, source, feature, start, end, score, strand, phase, attributes = parts
                    if feature != 'CDS':
                        continue
                    contig_simple = simplify(contig)
                    attrdict = {}
                    for entry in attributes.split(';'):
                        if '=' in entry:
                            k, v = entry.split('=', 1)
                            attrdict[k] = v
                    gene_id = attrdict.get('ID', '')
                    partial = 1 if attrdict.get('partial', '00') != '00' else 0
                    direction = 'f' if strand == '+' else 'r'
                    start0 = int(start) - 1
                    stop1 = int(end)
                    call_type = 1
                    gene_number = gene_id.split('_')[-1] if '_' in gene_id else gene_id
                    contig_gene_id = f"{contig_simple}_{gene_number}"
                    aa_sequence = seqs.get(contig_gene_id, '')
                    out_f.write(
                        f"{gene_counter}\t{contig_simple}\t{start0}\t{stop1}\t{direction}\t{partial}\t{call_type}\tpprodigal\t{version}\t{aa_sequence}\n"
                    )
                    gene_counter += 1
                    kept += 1
            lg.write(f"Wrote {kept} external gene calls to {out_tsv}\n")

rule anvio_gen_contigs_db:
    input:
        fasta = f"assembly/coassembly/{PROJECT}_coassembly_29_127.fixed.fa",
        egc   = f"anvio/contigs/{PROJECT}_external_gene_calls.tsv"
    output:
        db = "anvio/contigs/contigs.db"
    params:
        label = f"{PROJECT}_contigs"
    log:
        "logs/anvio/gen_contigs_db.log"
    shell:
        r"""
        
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    ml Conda
        
        conda activate anvio-8
        mkdir -p anvio/contigs logs/anvio
        anvi-gen-contigs-database \
          -f {input.fasta} \
          -o {output.db} \
          -n "{params.label}" \
          --external-gene-calls {input.egc} \
          --skip-predict-frame --ignore-internal-stop-codons \
          -T {resources.cpus_per_task} \
          &> {log}
        """

rule anvio_profile_single:
    input:
        bam       = lambda wc: os.path.join("mapping", "bowtie2", "coassembly", "bam", f"{wc.sample}.bam"),
        contigsdb = "anvio/contigs/contigs.db"
    output:
        "anvio/profile/single/{sample}/PROFILE.db"
    params:
        outdir = lambda wc: f"anvio/profile/single/{wc.sample}"
    log:
        "logs/anvio/profile_single/{sample}.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda
        conda activate anvio-8
        mkdir -p {params.outdir} $(dirname {log})
        anvi-profile \
          -i {input.bam} \
          -c {input.contigsdb} \
          -o {params.outdir} \
          -S pr_{PROJECT}_{wildcards.sample} \
          --min-contig-length 1000 \
          --force-overwrite \
          -T {resources.cpus_per_task} \
          &> {log}
        """

rule anvio_profile_merge:
    input:
        db = "anvio/contigs/contigs.db",
        profiles = expand("anvio/profile/single/{sample}/PROFILE.db", sample=SAMPLES)
    output:
        merged = "anvio/profile/merged/PROFILE.db"
    log:
        "logs/anvio/merge_profiles.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda
        conda activate anvio-8
        mkdir -p anvio/profile/merged
        anvi-merge {input.profiles} -c {input.db} -o anvio/profile/merged --force-overwrite &> {log}
        """

rule anvio_copy_rename_bins:
    input:
        drep_done = "binning/dRep/binette/.done"
    output:
        flag = "anvio/bins/non_refined/.done"
    run:
        import os, glob, shutil, pathlib
        src_dir = pathlib.Path("binning/dRep/binette/dereplicated_genomes")
        dst_dir = pathlib.Path("anvio/bins/non_refined")
        dst_dir.mkdir(parents=True, exist_ok=True)
        for old in dst_dir.glob("Bin_*.fa*"):
            try: old.unlink()
            except FileNotFoundError: pass
        files = sorted([f for f in src_dir.glob("*") if f.suffix in [".fa",".fna",".fasta"]])
        for i, src in enumerate(files, start=1):
            shutil.copy(src, dst_dir / f"Bin_{i}.fa")
        pathlib.Path(output.flag).touch()

rule anvio_bins_to_import:
    input:
        flag = "anvio/bins/non_refined/.done"
    output:
        table = "anvio/bins/bins_to_import.txt"
    log:
        "logs/anvio/bins_to_import.log"
    run:
        import os, glob
        os.makedirs(os.path.dirname(output.table), exist_ok=True)
        bin_dir = "anvio/bins/non_refined"
        written = 0
        with open(output.table, "w") as out_f:
            for filepath in sorted(glob.glob(os.path.join(bin_dir, "Bin_*.fa"))):
                bin_name = os.path.splitext(os.path.basename(filepath))[0]
                with open(filepath, "r") as fa:
                    for line in fa:
                        if line.startswith(">"):
                            contig_id = line[1:].strip().split()[0]
                            out_f.write(f"{contig_id}\t{bin_name}\n")
                            written += 1
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        with open(log[0], "w") as lg:
            lg.write(f"Wrote {written} contig-to-bin mappings to {output.table}\n")

rule anvio_import_collection:
    input:
        merged = "anvio/profile/merged/PROFILE.db",
        db     = "anvio/contigs/contigs.db",
        table  = "anvio/bins/bins_to_import.txt"
    output:
        flag = "anvio/profile/merged/collection_OriginalBins.done"
    log:
        "logs/anvio/import_collection.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda
        conda activate anvio-8
        anvi-import-collection \
          -p {input.merged} \
          -c {input.db} \
          -C Binette \
          --contigs-mode {input.table} \
          &> {log}
          
          anvi-import-collection \
          -p {input.merged} \
          -c {input.db} \
          -C Binette_backup \
          --contigs-mode {input.table} \
          &> {log}
        touch {output.flag}
        """

rule anvio_run_hmms:
    input:
        db = "anvio/contigs/contigs.db"
    output:
        flag = "anvio/contigs/hmms.done"
    log:
        "logs/anvio/run_hmms.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda
        conda activate anvio-8
        anvi-run-hmms -c {input.db} \
        --just-do-it \
        -T {resources.cpus_per_task} &> {log}
        touch {output.flag}
        """

rule anvio_run_scg_taxonomy:
    input:
        hmms = "anvio/contigs/hmms.done",
        db   = "anvio/contigs/contigs.db"
    output:
        flag = "anvio/contigs/scg_taxonomy.done"
    log:
        "logs/anvio/run_scg_taxonomy.log"
    shell:
        r"""
        
        # Remove path pollution from loaded modules
    module purge || true
    
    # Sanitize Python import behavior
    unset PYTHONPATH
    export PYTHONNOUSERSITE=1
    
    
        ml Conda
        conda activate anvio-8
        anvi-run-scg-taxonomy -c {input.db} -T {resources.cpus_per_task} &> {log}
        touch {output.flag}
        """
