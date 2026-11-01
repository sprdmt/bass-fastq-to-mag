"""Shared helpers: configuration, sample sheet, software and the size-aware resource model."""

import json
import math
import os
import re
import shlex

import pandas as pd
import yaml
from snakemake.exceptions import WorkflowError
from snakemake.utils import validate

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
validate(config, schema="../schemas/config.schema.yaml")

OUT = config["outdir"].rstrip("/")
PROJECT = config["project"]
LOGS = "logs"
BENCH = "benchmarks"


def _deep_merge(base, override):
    """Recursively merge two dicts; values in `override` win. Returns a new dict."""
    merged = dict(base)
    for key, value in (override or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _load_yaml(name):
    with open(os.path.join(workflow.basedir, name)) as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Sample sheet
# ---------------------------------------------------------------------------
def _load_samples(path):
    df = pd.read_csv(path, sep="\t", dtype=str, comment="#", keep_default_na=False)
    df.columns = [c.strip() for c in df.columns]
    df = df.apply(lambda col: col.str.strip())
    if "assembly_group" not in df.columns:
        df["assembly_group"] = ""
    df.loc[df["assembly_group"] == "", "assembly_group"] = PROJECT
    for optional in ("run_accession", "sample_accession"):
        if optional not in df.columns:
            df[optional] = ""
    validate(df, schema="../schemas/samples.schema.yaml")
    duplicated = df["sample"][df["sample"].duplicated()].tolist()
    if duplicated:
        raise WorkflowError(f"Duplicated sample names in {path}: {duplicated}")
    return df.set_index("sample", drop=False)


SAMPLES_DF = _load_samples(config["samples"])
SAMPLES = SAMPLES_DF["sample"].tolist()
GROUPS = sorted(SAMPLES_DF["assembly_group"].unique())
BINNERS = list(config["binning"]["binners"])
GRAPH_REFINE = bool(config["graphbin2"]["enabled"])


def group_samples(group):
    """Samples that are co-assembled (and co-binned) together."""
    return SAMPLES_DF.loc[SAMPLES_DF["assembly_group"] == group, "sample"].tolist()


def raw_reads(sample):
    return [SAMPLES_DF.loc[sample, "fq1"], SAMPLES_DF.loc[sample, "fq2"]]


def clean_reads(sample):
    """Reads used downstream: trimmed reads when QC is enabled, otherwise the raw files."""
    if config["qc"]["enabled"]:
        return [
            f"{OUT}/qc/trimmed/{sample}_R1.fastq.gz",
            f"{OUT}/qc/trimmed/{sample}_R2.fastq.gz",
        ]
    return raw_reads(sample)


def group_reads(group, mate):
    idx = 0 if mate == 1 else 1
    return [clean_reads(s)[idx] for s in group_samples(group)]


def group_bams(group):
    return [f"{OUT}/mapping/{group}/bam/{s}.bam" for s in group_samples(group)]


def contigs(group):
    return f"{OUT}/assembly/{group}/contigs.fa"


wildcard_constraints:
    sample="|".join(re.escape(s) for s in SAMPLES),
    group="|".join(re.escape(g) for g in GROUPS),
    binner="metabat2|genomeface|semibin2|vamb",
    chunk=r"\d+",


# ---------------------------------------------------------------------------
# Software: one conda env file and one pinned container image per tool
# ---------------------------------------------------------------------------
DEFAULT_CONTAINERS = _load_yaml("containers.yaml")
CONTAINERS = _deep_merge(DEFAULT_CONTAINERS, config.get("containers") or {})


def env(name):
    """Conda environment file for --software-deployment-method conda."""
    return f"../envs/{name}.yaml"


def image(name):
    """Container image for --software-deployment-method apptainer."""
    try:
        return CONTAINERS[name]
    except KeyError:
        raise WorkflowError(f"No container image configured for '{name}'.")


# ---------------------------------------------------------------------------
# Helper scripts (workflow/scripts). workflow.source_path() copies them into a per-run cache
# that Snakemake also makes available inside the containers. The cache directory changes with
# every Snakemake invocation, so the rules reach the scripts through the environment variable
# BFM_SCRIPTS: as a `params` entry the changing path would re-run the jobs on every restart.
# ---------------------------------------------------------------------------
for _script in sorted(os.listdir(os.path.join(workflow.basedir, "scripts"))):
    if _script.endswith(".py"):
        _cached_script = workflow.source_path(f"../scripts/{_script}")
os.environ["BFM_SCRIPTS"] = os.path.dirname(_cached_script)


envvars:
    "BFM_SCRIPTS",


# ---------------------------------------------------------------------------
# GPU handling
# ---------------------------------------------------------------------------
GPU = config["gpu"]
USE_GPU = bool(GPU["enabled"])
if not USE_GPU and "genomeface" in BINNERS:
    raise WorkflowError(
        "GenomeFace needs a GPU: its CuPy and TensorFlow code does not start without the CUDA "
        "driver. With gpu.enabled: false, remove 'genomeface' from binning.binners."
    )


def gpu_resources():
    """Scheduler resources that request one GPU (none when GPUs are disabled)."""
    if not USE_GPU or GPU["request"] == "none":
        return {}
    if GPU["request"] == "gres":
        return {"gres": GPU["gres"]}
    return {"gpu": int(GPU["count"])}


# ---------------------------------------------------------------------------
# Size-aware resource model
#
# Every rule names an entry of workflow/resources.yaml (overridable under
# `resources:` in the config). For each of threads, mem_mb and runtime:
#
#   value = base + per_gb * GB + per_sample * N + per_gb_sample * GB * N + per_file * F
#
# GB  = total size of the inputs listed in `size_of` (all inputs if omitted;
#       directories are walked), N = samples in the assembly group,
#       F = number of input files. The value is raised to `min`, multiplied by
#       retry_factor^(attempt - 1) for mem_mb/runtime, and clamped to the rule
#       `max` and to the global limits of your cluster.
# ---------------------------------------------------------------------------
RESOURCES = _deep_merge(_load_yaml("resources.yaml"), config.get("resources") or {})
LIMITS = RESOURCES["limits"]
_LIMIT_KEY = {
    "threads": "max_threads",
    "mem_mb": "max_mem_mb",
    "runtime": "max_runtime",
}


def _selected_paths(input, names):
    if not names:
        return [str(p) for p in input]
    paths = []
    for name in names:
        value = getattr(input, name)
        if isinstance(value, (list, tuple)):
            paths.extend(str(p) for p in value)
        else:
            paths.append(str(value))
    return paths


def _measure(paths):
    """Return (size in GB, number of files). Missing paths raise FileNotFoundError,
    which Snakemake reports as <TBD> in dry runs."""
    total, nfiles = 0, 0
    for path in paths:
        stat = os.stat(path)
        if os.path.isdir(path):
            for root, _dirs, files in os.walk(path):
                for fname in files:
                    total += os.path.getsize(os.path.join(root, fname))
                    nfiles += 1
        else:
            total += stat.st_size
            nfiles += 1
    return total / 1e9, nfiles


def _n_samples(wildcards):
    if "group" in wildcards.keys():
        return len(group_samples(wildcards.group))
    if "sample" in wildcards.keys():
        return 1
    return len(SAMPLES)


def _spec(entry, key):
    rules = RESOURCES["rules"]
    if entry not in rules or key not in rules[entry]:
        raise WorkflowError(
            f"No '{key}' resource entry for '{entry}' (workflow/resources.yaml)."
        )
    spec = rules[entry][key]
    return {"base": spec} if isinstance(spec, (int, float)) else spec


def estimate(entry, key, wildcards, input, attempt=1, strict=True):
    """Estimate one resource value for a job (see module comment)."""
    spec = _spec(entry, key)
    gb, nfiles = 0.0, 0
    if any(spec.get(k) for k in ("per_gb", "per_gb_sample", "per_file")):
        try:
            gb, nfiles = _measure(_selected_paths(input, spec.get("size_of")))
        except FileNotFoundError:
            if strict:
                raise
    n = _n_samples(wildcards)
    value = (
        spec.get("base", 0)
        + spec.get("per_gb", 0) * gb
        + spec.get("per_sample", 0) * n
        + spec.get("per_gb_sample", 0) * gb * n
        + spec.get("per_file", 0) * nfiles
    )
    value = max(value, spec.get("min", 0))
    if key != "threads" and attempt > 1:
        value *= RESOURCES["retry_factor"] ** (attempt - 1)
    cap = min(spec.get("max", math.inf), LIMITS[_LIMIT_KEY[key]])
    return max(1, int(math.ceil(min(value, cap))))


def threads_for(entry):
    """Thread count; falls back to the base value while inputs do not exist yet."""

    def _threads(wildcards, input, attempt):
        return estimate(entry, "threads", wildcards, input, attempt, strict=False)

    return _threads


def _pick_partition(entry, wildcards, input, attempt, partitions):
    mem = estimate(entry, "mem_mb", wildcards, input, attempt)
    runtime = estimate(entry, "runtime", wildcards, input, attempt)
    threads = estimate(entry, "threads", wildcards, input, attempt, strict=False)
    ordered = sorted(partitions, key=lambda p: p.get("max_mem_mb", math.inf))
    for part in ordered:
        if (
            mem <= part.get("max_mem_mb", math.inf)
            and runtime <= part.get("max_runtime", math.inf)
            and threads <= part.get("max_threads", math.inf)
        ):
            return part["name"]
    return ordered[-1]["name"]


def resources_for(entry, gpu=False):
    """Scheduler resources for a rule: mem_mb, runtime, optional partition and GPU request."""
    res = {
        "mem_mb": lambda wildcards, input, attempt: estimate(
            entry, "mem_mb", wildcards, input, attempt
        ),
        "runtime": lambda wildcards, input, attempt: estimate(
            entry, "runtime", wildcards, input, attempt
        ),
    }
    partitions = RESOURCES.get("gpu_partitions" if gpu else "partitions") or []
    if partitions:
        res["slurm_partition"] = lambda wildcards, input, attempt: _pick_partition(
            entry, wildcards, input, attempt, partitions
        )
    if gpu:
        res.update(gpu_resources())
    return res


def _mem_mb(resources):
    """Job memory in MB, or None while it is still undetermined (<TBD> in dry runs)."""
    try:
        return float(str(resources.mem_mb))
    except (TypeError, ValueError):
        return None


def mem_bytes(resources, fraction=0.9):
    """Memory a tool may use, in bytes (a fraction of the job allocation)."""
    mb = _mem_mb(resources)
    return "<TBD>" if mb is None else int(mb * 1_000_000 * fraction)


def per_thread_mb(resources, threads, fraction=0.5):
    """Memory per thread in MB, e.g. for samtools sort -m."""
    mb = _mem_mb(resources)
    return "<TBD>" if mb is None else max(256, int(mb * fraction / max(1, threads)))


# ---------------------------------------------------------------------------
# Rule helpers (input functions and parameter builders)
# ---------------------------------------------------------------------------
def cutadapt_args():
    q = config["qc"]
    fwd, rev = ("-b", "-B") if q["adapter_mode"] == "anywhere" else ("-a", "-A")
    args = [f"{fwd} {q['adapter_r1']}", f"{rev} {q['adapter_r2']}"]
    if q["nextseq_trim"] is not None:
        args.append(f"--nextseq-trim={q['nextseq_trim']}")
    if q["quality_cutoff"]:
        args.append(f"-q {q['quality_cutoff']}")
    args += [f"--overlap {q['overlap']}", f"--pair-filter={q['pair_filter']}"]
    if q["cut_r1"]:
        args.append(f"-u {q['cut_r1']}")
    if q["cut_r2"]:
        args.append(f"-U {q['cut_r2']}")
    args.append(f"-m {q['minimum_length']}")
    if q["extra"]:
        args.append(q["extra"])
    return " ".join(args)


def multiqc_inputs(wildcards):
    files = []
    if config["qc"]["enabled"]:
        files += expand(f"{OUT}/qc/cutadapt/{{sample}}.cutadapt.txt", sample=SAMPLES)
    if config["qc"]["fastqc"]:
        stages = ["raw", "trimmed"] if config["qc"]["enabled"] else ["raw"]
        files += expand(
            f"{OUT}/qc/fastqc/{{stage}}/{{sample}}_R{{mate}}_fastqc.zip",
            stage=stages,
            sample=SAMPLES,
            mate=[1, 2],
        )
    return files


def fastqc_reads(wc):
    return raw_reads(wc.sample) if wc.stage == "raw" else clean_reads(wc.sample)


def megahit_k_params():
    a = config["assembly"]
    if a["k_list"]:
        return f"--k-list {a['k_list']}"
    return f"--k-min {a['k_min']} --k-max {a['k_max']} --k-step {a['k_step']}"


def group_member_reads(wc):
    if SAMPLES_DF.loc[wc.sample, "assembly_group"] != wc.group:
        raise WorkflowError(f"Sample {wc.sample} is not in assembly group {wc.group}.")
    return clean_reads(wc.sample)


def binner_rundir(group, binner):
    return f"{OUT}/binning/{group}/{binner}/run"


def binner_binsdir(group, binner):
    """Where each binner writes its bin FASTA files inside its run directory."""
    sub = {"metabat2": "", "genomeface": "", "semibin2": "output_bins", "vamb": "bins"}[
        binner
    ]
    return (
        os.path.join(binner_rundir(group, binner), sub)
        if sub
        else binner_rundir(group, binner)
    )


def refined_table(group, binner):
    if GRAPH_REFINE:
        return f"{OUT}/graphbin2/{group}/{binner}/contig2bin.tsv"
    return f"{OUT}/binning/{group}/{binner}/contig2bin.tsv"


def binette_tables(group):
    """(label, contig2bin table) pairs given to Binette for one group.

    graphbin2.binette_input: both = the original and the GraphBin2 bins of every binner (Binette
    keeps the better version of each genome); refined = GraphBin2 bins only.
    """
    original = [(b, f"{OUT}/binning/{group}/{b}/contig2bin.tsv") for b in BINNERS]
    if not GRAPH_REFINE:
        return original
    refined = [
        (f"{b}_graphbin2", f"{OUT}/graphbin2/{group}/{b}/contig2bin.tsv")
        for b in BINNERS
    ]
    if config["graphbin2"]["binette_input"] == "refined":
        return refined
    return original + refined


def candidate_inputs(wildcards):
    """Bins and quality reports that become MAG candidates."""
    if config["mags"]["source"] == "anvio":
        return {
            "bins": expand(f"{OUT}/anvio/{{group}}/candidates", group=GROUPS),
            "reports": expand(
                f"{OUT}/anvio/{{group}}/candidates_quality.tsv", group=GROUPS
            ),
        }
    return {
        "bins": expand(f"{OUT}/binette/{{group}}/final_bins", group=GROUPS),
        "reports": expand(
            f"{OUT}/binette/{{group}}/final_bins_quality_reports.tsv", group=GROUPS
        ),
    }


def mag_domains_input(wildcards):
    return (
        f"{OUT}/mags/gtdbtk/gtdbtk.summary.tsv" if config["gtdbtk"]["enabled"] else []
    )


def anvio_name(text):
    """anvi'o names: letters, digits and underscores, starting with a letter."""
    name = re.sub(r"[^A-Za-z0-9_]", "_", text)
    return name if name[:1].isalpha() else f"s_{name}"


# ---------------------------------------------------------------------------
# Target lists
# ---------------------------------------------------------------------------
def qc_targets():
    if config["qc"]["enabled"] or config["qc"]["fastqc"]:
        return [f"{OUT}/reports/multiqc_reads.html"]
    return []


def assembly_targets():
    targets = expand(
        f"{OUT}/assembly/{{group}}/{{f}}",
        group=GROUPS,
        f=["contigs.fa", "assembly_graph.gfa", "assembly_graph_stats.tsv"],
    )
    if config["assembly"]["quast"]:
        targets += expand(f"{OUT}/assembly/{{group}}/quast/report.tsv", group=GROUPS)
    return targets


def bins_targets():
    return expand(
        f"{OUT}/binette/{{group}}/final_bins_quality_reports.tsv", group=GROUPS
    )


def anvio_targets():
    return expand(f"{OUT}/anvio/{{group}}/collections.done", group=GROUPS)


def mags_targets():
    return [f"{OUT}/mags/mag_summary.tsv"]


def mimag_rna_args():
    rna = config["mags"]["rna"]
    args = [f"--trna-mode {rna['trna_mode']}", f"--min-trna {rna['min_trna']}"]
    if rna["exclude_partial_rrna"]:
        args.append("--exclude-partial-rrna")
    return " ".join(args)


# ---------------------------------------------------------------------------
# ENA submission helpers
# ---------------------------------------------------------------------------
def _tag_version(uri):
    """Version from an image tag (BioContainers 'tool:2.17--h1234_0' -> '2.17'); '' if no tag."""
    last = os.path.basename(str(uri or ""))
    tag = last.rsplit(":", 1)[1] if ":" in last else ""
    return tag.split("--")[0]


def tool_version(name):
    """Tool version for the ENA metadata, from the tag of the image used for `name`.

    A local image file (for example a .sif) has no tag; then the version of the workflow's
    default image is reported. Set the ena.*_software fields if your local images differ.
    """
    return (
        _tag_version(image(name))
        or _tag_version(DEFAULT_CONTAINERS.get(name))
        or "unknown"
    )


BINNER_LABELS = {
    "metabat2": ("MetaBAT2", "metabat2"),
    "genomeface": ("GenomeFace", "genomeface"),
    "semibin2": ("SemiBin2", "semibin2"),
    "vamb": ("TaxVAMB", "vamb"),
}


def ena_effective_config():
    """The ena: block with empty software fields filled from the workflow's images."""
    ena = json.loads(json.dumps(config["ena"]))
    if not ena.get("assembly_software"):
        ena["assembly_software"] = f"MEGAHIT_v{tool_version('megahit')}"
    if not ena.get("stats_generation_software"):
        ena["stats_generation_software"] = f"CheckM2_v{tool_version('checkm2')}"
    if not ena.get("binning_software"):
        tools = [
            f"{BINNER_LABELS[b][0]}_v{tool_version(BINNER_LABELS[b][1])}"
            for b in BINNERS
        ]
        if GRAPH_REFINE:
            tools.append(f"GraphBin2_v{tool_version('graphbin2')}-bfm")
        tools.append(f"Binette_v{tool_version('binette')}")
        ena["binning_software"] = ",".join(tools)
    if not ena["coassembly"].get("program"):
        ena["coassembly"]["program"] = f"MEGAHIT v{tool_version('megahit')}"
    return ena


# parameters the ENA genome table needs for every MAG (workflow/scripts/ena_genome_table.py)
ENA_PARAM_KEYS = [
    "metagenome",
    "broad_environment",
    "local_environment",
    "environmental_medium",
    "assembly_software",
    "binning_software",
    "binning_parameters",
    "stats_generation_software",
]


def ena_missing_params():
    """ENA parameters that are empty for at least one assembly group (with ena.groups overrides)."""
    ena = ena_effective_config()
    overrides = ena.get("groups") or {}
    missing = {}
    for group in GROUPS:
        for key in ENA_PARAM_KEYS:
            value = (overrides.get(group) or {}).get(key) or ena.get(key)
            if value is None or not str(value).strip() or "_vunknown" in str(value):
                missing.setdefault(key, []).append(group)
    return [f"ena.{k} (groups: {', '.join(g)})" for k, g in missing.items()]


def require_ena_params():
    """Stop ENA preparation at the start of a run if the ena: block is incomplete."""
    missing = ena_missing_params()
    if missing:
        raise WorkflowError(
            "Empty ENA parameters or unknown tool versions, set them in the ena: block of "
            "the config: " + "; ".join(missing)
        )


def gtdb_metadata_inputs(wildcards):
    meta = config["ena"]["gtdb_metadata"]
    return [p for p in (meta["bac120"], meta["ar53"]) if p]


def ena_lineage_inputs(wildcards):
    ena = config["ena"]
    method = ena["lineage_method"]
    inputs = {"summary": f"{OUT}/mags/gtdbtk/gtdbtk.summary.tsv"}
    if method in ("majority_vote", "consensus"):
        inputs["mv"] = f"{OUT}/ena/gtdbtk_majority_vote.tsv"
    if method in ("mapping_table", "consensus") or ena["lineage_cross_check"]:
        inputs["workbooks"] = ena["gtdb_ncbi_workbooks"]
    if ena["taxdump"]:
        inputs["taxdump"] = ena["taxdump"]
    return inputs


def ena_lineage_args(wildcards, input):
    ena = config["ena"]
    args = [
        f"--method {ena['lineage_method']}",
        f"--max-rank {ena['lineage_max_rank']}",
    ]
    if "mv" in input.keys():
        args.append(f"--mv-table {input.mv}")
    if ena["lineage_cross_check"] and ena["lineage_method"] == "majority_vote":
        args.append("--cross-check")
    if "workbooks" in input.keys() and len(input.workbooks):
        args.append("--mapping-table " + " ".join(input.workbooks))
    if "taxdump" in input.keys():
        args.append(f"--taxdump {input.taxdump}")
    return " ".join(args)


def ena_coassembly_groups():
    """Groups whose primary assembly needs a registered co-assembly sample."""
    if not config["ena"]["coassembly"]["submit_primary"]:
        return []
    return [g for g in GROUPS if len(group_samples(g)) > 1]


def ena_prepare_targets():
    return [f"{OUT}/ena/genome_uploader_input.tsv"] + expand(
        f"{OUT}/ena/coassembly/{{group}}/sample", group=ena_coassembly_groups()
    )


def ena_submit_inputs(wildcards):
    inputs = {
        "table": f"{OUT}/ena/genome_uploader_input.tsv",
        "validation": f"{OUT}/ena/validation/validate_table.txt",
    }
    for g in ena_coassembly_groups():
        inputs[f"xml_{g}"] = f"{OUT}/ena/coassembly/{g}/sample"
        inputs[f"contigs_{g}"] = f"{OUT}/assembly/{g}/contigs.fa"
        inputs[f"coverm_{g}"] = f"{OUT}/mapping/{g}/coverm_contig_mean.tsv"
    return inputs


def ena_coassembly_commands(wildcards, input):
    """Shell lines that submit each co-assembly to the ENA TEST server."""
    ena = ena_effective_config()
    lines = []
    for g in ena_coassembly_groups():
        cmd = [
            "--xml-dir",
            getattr(input, f"xml_{g}"),
            "--samples",
            config["samples"],
            "--group",
            g,
            "--study",
            ena["study"],
            "--assembly-fasta",
            getattr(input, f"contigs_{g}"),
            "--assembly-name",
            f"{PROJECT}_{g}_coassembly",
            "--program",
            ena["coassembly"]["program"],
            "--platform",
            ena["coassembly"]["platform"],
            "--coverm",
            getattr(input, f"coverm_{g}"),
            "--coverage-method",
            config["mags"]["coverage_method"],
            "--outdir",
            f"{OUT}/ena/coassembly/{g}/manifest",
            "--workdir",
            f"{OUT}/ena/submit",
        ]
        args = " ".join(shlex.quote(str(x)) for x in cmd)
        lines.append(
            f"python $BFM_SCRIPTS/ena_submit.py coassembly {args}"
            f" >>{LOGS}/ena/submit_test.log 2>&1"
        )
    return "\n".join(lines) if lines else "true"


def default_targets():
    targets = qc_targets() + assembly_targets() + bins_targets()
    if config["anvio"]["enabled"] or config["mags"]["source"] == "anvio":
        targets += anvio_targets()
    targets += mags_targets()
    if config["ena"]["enabled"]:
        targets += ena_prepare_targets()
    return targets
