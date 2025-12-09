BASS fastq-to-MAG Snakemake Workflow (SLURM profile)

Overview
--------
Snakemake workflow for processing marine shotgun metagenomes from paired-end FASTQ files through coassembly, graph construction, mapping, binning, graph-based refinement, dereplication, quality control, and optional anvi’o integration (fastq → MAGs). The repository includes a Snakemake 9 workflow profile designed for SLURM-based HPC systems. Portability to other schedulers is possible by adapting the profile and software stack.

Status / Scope
--------------
- Target environment: SLURM-managed HPC with a module system and/or conda.
- Execution is orchestrated through a Snakemake “workflow profile”.
- IMPORTANT: A temporarily patched GraphBin2 environment (“graphbin2-patched”) is required to avoid a known runtime bug; this note will be removed once the upstream fix is available.

Repository Layout
-----------------
.
├── README.md
├── Snakefile
├── snakefiles/
│   └── fastq_to_bin.snakefile
├── workflow-profile/
│   └── config.yaml
├── snakemake-config/
│   └── project_config.yaml.example
└── .gitignore

Input Expectations
------------------
Project working directory (example):
  project_root/
    fastq/
      <SAMPLE>_1.fastq.gz
      <SAMPLE>_2.fastq.gz
    snakemake-config/
      project_config.yaml
    logs/

FASTQ naming (paired-end):
  <SAMPLE>_1.fastq.gz
  <SAMPLE>_2.fastq.gz

Minimum project configuration:
  snakemake-config/project_config.yaml
    project: YOUR_PROJECT_ID
  (Optional) fastq_dir: fastq

Running the Workflow
--------------------
1) Point to this repository:
   export SNAKEMAKE_WF_FOLDER=/path/to/this/repo

2) Dry-run (DAG and rule validation):
   snakemake \
     --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" \
     --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" \
     -n

3) Submit jobs with the profile (adjust concurrency to your system):
   snakemake \
     --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" \
     --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" \
     -j 50 -c 400 .

Software Requirements
---------------------
The workflow uses environment modules and conda environments. Names below reflect what the Snakefile expects and can be adapted to local conventions.

Modules (loaded via `ml <Name>` in rules):
- MEGAHIT
- Bowtie2
- SAMtools
- MMseqs2
- SeqKit
- dRep
- CheckM2
- vamb
- Python

Conda environments (activated via `conda activate <env>`):
- agtools                   (FASTG → GFA conversion)
- pprodigal                 (parallel prodigal wrapper)
- CoverM-0.7.0              (contig coverage)
- genomeface                (marker search and genomeface binning)
- SemiBin-2.2.0             (SemiBin2; GPU supported with `--engine gpu`)
- taxconv                   (taxonomy table conversion)
- graphbin2-patched         (patched GraphBin2 runner and utilities)
- binette                   (bin consolidation and QC with CheckM2)
- anvio-8                   (anvi’o database/profile steps)

Optional GPU usage:
- SemiBin2: `--engine gpu` when GPUs are available.
- VAMB: `vamb bin taxvamb --cuda` when GPUs are available.
- MMseqs2: `--gpu` optional per local build and hardware.

Patched GraphBin2
-----------------
- A patched GraphBin2 environment is required to prevent a NameError/logging crash under certain inputs.
- The workflow includes input sanitation and defensive shims, but a patched installation remains necessary until the upstream release resolves the issue.

Privacy & Paths
---------------
- Examples use placeholders (e.g., `/path/to/this/repo`, `YOUR_PROJECT_ID`).

Contributing
------------
Contributions are welcome via pull requests and issues. For substantial changes (new rules, containerization, scheduler backends), open an issue to discuss design and interoperability. Please keep changes modular and update the profile/config documentation when behavior changes.

Issues
------
Use GitHub Issues for:
- Reproducible failures (attach rule logs, Snakemake version, profile settings, and minimal configs).
- Feature requests (additional binners, container images, portability to other schedulers).
- Environment fixes (module names, conda pinning, GPU flags).

How to Cite
-----------


License
-------
This project is distributed under the MIT License. See the `LICENSE` file for details.

Acknowledgements
----------------
This workflow orchestrates community tools (MEGAHIT, Bowtie2, SAMtools, MMseqs2, Prodigal/pprodigal, SeqKit, SemiBin2, VAMB, GraphBin2, Binette, dRep, CheckM2, anvi’o). Please cite individual tools according to their licenses/publications in addition to citing the workflow release/DOI.

Contact
-------
Please use GitHub Issues for questions, bug reports, and feature requests.

Changelog
---------
- v0.1.0 — Initial public release with SLURM profile, patched GraphBin2 requirement, and documented software stack.
