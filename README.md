# BASS fastq-to-MAG Snakemake Workflow

This repository contains a Snakemake workflow for processing marine shotgun metagenomes from paired-end FASTQ files through coassembly, mapping, binning, graph-based refinement, dereplication, quality control, and optional anvi’o integration (fastq → MAGs).

The workflow is implemented in:
- snakefiles/fastq_to_bin.snakefile — main workflow file.
- workflow-profile/config.yaml — Snakemake 9 workflow profile for LiSC (SLURM).
- snakemake-config/project_config.yaml.example — template for per-project configuration.

The pipeline is currently tailored to the LiSC HPC cluster (University of Vienna) and depends on its module system, SLURM, and shared databases. Porting to other environments is possible but not provided here.

## Platform and dependencies

- Cluster: LiSC HPC
- Scheduler: SLURM
- Workflow engine: snakemake/9.13 (module)
- Execution model: Snakemake 9 workflow profile (--workflow-profile)

## Repository layout

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

## Project directory structure

/project/path/
  fastq/
    SAMPLE_1.fastq.gz
    SAMPLE_2.fastq.gz
  snakemake-config/
    project_config.yaml
  logs/

FASTQs must follow naming:
<SAMPLE>_1.fastq.gz
<SAMPLE>_2.fastq.gz

## Project configuration

mkdir -p snakemake-config
cp /lisc/data/scratch/oceanography/spriahailo/bass-fastq-to-mag/snakemake-config/project_config.yaml.example snakemake-config/project_config.yaml

Edit project_config.yaml:
project: 1_HE327_1_10m

## Running the workflow on LiSC

export SNAKEMAKE_WF_FOLDER=/lisc/data/scratch/oceanography/spriahailo/bass-fastq-to-mag

### SLURM script example (run_snakemake.sh)

#!/bin/bash
#SBATCH --job-name=skmk_bs
#SBATCH --cpus-per-task=1
#SBATCH --mem=1GB
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=dmytro.spriahailo@univie.ac.at
#SBATCH --output=logs/skmk-%j.out
#SBATCH --error=logs/skmk-%j.err

export SNAKEMAKE_WF_FOLDER=/lisc/data/scratch/oceanography/spriahailo/bass-fastq-to-mag
ml snakemake/9.13
cd /project/path

snakemake --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" --unlock
snakemake --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" --touch
snakemake --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" -j 50 -c 400 .

Submit with:
sbatch run_snakemake.sh

## Dry-run mode

snakemake --snakefile "$SNAKEMAKE_WF_FOLDER/snakefiles/fastq_to_bin.snakefile" --workflow-profile "$SNAKEMAKE_WF_FOLDER/workflow-profile/" -n

## Citation

Spriahailo, D. (2025). BASS fastq-to-MAG Snakemake workflow. GitHub: https://github.com/sprdmt/bass-fastq-to-mag
