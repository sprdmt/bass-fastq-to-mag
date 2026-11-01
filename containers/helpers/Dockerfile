# syntax=docker/dockerfile:1
# bass-fastq-to-mag-helpers: Python 3.12 with pandas, biopython, pyyaml and openpyxl for workflow/scripts
# Build from the repository root (the pin file lives in workflow/envs/):
#   docker build -f containers/helpers/Dockerfile -t ghcr.io/sprdmt/bass-fastq-to-mag-helpers:1.0.0 .
# The env is installed from the explicit pin file (fully solved, md5-checked package URLs) generated
# from workflow/envs/helpers.yaml, so the image matches the conda env used with --sdm conda.
FROM ghcr.io/mamba-org/micromamba:2.9.0

ARG TOOL_VERSION=1.0.0
ARG CONDA_OVERRIDE_CUDA=

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/helpers.yaml /opt/env/helpers.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/helpers.linux-64.pin.txt /opt/env/helpers.linux-64.pin.txt
RUN micromamba install -y -n base -f /opt/env/helpers.linux-64.pin.txt \
 && micromamba clean --all --yes

# Apptainer `exec` bypasses the micromamba entrypoint, so the env must be on PATH via ENV.
ENV PATH=/opt/conda/bin:${PATH} \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1

# Build-time smoke test (CPU only): fails the build if the tools cannot start.
RUN python -c "import pandas, Bio, yaml, openpyxl; print(pandas.__version__, Bio.__version__)"

LABEL org.opencontainers.image.title="bass-fastq-to-mag-helpers" \
      org.opencontainers.image.description="Python 3.12 with pandas, biopython, pyyaml and openpyxl for workflow/scripts" \
      org.opencontainers.image.version="${TOOL_VERSION}" \
      org.opencontainers.image.source="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.url="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.documentation="https://github.com/sprdmt/bass-fastq-to-mag/blob/main/containers/README.md" \
      org.opencontainers.image.base.name="ghcr.io/mamba-org/micromamba:2.9.0"
