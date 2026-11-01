# syntax=docker/dockerfile:1
# bass-fastq-to-mag-anvio: anvi'o 9 (anvio-minimal 9 + HMMER, DIAMOND, Prodigal, SAMtools; SCG taxonomy databases set up at build time)
# Build from the repository root (the pin file lives in workflow/envs/):
#   docker build -f containers/anvio/Dockerfile -t ghcr.io/sprdmt/bass-fastq-to-mag-anvio:9 .
# The env is installed from the explicit pin file (fully solved, md5-checked package URLs) generated
# from workflow/envs/anvio.yaml, so the image matches the conda env used with --sdm conda.
FROM ghcr.io/mamba-org/micromamba:2.9.0

ARG TOOL_VERSION=9
ARG CONDA_OVERRIDE_CUDA=

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/anvio.yaml /opt/env/anvio.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/anvio.linux-64.pin.txt /opt/env/anvio.linux-64.pin.txt
RUN micromamba install -y -n base -f /opt/env/anvio.linux-64.pin.txt \
 && micromamba clean --all --yes

# Apptainer `exec` bypasses the micromamba entrypoint, so the env must be on PATH via ENV.
ENV PATH=/opt/conda/bin:${PATH} \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1

# anvi'o writes the SCG-taxonomy DIAMOND databases into its package dir, which is read-only in a SIF
RUN anvi-setup-scg-taxonomy -T 4

# Build-time smoke test (CPU only): fails the build if the tools cannot start.
RUN anvi-gen-contigs-database --version && anvi-run-scg-taxonomy --help > /dev/null && hmmsearch -h > /dev/null && diamond version

LABEL org.opencontainers.image.title="bass-fastq-to-mag-anvio" \
      org.opencontainers.image.description="anvi'o 9 (anvio-minimal 9 + HMMER, DIAMOND, Prodigal, SAMtools; SCG taxonomy databases set up at build time)" \
      org.opencontainers.image.version="${TOOL_VERSION}" \
      org.opencontainers.image.source="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.url="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.documentation="https://github.com/sprdmt/bass-fastq-to-mag/blob/main/containers/README.md" \
      org.opencontainers.image.base.name="ghcr.io/mamba-org/micromamba:2.9.0"
