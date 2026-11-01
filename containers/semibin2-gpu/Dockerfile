# syntax=docker/dockerfile:1
# bass-fastq-to-mag-semibin2-gpu: SemiBin2 2.4.1 with CUDA 12.9 PyTorch 2.8.0
# Build from the repository root (the pin file lives in workflow/envs/):
#   docker build -f containers/semibin2-gpu/Dockerfile -t ghcr.io/sprdmt/bass-fastq-to-mag-semibin2-gpu:2.4.1 .
# The env is installed from the explicit pin file (fully solved, md5-checked package URLs) generated
# from workflow/envs/semibin2-gpu.yaml, so the image matches the conda env used with --sdm conda.
FROM ghcr.io/mamba-org/micromamba:2.9.0

ARG TOOL_VERSION=2.4.1
ARG CONDA_OVERRIDE_CUDA=12.9
# Only used if the YAML is solved instead of the pin file (GPU-less build host).
ENV CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA}

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/semibin2-gpu.yaml /opt/env/semibin2-gpu.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/semibin2-gpu.linux-64.pin.txt /opt/env/semibin2-gpu.linux-64.pin.txt
RUN micromamba install -y -n base -f /opt/env/semibin2-gpu.linux-64.pin.txt \
 && micromamba clean --all --yes

# Apptainer `exec` bypasses the micromamba entrypoint, so the env must be on PATH via ENV.
ENV PATH=/opt/conda/bin:${PATH} \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1

# Build-time smoke test (CPU only): fails the build if the tools cannot start.
RUN SemiBin2 --version && python -c "import torch; print(torch.__version__, torch.version.cuda)"

LABEL org.opencontainers.image.title="bass-fastq-to-mag-semibin2-gpu" \
      org.opencontainers.image.description="SemiBin2 2.4.1 with CUDA 12.9 PyTorch 2.8.0" \
      org.opencontainers.image.version="${TOOL_VERSION}" \
      org.opencontainers.image.source="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.url="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.documentation="https://github.com/sprdmt/bass-fastq-to-mag/blob/main/containers/README.md" \
      org.opencontainers.image.base.name="ghcr.io/mamba-org/micromamba:2.9.0"
