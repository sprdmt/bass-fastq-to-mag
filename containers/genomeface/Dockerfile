# syntax=docker/dockerfile:1
# bass-fastq-to-mag-genomeface: GenomeFace 1.0.1 (markersgf, genomeface) with TensorFlow 2.9.1 GPU, CUDA 11.7, RAPIDS 22.12
# Build from the repository root (the pin file lives in workflow/envs/):
#   docker build -f containers/genomeface/Dockerfile -t ghcr.io/sprdmt/bass-fastq-to-mag-genomeface:1.0.1 .
# The env is installed from the explicit pin file (fully solved, md5-checked package URLs) generated
# from workflow/envs/genomeface.yaml, so the image matches the conda env used with --sdm conda.
FROM ghcr.io/mamba-org/micromamba:2.9.0

ARG TOOL_VERSION=1.0.1
ARG CONDA_OVERRIDE_CUDA=11.7
# Only used if the YAML is solved instead of the pin file (GPU-less build host).
ENV CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA}

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/genomeface.yaml /opt/env/genomeface.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/genomeface.linux-64.pin.txt /opt/env/genomeface.linux-64.pin.txt
RUN micromamba install -y -n base -f /opt/env/genomeface.linux-64.pin.txt \
 && micromamba clean --all --yes

# Apptainer `exec` bypasses the micromamba entrypoint, so the env must be on PATH via ENV.
ENV PATH=/opt/conda/bin:${PATH} \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1 \
CUDA_PATH=/opt/conda

# Build-time smoke test (CPU only): `genomeface -h` imports CuPy and needs the GPU driver (libcuda),
# so it is the GPU smoke test (apptainer exec --nv) rather than a build step.
RUN markersgf --help > /dev/null && python -c "import tensorflow as tf; print(tf.__version__)"

LABEL org.opencontainers.image.title="bass-fastq-to-mag-genomeface" \
      org.opencontainers.image.description="GenomeFace 1.0.1 (markersgf, genomeface) with TensorFlow 2.9.1 GPU, CUDA 11.7, RAPIDS 22.12" \
      org.opencontainers.image.version="${TOOL_VERSION}" \
      org.opencontainers.image.source="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.url="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.documentation="https://github.com/sprdmt/bass-fastq-to-mag/blob/main/containers/README.md" \
      org.opencontainers.image.base.name="ghcr.io/mamba-org/micromamba:2.9.0"
