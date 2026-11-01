# syntax=docker/dockerfile:1
# bass-fastq-to-mag-vamb-gpu: VAMB/TaxVamb 5.0.4 with CUDA 12.9 PyTorch 2.8.0, plus taxconverter 1.1.0 in a second env
# Build from the repository root (the pin file lives in workflow/envs/):
#   docker build -f containers/vamb-gpu/Dockerfile -t ghcr.io/sprdmt/bass-fastq-to-mag-vamb-gpu:5.0.4 .
# The env is installed from the explicit pin file (fully solved, md5-checked package URLs) generated
# from workflow/envs/vamb-gpu.yaml, so the image matches the conda env used with --sdm conda.
FROM ghcr.io/mamba-org/micromamba:2.9.0

ARG TOOL_VERSION=5.0.4
ARG CONDA_OVERRIDE_CUDA=12.9
# Only used if the YAML is solved instead of the pin file (GPU-less build host).
ENV CONDA_OVERRIDE_CUDA=${CONDA_OVERRIDE_CUDA}

COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/vamb-gpu.yaml /opt/env/vamb-gpu.yaml
COPY --chown=$MAMBA_USER:$MAMBA_USER workflow/envs/vamb-gpu.linux-64.pin.txt /opt/env/vamb-gpu.linux-64.pin.txt
RUN micromamba install -y -n base -f /opt/env/vamb-gpu.linux-64.pin.txt \
 && micromamba clean --all --yes
# taxconverter 1.1.0 pins python<3.12, pandas==2.0.0, numpy==1.26.4: separate env, entry point linked into PATH
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/vamb-gpu/taxconverter.conda.linux-64.pin.txt /opt/env/taxconverter.conda.linux-64.pin.txt
COPY --chown=$MAMBA_USER:$MAMBA_USER containers/vamb-gpu/taxconverter-requirements.txt /opt/env/taxconverter-requirements.txt
RUN micromamba create -y -n taxconverter -f /opt/env/taxconverter.conda.linux-64.pin.txt \
 && /opt/conda/envs/taxconverter/bin/python -m pip install --no-deps --no-cache-dir --require-hashes -r /opt/env/taxconverter-requirements.txt \
 && ln -s /opt/conda/envs/taxconverter/bin/taxconverter /opt/conda/bin/taxconverter \
 && micromamba clean --all --yes

# Apptainer `exec` bypasses the micromamba entrypoint, so the env must be on PATH via ENV.
ENV PATH=/opt/conda/bin:${PATH} \
    MAMBA_ROOT_PREFIX=/opt/conda \
    PYTHONNOUSERSITE=1

# Build-time smoke test (CPU only): fails the build if the tools cannot start.
RUN vamb --version && python -c "import torch; print(torch.__version__, torch.version.cuda)" && taxconverter --help > /dev/null

LABEL org.opencontainers.image.title="bass-fastq-to-mag-vamb-gpu" \
      org.opencontainers.image.description="VAMB/TaxVamb 5.0.4 with CUDA 12.9 PyTorch 2.8.0, plus taxconverter 1.1.0 in a second env" \
      org.opencontainers.image.version="${TOOL_VERSION}" \
      org.opencontainers.image.source="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.url="https://github.com/sprdmt/bass-fastq-to-mag" \
      org.opencontainers.image.documentation="https://github.com/sprdmt/bass-fastq-to-mag/blob/main/containers/README.md" \
      org.opencontainers.image.base.name="ghcr.io/mamba-org/micromamba:2.9.0"
