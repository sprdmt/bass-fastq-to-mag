# Containers and conda environments (bass-fastq-to-mag v1.0.0)

Every rule group has one pinned conda environment, `workflow/envs/<env>.yaml`, and one container image.
`containers/manifest.tsv` maps env -> tool/version -> image URI and records how each image was verified on LiSC
(2026-10-31). `config["containers"][<env>]` takes the `image_uri` column.

## Layout

```
workflow/envs/<env>.yaml                 pinned env (channels conda-forge, bioconda, nodefaults; strict priority)
workflow/envs/<env>.linux-64.pin.txt     explicit, md5-checked package list (Snakemake uses it before the YAML)
containers/<image>/Dockerfile            custom images (micromamba base, installs the pin file)
containers/<image>/<image>.def           Apptainer equivalent (same base, same pin file)
containers/<image>/build.env             TOOL_VERSION and CONDA_OVERRIDE_CUDA for CI
containers/manifest.tsv                  env, tool, version, env file, image URI, verification
containers/manifest_alternatives.tsv     alternative images (validated older versions, single-tool ENA images)
.github/workflows/containers.yml         builds/pushes custom images to ghcr.io/sprdmt/bass-fastq-to-mag-<image>:<version>
```

Image sources: `biocontainers` = `quay.io/biocontainers/<pkg>:<version>--<build>` (the build of the pinned conda
package); `mulled` = the multi-package Seqera/Wave image that nf-core/modules `bowtie2/align` pinned from 2025-03-27
to 2026-10-20 (bowtie2 2.5.4, samtools/htslib 1.21, pigz 2.8); `custom` = built here; `official` = upstream image.

Custom images: `genomeface`, `semibin2-gpu`, `vamb-gpu` (+ taxconverter), `anvio` (no official v9 image exists and the
Bioconda `anvio-minimal` image lacks HMMER/DIAMOND/Prodigal/SAMtools) and `helpers` (Python for `workflow/scripts`).
The ENA rules use the BioContainers image of genome-uploader 2.5.1 (with webin_cli_handler and Webin-CLI 9.0.1),
the version the ENA scripts were tested with.

## Version decisions

Policy: start from what the user validated on LiSC; move to a newer release only for a documented fix.

| Tool | User-validated | Pinned | Decision | Reason |
|---|---|---|---|---|
| cutadapt | 5.2 (LiSC module cutadapt/5.2) | 5.2 | kept | latest release |
| FastQC | 0.12.1 (inferred: LiSC FastQC/0.12.1 era of the bass3 script) | 0.12.1 | kept | 0.13.0 (2026-10) moves the launcher to python3 and is documented as not drop-in |
| MultiQC | not used before | 1.35 | new | latest release; LiSC central env MultiQC-1.35 |
| MEGAHIT | 1.2.9 | 1.2.9 | kept | only release |
| agtools | 0.1.2 (env agtools) | 0.1.2 | kept | 1.x fixes duplicate links only in the GFA reader (commit 654ba27), not fastg2gfa; 1.1.0 changes --output from folder to file |
| QUAST | 5.3.0 | 5.3.0 | kept | latest release |
| pprodigal / Prodigal | 1.0.1 / 2.6.3 | 1.0.1 / 2.6.3 | kept | only releases |
| Bowtie2 | 2.5.4 (anvio-9 env; LiSC default module in Jan 2026, inferred) | 2.5.4 | kept | 2.5.5 only adds -d and runtime AVX dispatch |
| SAMtools/HTSlib | 1.21 (anvio-9 env; LiSC SAMtools/1.21 module) | 1.21 | kept | nf-core-tested combination; 1.21.1-1.21.3 security fixes (mpileup ref cache, aux types, cram-size) are not packaged in Bioconda and do not affect view/sort/index/idxstats on workflow-generated BAMs |
| CoverM | 0.7.0 (LiSC CoverM-0.7.0) | 0.7.0 | kept | 0.8.0 adds mappers/features only |
| MetaBAT2 | 2.18 (inferred: LiSC default MetaBAT module) | 2.18 | kept | exact pin (==2.18) avoids the 2.18_23_gc869c52 development snapshot |
| SeqKit | 2.12.0 (inferred: newest LiSC SeqKit module in Jan 2026) | 2.12.0 | kept | 2.12.0 fixes the 2.11.0 split2 --by-part-prefix regression; 2.13/2.14 add options/speed |
| GenomeFace | 1.0.1 (env genomeface) | 1.0.1 | kept (deps fixed) | only release; image = the user's 218-package env except pyrodigal 2.3.0 / pyhmmer 0.10.15: with the user's 3.6.3.post1 / 0.11.1, markersgf fails at pyrodigal.OrfFinder (removed in 3.0; observed on LiSC); it also uses TopHits.query_name, which pyhmmer 0.10.15 reports as deprecated and removed in 0.11.0 |
| SemiBin2 | 2.2.0 (LiSC SemiBin-2.2.0) | 2.4.1 | upgraded | 2.2.1 igraph compatibility (#208); 2.3.0 Prodigal errors silently ignored and stale markers.hmmout reused on reruns; 2.4.0 more fixes. 2.5/2.6 not taken (features, pandas 3) |
| MMseqs2 | LiSC MMseqs2 17/18 modules (default in Jan 2026 unknown) | 18.8cc5c | kept | Bioconda build is CUDA-enabled since 16.747c6; `--gpu 1` search verified on a T4 (20/20 correct top hits) |
| VAMB | 5.0.4 (LiSC vamb-5.0.4) | 5.0.4 | kept | 5.0.4 fixes two critical TaxVamb bugs |
| taxconverter | 1.1.0 (env taxconv) | 1.1.0 | kept | only release (PyPI); needs python<3.12, pandas 2.0.0, numpy 1.26.4 -> separate env |
| GraphBin2 | 1.3.3 (env graphbin2-patched, Python 3.13.5) | 1.3.3 (build pyhdfd78af_1) | kept | corrected script runs inside; build _1 image ships Python 3.12.13 (build _0 image: Python 3.9.20) |
| Binette | 1.1.2 (env binette) | 1.2.1 | upgraded | 1.2.0 no longer crashes on empty DIAMOND results (#72/#76); 1.2.1 fixes a severe contig-filtering slowdown (#89/#90). CLI (Typer) and final TSV columns changed in 1.2.0 |
| dRep | 3.6.2 (LiSC module) | 3.6.2 | kept | 4.x changes default algorithms and results; its skani coverage fix does not apply to the 3.6.2 default |
| CheckM2 | 1.1.0 (module, binette env) | 1.1.0 (build pyh7e72e81_1) | kept | build of the user binette env / LiSC Binette-1.2.1 (TF 2.17) |
| GTDB-Tk | 2.5.2 (LiSC GTDBTk-2.5.2, R226) | 2.6.1 | upgraded | 2.6.0 removes the 100-genome cap in dense genera (#671, wrong species calls) and the pplacer alpha20 crash; 2.6.1 fixes check_install; still R226 (2.7+ needs R232) |
| anvi'o | 9 (LiSC anvio-9) | 9 | kept | no official v9 image (Docker Hub tags stop at 8) and the Bioconda anvio-minimal image lacks HMMER/DIAMOND/Prodigal/SAMtools: custom image (anvio-minimal 9 + tools + nodejs 22, SCG DMNDs built in) |
| barrnap | 0.9 (LiSC module) | 0.9 | kept | 1.10.x is a rewrite with changed feature names |
| tRNAscan-SE | 2.0.12 (Bakta/anvio envs) | 2.0.12 | kept | 2.0.13 only fixes legacy-mode compilation with GCC 15 |
| genome_uploader | 2.5.1 (env ena_tools) | 2.5.1 | kept | the ENA scripts were tested with 2.5.1; 3.x registers samples through the webin-v2 queue and needs re-testing |
| Webin-CLI | 9.0.1 (env ena_tools) | 9.0.1 | kept | shipped with genome-uploader 2.5.1; ENA enforces a minimum version server-side, so switch to a newer image if Webin-CLI reports that its version is no longer supported |
| assembly_uploader | not used before | 1.3.6 | new | latest; provides assembly_manifest/study_xmls/submit_study |
| mgnify-pipelines-toolkit | 1.4.8 (env ena_tools) | 1.4.8 | kept | provides webin_cli_handler for genome-uploader 2.5.1 |

## Pin files

`<env>.linux-64.pin.txt` lists exact package URLs with md5 sums (format of `conda list --explicit --md5`).
Snakemake (`--sdm conda`) creates the env from the pin file when one exists next to the YAML, and the custom
Dockerfiles/.def files install the same file, so the conda and container routes have identical packages.
They were generated with `micromamba create --dry-run --json -f <env>.yaml --platform linux-64 --strict-channel-priority`
(GPU envs with `CONDA_OVERRIDE_CUDA=12.9`, GenomeFace with `CONDA_OVERRIDE_CUDA=11.7` and flexible priority).
Regenerate after editing a YAML (e.g. `snakedeploy pin-conda-envs workflow/envs/<env>.yaml`) and commit both files.
`taxconverter.yaml` has no pin file on purpose: it needs a pip step that pin files cannot express.

Special cases:
- `genomeface` needs **flexible** channel priority (the nvidia channel shadows conda-forge `cuda-version`/`cudnn`/`nccl`
  under strict priority). Its pin file reproduces the user's LiSC env package-for-package except
  pyrodigal 2.3.0 / pyhmmer 0.10.15 (instead of 3.6.3.post1 / 0.11.1): `markersgf` calls `pyrodigal.OrfFinder`
  (renamed `GeneFinder` in pyrodigal 3.0) and the bytes-name pyhmmer API. The image that reproduced the user's env
  exactly failed with `AttributeError: module 'pyrodigal' has no attribute 'OrfFinder'`.
- `vamb-gpu` carries a second env `/opt/conda/envs/taxconverter` (taxconverter 1.1.0 pins python<3.12, pandas 2.0.0,
  numpy 1.26.4); `taxconverter` is linked into `/opt/conda/bin`.
- `checkm2`/`binette` pin CheckM2 build `pyh7e72e81_1` (TensorFlow 2.17); `gtdbtk` pins build `pyh1f0d9b5_2` and skani 0.3.1.

## Using the images with Snakemake on LiSC

- Run with `--sdm apptainer`; GPU rules need `--nv`, e.g. `--apptainer-args "--nv -B /lisc/opt/mirror,/lisc/opt/sw/software/CheckM2,/lisc/data/scratch"`.
  On nodes without an NVIDIA driver Apptainer only warns that no nv libraries were found (documented behaviour;
  not tested on LiSC CPU nodes here), so `--nv` can be passed globally or set per rule.
- Set `APPTAINER_CACHEDIR` (and `--apptainer-prefix`) to a scratch directory, not `$HOME`.
- Pre-pulled SIFs of every image in this manifest are in `<sif_dir>/`
  (custom builds in `sif/custom/`); a LiSC config can point `config["containers"]` at these files until the custom
  images are published.
- GPU jobs request `gpu:N` (Slurm `GresTypes = gpu,localtmp`; the `shard` GRES used by the old profile no longer exists).
- GTDB-Tk: `GTDBTK_DATA_PATH=/lisc/opt/mirror/GTDBTk/r226` (bind `/lisc/opt/mirror`). Binette/CheckM2: pass the
  DIAMOND DB path explicitly. dRep 3.6.2: pass `--genomeInfo` from CheckM2 (genome,completeness,contamination CSV)
  or it tries to run CheckM1.

## Building and publishing the custom images

1. Commit the tree (recipes, env YAMLs, pin files, workflow) and push to `main` of github.com/sprdmt/bass-fastq-to-mag.
   `.github/workflows/containers.yml` builds only the images whose `containers/<image>/**`, env YAML or pin file
   changed (all images on changes to the workflow file or via *Run workflow*), and pushes
   `ghcr.io/sprdmt/bass-fastq-to-mag-<image>:<TOOL_VERSION>` and `:<TOOL_VERSION>-<short sha>`. Pull requests build
   without pushing. The runner first removes the preinstalled .NET, Android, GHC, CodeQL, Boost and Swift trees and prunes Docker images,
   because the GPU images are 6-7 GB even as compressed SIFs (sizes in manifest notes / LISC_NOTES.md).
2. After the first push, open each package on GitHub (Profile -> Packages -> bass-fastq-to-mag-<image> ->
   Package settings), link it to the repository if needed and set visibility to **Public** (GHCR packages start private).
3. Record the pushed digests (`docker buildx imagetools inspect ghcr.io/sprdmt/bass-fastq-to-mag-<image>:<version>`)
   in `manifest.tsv` and, for full immutability, reference `@sha256:<digest>` in the config.
4. Verify on LiSC: `apptainer pull <image>.sif docker://ghcr.io/sprdmt/bass-fastq-to-mag-<image>:<version>` and rerun
   the smoke tests listed in the manifest (`verify_cmd`; GPU images with `apptainer exec --nv`).
Local builds without CI: `docker build -f containers/<image>/Dockerfile -t <tag> .` from the repo root, or on LiSC
`apptainer build --fakeroot <image>.sif containers/<image>/<image>.def` (works on compute nodes, see LISC_NOTES.md).

## Verification on LiSC (2026-10-31)

- All 28 envs (29 at the time, incl. TaxonKit, since removed as unused): image pulled (or custom image built with `apptainer build --fakeroot`) and the version command matched the
  pinned version (`status` and `verified_version_string` in `manifest.tsv`); CLI help of every tool captured from inside the images.
- GPU node (Tesla T4): `genomeface -h` lists the GPU for CuPy and TensorFlow; TF and CuPy kernels ran on the GPU;
  `torch.cuda.is_available()` and a GPU matmul passed in semibin2-gpu and vamb-gpu (PyTorch 2.8.0, CUDA 12.9);
  the Bioconda MMseqs2 18 image ran `mmseqs search --gpu 1` (20/20 correct top hits).
- Functional checks: `markersgf` on a 213-contig MAG (22 marker hits, 20 COGs); anvi'o contigs DB + HMMs + SCG taxonomy on the
  same MAG and `anvi-self-test --suite mini`; `gtdbtk check_install` against `/lisc/opt/mirror/GTDBTk/r226`.
- Not covered: full runs of the binners on real data (left to the end-to-end test run), the GHCR images themselves (not yet
  pushed), GPUs other than the T4 (Slurm routing, see LISC_NOTES.md).
