# nanoCASTLE

**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer

Structural variant (SV) and copy-number variant (CNV) calling and annotation
for Oxford Nanopore whole-genome sequencing data aligned to GRCh38.

nanoCASTLE runs four callers on one sorted BAM and merges their calls into a
single annotated table:

| Step | Tool | Output |
|---|---|---|
| 1 | NanoGLADIATOR 1.0 (1 Mb windows) | CNV segments |
| 2 | Sniffles2 2.8.0 | SV VCF |
| 3 | cuteSV 2.1.4 | SV VCF |
| 4 | GASOLINE 2.0 (Parse + GermlineDetector) | SV tables, inter-chromosomal translocations |
| 5 | Annotation + aggregation | per-tool CSVs and `Nanomatizer_<SAMPLE>_aggregated.xlsx` |

The pipeline is distributed as a Docker image with all tools, conda
environments and redistributable annotation databases already installed.

## Note on names

nanoCASTLE (**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer) was developed under the code name *Nanomatizer*: folders, scripts and output files inside the image still use that name (`/opt/nanomatizer`, `Nanomatizer_<SAMPLE>_aggregated.xlsx`, ...).

## Requirements

- Docker (any user that can run `docker run`; the container runs with your
  UID via `--user $(id -u):$(id -g)`).
- About 8 GB of disk for the image and 1.3 GB for the external databases.
- A coordinate-sorted, indexed BAM of Oxford Nanopore reads aligned to
  **GRCh38 no-alt analysis set** (chromosome names with `chr`), e.g. with
  minimap2.
- The same reference FASTA and its `.fai` index. It can be downloaded from NCBI:
  ```bash
  wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
  gunzip GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
  mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
  samtools faidx GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta
  ```
  (if the `.bai` or `.fai` is missing the pipeline builds it, but it is lost
  at the end of the run unless you mount it.)
- A full run takes several hours: in our test (a 240 GB BAM, `-p 20`) about
  7.5 h, with GASOLINE as the slowest step, and up to ~60 GB of RAM.

## Quick start

```bash
docker pull gmatteigmattei/nanocastle:1.0.0

# 1. one-time download of the databases that are not shipped with the image
mkdir -p ~/nanocastle_db
docker run --rm --user $(id -u):$(id -g) \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 bash Nanomatizer/download_databases.sh

# 2. full run (callers + annotation) on one BAM
mkdir -p ~/nanocastle_results
docker run --rm --user $(id -u):$(id -g) -e SAMPLE=MySample \
  -v /path/to/sample.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam:ro \
  -v /path/to/sample.srt.bam.bai:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam.bai:ro \
  -v /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:ro \
  -v /path/to/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.fai:ro \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full
```

Modes: `-m full` (callers + annotation), `-m callers`, `-m annotate` (on
results already present). `-t` selects a subset of callers, `-T` adds the
DepMap annotation for a tumour lineage, `-p` sets the number of cores.

The complete guide, with every use case, is in [`docs/`](docs):
[TUTORIAL_EN.md](docs/TUTORIAL_EN.md) (full),
[TUTORIAL_essenziale_EN.md](docs/TUTORIAL_essenziale_EN.md) (short), and the
Italian versions [TUTORIAL.md](docs/TUTORIAL.md) /
[TUTORIAL_essenziale.md](docs/TUTORIAL_essenziale.md).

## External databases

OncoKB, DECIPHER, TRRUST, DGV and DepMap are not included in the image: they
are downloaded from the original providers by
[`download_databases.sh`](download_databases.sh) (it also runs on the host;
it needs `curl`, `gzip` and `md5sum`). By downloading them you accept the
terms of use of each provider. Sources and terms of every bundled or
downloaded resource are listed in [THIRD_PARTY.md](THIRD_PARTY.md).

## Repository content

| Path | Content |
|---|---|
| `Scripts/NewScripts/` | pipeline driver (`run_pipeline_sequential.sh`) and annotation module (`step6/`) — the same files installed in the image |
| `download_databases.sh` | downloader for the external databases |
| `docs/` | tutorials |
| `example/` | small example DepMap tables |
| `docker/envspecs/` | exact package list (conda explicit specs) of the four environments in the image |
| `THIRD_PARTY.md` | licensing of every component, third-party software and data with their terms |
| `LICENSE` | MIT licence of the nanoCASTLE code |

## Licence

The code in this repository is released under the [MIT licence](LICENSE)
(Copyright (c) 2026 Azienda Ospedaliero-Universitaria Meyer IRCCS).

The Docker image also contains NanoGLADIATOR and GASOLINE, which are licensed
for **noncommercial use only** (PolyForm Noncommercial 1.0.0; the licence
text is in `/opt/nanomatizer/NanoGLADIATOR_1.0/LICENSE` and
`/opt/nanomatizer/GASOLINE_2.0/LICENSE` inside the image). Commercial users
can run nanoCASTLE without them by selecting the other callers:

```bash
... run_pipeline_sequential.sh -m full -t sniffles,cutesv
```

All other bundled tools and databases keep their own licences, listed in
[THIRD_PARTY.md](THIRD_PARTY.md).

## Reproducibility note

Sniffles2, cuteSV and NanoGLADIATOR are deterministic on the same BAM.
GASOLINE extracts the SV signatures in parallel and their order is not fixed
between runs; the detector is order-dependent, so two complete runs on the
same BAM can differ in a small fraction of the GASOLINE calls (about 2% in
our test). Re-running the detector on the same extracted signatures gives
identical results.
