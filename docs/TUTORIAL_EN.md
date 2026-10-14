# nanoCASTLE — image user guide

> nanoCASTLE (**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer) was developed under the code name *Nanomatizer*: folders, scripts and output files inside the image still use that name (`/opt/nanomatizer`, `Nanomatizer_<SAMPLE>_aggregated.xlsx`, ...).

This guide is for people who **download the Docker image** and want to run
the pipeline on their own sample. Nothing to install: callers, pipeline,
annotation databases and R environments are already in the image —
except three external databases you download once (§0).

**The pipeline is GRCh38-only.** The genome must be GRCh38 (hg38): target,
GFF, ENCODE blacklist, EnsDb and the chromosome names used by the callers
are all hg38. It will not work on hg37/GRCh37.

**Chromosome prefix**: chromosome names must be **consistent between the
BAM and the reference** — the format used in production is the
UCSC/Ensembl style with the `chr` prefix (`chr1`, `chr2`, … `chrX`,
`chrY`) both in the BAM `@SQ` headers and in the reference `.fai` index.
If your BAM or reference uses names without `chr` (NCBI style), the
pipeline does not support that: regenerate the files with `chr…` names.

What you need on your machine (file names on your disk are free; what
matters is the path **inside** the container, §1):

- your **WGS BAM** (primary chromosomes, GRCh38) — for runs with callers
  (cases 1-4);
- the **GRCh38 genome** `.fa` (e.g. the no-alt
  `GCA_000001405.15_GRCh38_no_alt_analysis_set` from UCSC/Ensembl) — for
  runs with callers (cases 1-4). The `.fai` index is **not required**: the
  pipeline checks for it and builds it automatically if missing (§1);
- a **results folder** to mount into the container (§1) — always (it does
  **not** have to be empty: the pipeline is idempotent, see case 1);
- the two **DepMap** CSVs (`CRISPRGeneEffect.csv` + `Model.csv`) — only
  for runs with STEP6 (cases 1, 2, 5,6) **if you pass `-T`** (see §4);
  without `-T` the DepMap annotation is skipped and they are not needed.

---

## 0. External databases (one-time download)

Some resources are **not shipped with the image** because their terms of
use do not allow redistribution (or because they are large and optional).
You download them once from the original providers, into a folder of your
choice (`~/nanocastle_db` in the examples), and mount that folder at
`/opt/nanomatizer/userdb` in every run that annotates:

| File | Resource | Needed for |
|---|---|---|
| `OncoKB.tsv` | OncoKB Cancer Gene List | annotation (**required**) |
| `population_cnv_grch38.txt` | DECIPHER population CNVs | annotation (**required**) |
| `trrust_rawdata.human.tsv` | TRRUST v2 | annotation (**required**) |
| `DGV.GS.hg38.gff3` | DGV Gold Standard (hg38, ~380 MB) | annotation (**required**) |
| `CRISPRGeneEffect.csv`, `Model.csv` | DepMap Public | only with `-T` (optional) |

```bash
mkdir -p ~/nanocastle_db
docker run --rm --user $(id -u):$(id -g) \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 bash Nanomatizer/download_databases.sh
```

Add `--no-depmap` to skip the two DepMap CSVs (~430 MB) if you will not
use `-T`. For DepMap the script looks for the most recent "DepMap Public"
release on the portal; if the portal does not answer to a script (anti-bot
verification) it downloads 24Q4, the last release with a stable public
URL. The downloaded release is written in `DEPMAP_RELEASE.txt`. To use
another release just copy your own `CRISPRGeneEffect.csv` + `Model.csv`
(from https://depmap.org/portal/data_page/?tab=allData) into the folder.
By downloading the files you accept the terms of use of
each provider (see `/opt/nanomatizer/THIRD_PARTY.md`). Caller-only runs
(`-m callers`) do not need this folder.

---

## 1. The mounts

The paths **inside** the container are fixed (the scripts use them as-is);
the path **on your host** can be anything — just link them with
`-v host:container`.

For runs with callers (cases 1-4) you need:

- your BAM at the **fixed** internal path
  `/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam` — the name is always
  this one, independent of the sample (§2);
- your GRCh38 genome at the internal path
  `/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta`;
- the results folder at the internal path
  `/opt/nanomatizer/Nanomatizer/Results`.

**Indexes are not an option.** The pipeline checks that
`BAM.srt.bam.bai` and `<genome>.fai` exist; if they are missing, it builds
them automatically (`samtools index` / `samtools faidx`) and continues.
If your BAM already has a `.bai`, you can mount it next to it (same name +
`.bai`) to skip the on-the-fly index; if you do not have one, the pipeline
generates it. Same for the genome `.fai`.

For runs with STEP6 (cases 1, 2, 5, 6) add the database folder of §0 at
the internal path `/opt/nanomatizer/userdb` (it holds `OncoKB.tsv`,
`population_cnv_grch38.txt`, `trrust_rawdata.human.tsv`,
`DGV.GS.hg38.gff3` and, for `-T`, `CRISPRGeneEffect.csv` +
`Model.csv`).

Custom depmap files can be mounted **anywhere** inside the
container; pass the internal path to `-d` (e.g. mount under
`/opt/nanomatizer/depmap/` and pass
`-d /opt/nanomatizer/depmap/myfile1.csv,...`).

Notes:

- if your files have different names, mount them under the internal name
  anyway, e.g.
  `-v ~/data/myrun.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam`.
- All results (outputs, logs, `STATUS.txt`) land in the **results folder**
  you mount: when the container exits, the files are already on your
  machine.
- The results folder **does not have to be empty**: if it already holds
  valid results from a previous run, completed steps are skipped
  (idempotency, case 1).
- For cases 5-6 the results folder must already contain the caller results
  produced by a previous run (cases 1-4).

---

## 2. The sample name (`SAMPLE`) — mandatory

`-e SAMPLE=<name>` tells the pipeline the sample name. The name is used in
**all** result paths (per-tool folders, per-tool CSVs, final xlsx name) and
must be the same across every case you run on the same experiment.

Example: `-e SAMPLE=MySample`. **There is no default**: the pipeline is
generic and stops if `SAMPLE` is not provided.

The BAM name, in contrast, does **not** depend on the sample: inside the
container it is always `BAM.srt.bam` (§1).

---

## 3. The six use cases

The pipeline is launched like this (inside the image):

```
bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh \
  [-m MODE] [-t TOOLS] [-T TUMOR] [-d DEPMAP] [-p NPROC]
```

- `-m` = `full` (default) | `callers` | `annotate`
- `-t` = `all` (default) | list from: `sniffles` | `cutesv` | `gasoline` |
  `intertrans` | `gladiator`
- `-T` = tumor type for the DepMap annotation (§4) — **optional**; without
  `-T` (and without custom `-d` files) the DepMap annotation is skipped
- `-d` = `none` (default) | `FILE1,FILE2,...` (custom files, §4)
- `-p` = number of cores to use (default: all available cores)

In `full`/`callers`, `-t` selects which **callers** to run
(`gladiator` = NanoGLADIATOR; `intertrans` implies Gasoline, which it
derives from). In `annotate`, `-t` selects which **tools** to annotate and
aggregate.

### Case 1 — Full run

All 4 callers (NanoGLADIATOR → Sniffles → cuteSV → Gasoline) +
annotation + final aggregation. This is the complete production run.
Without `-T` the aggregate comes out **without** DepMap columns (§4).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/data/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh
```

To add the DepMap annotation for your tumor, add `-T` (the two DepMap CSVs
must be in the database folder, §0), e.g. for a lung tumor:

```bash
  ...
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -T Lung
```

- Time: **many hours** (CPU-dependent; cuteSV is the slowest one). For a
  new sample, the first step (NanoGLADIATOR) includes the "prepare" phase
  of ~35-40 minutes at 10 processors.
- RAM: keep plenty of free RAM (cuteSV is the most demanding tool; the
  pipeline requires ≥ 4 GB free before starting it).
- When it finishes, everything is in `~/nanocastle_results` (§5).
- The run is **idempotent**: if you relaunch it (or resume from a
  partially filled results folder), already-completed and valid steps are
  skipped and it resumes where it stopped.
- If a step fails the pipeline **stops** (it cleans nothing): read
  `STATUS.txt` and `logs/pipeline.log` in the results folder.
- Run the container **as your own user** by adding
  `--user $(id -u):$(id -g)` to `docker run`: the results are then owned
  by you and the results folder only has to be writable by you. Without
  `--user` the pipeline runs as uid/gid 1000 and the results folder must
  be writable by that uid.

### Case 2 — Full run on selected tools only

Like case 1, but only the callers chosen with `-t` are run and STEP6
annotates/aggregates just those. Example: only Sniffles 2 + NanoGLADIATOR:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (same mounts as case 1) ... \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full -t sniffles,gladiator
```

Other examples: `-m full -t cutesv`, `-m full -t gasoline,intertrans`.
Unselected caller steps are skipped; STEP6 produces the aggregate with
only the selected tools (the other slots stay empty).

### Case 3 — Caller-only run (all callers)

Runs the 4 callers and stops **before** STEP6: no annotation, no
aggregation. DepMap files are not needed.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/data/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers
```

### Case 4 — Caller-only run on selected callers

Like case 3, but only the callers chosen with `-t`. Example: only cuteSV
+ Gasoline:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (same mounts as case 3) ... \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers -t cutesv,gasoline
```

If a caller's result already exists and is valid, that caller is skipped
(idempotency).

### Case 5 — Annotation-only run (all tools)

No callers: takes the results **already produced** (cases 1-4) in the
results folder and runs the full STEP6 (annotation + aggregation) on all
tools. No BAM or genome needed.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate
```

Without `-T` the aggregate comes out without DepMap columns. To annotate
your tumor, add `-T <tumor>` (DepMap CSVs in the database folder, §0).

### Case 6 — Annotation-only run on selected tools

Like case 5, but annotates/aggregates **only the tools chosen** with
`-t`. Example: only Sniffles + NanoGLADIATOR:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate -t sniffles,gladiator
```

The per-tool CSVs of excluded tools, if present, are zeroed (0 rows)
before aggregation; if absent, the module substitutes an empty
placeholder. The final aggregate contains only the selected tools.

---

## 4. DepMap annotation (`-T` and `-d`)

STEP6 can add DepMap columns to the aggregate. **By default it does not**:
without `-T` (and without custom `-d` files) the DepMap annotation is
**skipped** and the aggregate comes out without DepMap columns. To enable
it, choose one of the two paths (they are **mutually exclusive**):

### 4a. `-T TUMOR` — the tumor type of your sample

Pass the tumor type:

- the exact `OncotreeLineage` value present in `Model.csv` of the DepMap
  release you are using (e.g. `Lung`, `CNS/Brain`,
  `Ovary/Fallopian Tube`, `Ovary`); or
- `<lineage>+Pediatric` for the **pediatric-only** subset of that lineage
  (e.g. `-T "CNS/Brain+Pediatric"` selects models with
  `OncotreeLineage == CNS/Brain` and `AgeCategory == Pediatric`); or
- `MBL` for hematological models (`DepmapModelType == MBL`).

The pipeline reads `Model.csv` + `CRISPRGeneEffect.csv` from the database
folder mounted at `/opt/nanomatizer/userdb` (§0), selects the models of that tumor and
generates **one** DepMap group with the columns:

```
DepMap.<TUMOR>.Min  DepMap.<TUMOR>.Q1  DepMap.<TUMOR>.Median
DepMap.<TUMOR>.Q3   DepMap.<TUMOR>.Max  <TUMOR>.Sum.Score  <TUMOR>.Sum.Score.match
```

where `<TUMOR>` is the sanitized name (e.g. `Ovary/Fallopian Tube` →
`Ovary_Fallopian_Tube`). If the tumor has no models in `Model.csv`, the
pipeline **stops** with a message listing the available lineages.

Example: full run with DepMap for a lung tumor:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/data/MySample_no_alt.srt.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/data/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full -T Lung
```

### 4b. `-d FILE1,FILE2,...` — custom files

DepMap annotation with **your own** files, in any number (also 1, also 5,
as many as you want). Each file has the same format as
`CRISPRGeneEffect.csv` (column 1 = ModelID, following columns = genes with
scores). The path you pass is the one **inside the container** (mount the
files wherever you like, §1). For each valid file the 5 columns
`DepMap.<name>.Min/Q1/Median/Q3/Max` are generated; invalid files are
skipped with a warning (if all are skipped the DepMap annotation is empty,
the run continues).

Example: full run with 2 custom files:

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  ... (BAM/GRCh38/results mounts of case 1) ... \
  -v ~/data/my_depmap_1.csv:/opt/nanomatizer/depmap/my_depmap_1.csv \
  -v ~/data/my_depmap_2.csv:/opt/nanomatizer/depmap/my_depmap_2.csv \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m full \
    -d /opt/nanomatizer/depmap/my_depmap_1.csv,/opt/nanomatizer/depmap/my_depmap_2.csv
```

`-T` and `-d FILE` **cannot** be combined: the pipeline stops with an
error.

---

## 5. Output (where to find it)

Everything lands in the results folder you mounted
(`~/nanocastle_results` in the examples; in the container
`/opt/nanomatizer/Nanomatizer/Results`):

```
Results/
├── NanoGLADIATOR/<SAMPLE>/R6_W1Mb/      # NanoGLADIATOR results (1 Mb windows)
├── Sniffles2/<SAMPLE>/                  # Sniffles results
├── CuteSV/<SAMPLE>/                     # cuteSV results
├── Gasoline/Parse/<SAMPLE>/             # Gasoline results
├── <per-tool annotation CSV/xlsx>       # cases 1, 2, 5, 6
├── Nanomatizer_<SAMPLE>_aggregated.xlsx # final aggregate (cases 1, 2, 5, 6)
├── STATUS.txt                           # state: last stop/completion
└── logs/                                # pipeline.log, step6_v2.log
```

- If the pipeline **succeeded**: `STATUS.txt` has the completion line and
  the aggregate xlsx is present (in the cases that produce it).
- If it **stopped**: `STATUS.txt` says why, `logs/pipeline.log` has the
  detail. Cleaning "dirty" folders is up to you: the script never removes
  anything.

---

## 6. Messages

All progress and error messages are written in English to the terminal
and to `Results/logs/`. `STATUS.txt` records the state of the last
stop/completion.

---

## 7. Notes

- **GRCh38-only** and **`chr` prefix** consistent between BAM and
  reference: see the top of this guide.
- **1 Mb coverage windows**: the `coverage` column of the aggregate and
  the `NanoGladiator.W1Mb.<SAMPLE>...exons.csv` file use 1-megabase
  windows. If the mounted genome is the standard GRCh38 (the pipeline
  recognizes it automatically from the `.fai` index), the target
  `hg38w1000000` shipped in the image is used as-is; if you mount a
  genome variant with the same chromosome layout (same 24 `chr`
  chromosomes, same lengths), the pipeline builds the matching 1 Mb
  windows automatically on the first run (a few minutes, one time) and
  reuses them on later runs. If the genome does not match this layout,
  the pipeline stops with an explicit message.
- **Automatic indexes**: `.bai` and `.fai` are checked and, if missing,
  built by the pipeline (§1): they are not an option.
- **Optional DepMap**: without `-T` (and without custom `-d` files) no
  DepMap columns are generated and the two DepMap CSVs are not needed
  (§4).
- The default internal user is `nanomatizer` (uid 1000); the image also
  runs with any `--user UID:GID` (recommended: your own, see case 1).
- Each caller has its own structural check (files present and non-empty,
  correct headers): if a check fails the pipeline STOPS and notifies,
  without cleaning anything.
- The run is idempotent: relaunching it never regenerates valid steps.
- `-p N` limits the run to N cores (useful on shared machines); the
  default is to use all available cores.
