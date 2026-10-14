# nanoCASTLE — how to run it

> nanoCASTLE (**nano**pore **C**opy-number **A**nd **ST**ructural variant **L**andscape **E**xplorer) was developed under the code name *Nanomatizer*: folders, scripts and output files inside the image still use that name (`/opt/nanomatizer`, `Nanomatizer_<SAMPLE>_aggregated.xlsx`, ...).

Short guide: the `docker run` commands and the options. For details
(depmap formats, output layout, notes) see `TUTORIAL_EN.md`.

**Before you start**

- Genome: **GRCh38 only**.
- BAM and genome must use the same chromosome names (`chr1` in both,
  not `chr1` in one and `1` in the other).
- Files you need: your sorted BAM and the genome FASTA (the `.bai` and
  `.fai` indexes are **not required**: the pipeline checks them and
  builds them automatically, see `TUTORIAL_EN.md` §1) and — **only if
  you use `-T`** for the DepMap annotation — the DepMap files
  `CRISPRGeneEffect.csv` + `Model.csv` (not needed without `-T`, see
  section 5).

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

## The one rule

Every `-v` line is `YOUR_FILE:CONTAINER_PATH`.

- **Left side**: your file. Any name, any folder.
- **Right side**: copy it **exactly** as written.

Replace `MySample` with your sample name everywhere it appears
(`SAMPLE=`). The BAM name **inside** the container is always
`BAM.srt.bam` (fixed, independent of the sample). Use the same
`SAMPLE` name for every run of that sample.

---

## 1. Run everything

All callers (NanoGLADIATOR, Sniffles, cuteSV, Gasoline) + annotation.
Without `-T` there are no DepMap columns (section 5).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/data/MySample.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/data/genome.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh
```

To add the DepMap annotation for your tumor: add `-T <tumor>` (e.g.
`-T Lung`, `-T "CNS/Brain+Pediatric"` for the pediatric-only subset of
that lineage, or `-T MBL` for hematological models) (the two
DepMap CSVs must be in the database folder, section 0).

> `Databsese` is not a typo in this guide: the folder inside the
> container is really called that way.

Results go to `~/nanocastle_results`. If the run is **interrupted**,
launch the same command again: callers already completed are
skipped. If it **keeps failing** at the same step, read `STATUS.txt`
and `logs/pipeline.log` in the results folder.

---

## 2. Callers only (no annotation)

No DepMap files. Choose the callers with `-t` (without `-t` = all
four).

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/data/MySample.bam:/opt/nanomatizer/Nanomatizer/Test/BAM/BAM.srt.bam \
  -v ~/data/genome.fasta:/opt/nanomatizer/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m callers -t sniffles,cutesv
```

---

## 3. Annotation only

Works on results already in the folder. No BAM, no genome. Use the
**same results folder** as the caller run.

```bash
docker run -it --rm --user $(id -u):$(id -g) \
  -e SAMPLE=MySample \
  -v ~/nanocastle_results:/opt/nanomatizer/Nanomatizer/Results \
  -v ~/nanocastle_db:/opt/nanomatizer/userdb \
  gmatteigmattei/nanocastle:1.0.0 \
  bash Nanomatizer/Scripts/NewScripts/run_pipeline_sequential.sh -m annotate -t sniffles,cutesv
```

With `-T <tumor>` (DepMap CSVs in the database folder, section 0) the
aggregate includes the DepMap columns for that tumor.

`-t` says which tool results to annotate:

- every tool you list **must** have its results in the folder,
  otherwise the pipeline stops;
- tools you don't list stay empty in the final xlsx, even if their
  results are in the folder.

Without `-t`, all five outputs must be present: `sniffles`, `cutesv`,
`gladiator`, `gasoline`, `intertrans`.

---

## 4. Adding callers one at a time

Same results folder every time:

```bash
... -m callers -t sniffles           # command from section 2
... -m callers -t cutesv             # command from section 2
... -m annotate -t sniffles,cutesv   # command from section 3
```

Callers already completed are skipped, so nothing is recomputed.

---

## 5. DepMap annotation (optional)

**By default there is none**: without `-T` (and without custom `-d`
files) the DepMap annotation is skipped and the aggregate has no
DepMap columns. The two DepMap CSVs are not needed.

**`-T <tumor>`**: DepMap annotation for the sample's tumor. Pass the
exact `OncotreeLineage` value from `Model.csv` (e.g. `Lung`,
`CNS/Brain`, `Ovary/Fallopian Tube`), or `<lineage>+Pediatric` for the
pediatric-only models of that lineage (e.g. `CNS/Brain+Pediatric`), or
`MBL` for hematological models. The two DepMap CSVs must be in the database folder (section 0):

```bash
  ... run_pipeline_sequential.sh -m annotate -T Lung
```

**Your own files** (`-d`, same format as `CRISPRGeneEffect.csv`):
mount them anywhere and pass the container paths to `-d`. The two
standard DepMap CSVs are then not needed. `-T` and `-d FILE` cannot be
combined.

```bash
  -v ~/data/my_screen.csv:/opt/nanomatizer/depmap/my_screen.csv \
  ...
  ... run_pipeline_sequential.sh -d /opt/nanomatizer/depmap/my_screen.csv
```

---

## Options

| Option | Values | Default |
|---|---|---|
| `-m` | `full` (callers + annotation), `callers`, `annotate` | `full` |
| `-t` | `sniffles`, `cutesv`, `gladiator`, `gasoline`, `intertrans`, `all` (comma-separated) | `all` |
| `-T` | tumor type (`OncotreeLineage` from Model.csv, `<lineage>+Pediatric`, or `MBL`) — §5 | (none: DepMap skipped) |
| `-d` | `none`, `FILE1,FILE2,...` (custom files) — §5 | `none` |
| `-p` | number of cores to use | all |

- `intertrans` is a second output of the Gasoline run: `-m callers
  -t intertrans` runs Gasoline.
- `-T`/`-d` are ignored with `-m callers`.

## Output

In your results folder:

- one folder per tool;
- annotation CSV/xlsx;
- `Nanomatizer_<SAMPLE>_aggregated.xlsx` (final table);
- `STATUS.txt`;
- `logs/`.
