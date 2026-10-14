#!/bin/bash
# =============================================================================
# SEQUENTIAL SV-CALLER PIPELINE (hg38/GRCh38 — GRCh38/hg38 ONLY) — one script,
# one tool at a time.
#
# This is a GENERIC pipeline: the sample name is provided by the user.
#   SAMPLE   sample name — REQUIRED (env or -e SAMPLE=...), no default
#   BAM      BAM path (default $T/Nanomatizer/Test/BAM/BAM.srt.bam). The .bai
#            index is NOT an option: the pipeline checks $BAM.bai and, if
#            missing, builds it with samtools (idempotent).
#   REF      reference FASTA GRCh38 path (default
#            $T/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta).
#            Same mechanism for the index: $REF.fai checked/built.
#
# Modes (-m):
#   full      (default) — all callers (or the -t subset) + full STEP6
#   callers           — CALLERS ONLY (all or the -t subset), no annotation/
#                       aggregation
#   annotate          — STEP6 ONLY (annotation + aggregation) on the results
#                       already in OUT (all tools or the -t subset); does not
#                       need BAM/REF
#   -t TOOLS   : all | list from sniffles|cutesv|gasoline|intertrans|gladiator
#   -T TUMOR   : tumor type for the DepMap annotation (STEP6).
#                Value = OncotreeLineage of Model.csv (e.g. Lung,
#                Ovary/Fallopian Tube) or MBL (hematologic models,
#                DepmapModelType==MBL). The +Pediatric suffix selects the
#                pediatric-only subset of that lineage (e.g.
#                "CNS/Brain+Pediatric"). The two DepMap CSVs
#                (CRISPRGeneEffect.csv + Model.csv) must be mounted in
#                /opt/nanomatizer/Nanomatizer/Databsese/.
#                WITHOUT -T (and without a custom -d file) DepMap is SKIPPED:
#                there is no longer a default CNS.
#   -d DEPMAP  : none | FILE1,FILE2,... (0..N custom depmap files in
#                CRISPRGeneEffect format, container paths). Mutually exclusive
#                with -T. Used only by full/annotate (STEP6).
#   -p NPROC   : number of cores to use (default: all available, nproc). The
#                per-caller threads stay as in the protocol, capped to this
#                number.
#
# Order (fixed by the user):
#   1. NanoGLADIATOR (run6, R server, 1 Mb windows — target hg38w1000000)
#   2. Sniffles2
#   3. cuteSV
#   4. Gasoline 2.0 (Parse + GermlineDetector + copy Parse->Germline)
#   5. Annotations + aggregation (module v2: module_v2/run_module.sh)
#
# BRIEF rules (source of truth):
#   - Never modify tool sources, never regenerate validated targets,
#     never rm/mv data. One tool at a time. STOP at the first failure.
#   - Idempotency: each step is SKIPPED if its marker/check says "valid".
#   - "Exists but not valid" state (e.g. a dirty Gasoline dir): STOP — cleaning
#     the dir is the USER's job (BRIEF 2a); the script removes nothing.
#   - Only one pipeline instance at a time (flock).
#   - Per-caller threads capped to the core count (-p, default nproc);
#     free RAM >= 4 GB before cuteSV.
#
# Usage:  bash run_pipeline_sequential.sh [-m full|callers|annotate] [-t TOOLS]
#        [-T TUMOR] [-d DEPMAP] [-p NPROC]
# Logs:   $OUT/logs/pipeline.log  +  $OUT/STATUS.txt
# =============================================================================
set -uo pipefail

# Generic layout: T = 3 levels above this folder (…/NewScripts).
# In the image: /opt/nanomatizer. No user path hardcoded.
T="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
# micromamba: host → $T/micromamba; image → /opt/micromamba/bin/micromamba.
# Overrides: export MM=/path/micromamba before launching.
MM=${MM:-$T/micromamba}
[ -x "$MM" ] || MM=/opt/micromamba/bin/micromamba
export MAMBA_ROOT_PREFIX=$T/mamba-root

SAMPLE=${SAMPLE:-}
# FIXED internal BAM name (that is how the in-container mount is named). The
# .bai index is NOT an option: it is checked/built automatically below.
BAM=${BAM:-$T/Nanomatizer/Test/BAM/BAM.srt.bam}
REF=${REF:-$T/reference/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta}
# Results directory: $T/Nanomatizer/Results (container: /opt/nanomatizer/Nanomatizer/Results)
OUT=$T/Nanomatizer/Results
LOGDIR=$OUT/logs
STATUS=$OUT/STATUS.txt
PLLOG=$LOGDIR/pipeline.log
LOCK=$OUT/logs/pipeline.lock

mkdir -p "$LOGDIR"
step()  { echo "[$(date '+%F %T')] [PIPELINE] $*" | tee -a "$STATUS" >> "$PLLOG"; }
die() { step "STOP: $1"; exit 1; }

# ---------- CLI options -------------------------------------------------------
#   -m MODE   : full (default) | callers | annotate
#                full     -> all callers (or the -t subset) + full STEP6
#                callers  -> CALLERS ONLY (all or the -t subset), no STEP6
#                annotate -> STEP6 ONLY (annotation+aggregation) on -t, no callers
#   -t TOOLS  : all | list (sniffles|cutesv|gasoline|intertrans|gladiator)
#               (default all) — in full/callers selects the callers, in annotate
#               the tools to annotate/aggregate
#   -T TUMOR  : DepMap tumor type (OncotreeLineage of Model.csv, or MBL).
#               +Pediatric suffix = only the pediatric models of that lineage
#               (e.g. "CNS/Brain+Pediatric").
#               Absent (and no -d file) -> DepMap SKIPPED. No default.
#   -d DEPMAP : none | FILE1[,FILE2,...] (custom files in CRISPRGeneEffect
#               format). Mutually exclusive with -T.
#   -p NPROC  : number of cores to use (default: nproc)
MODE=full
TOOLSEL=all
TUMOR=${TUMOR:-}
DEPMAP=none
NPROC_OPT=${NPROC:-}
while [ $# -gt 0 ]; do
  case "$1" in
    -m) MODE=${2:-}; shift 2 || die "-m requires a value";;
    -t) TOOLSEL=${2:-}; shift 2 || die "-t requires a value";;
    -T) TUMOR=${2:-}; shift 2 || die "-T requires a value";;
    -d) DEPMAP=${2:-}; shift 2 || die "-d requires a value";;
    -p) NPROC_OPT=${2:-}; shift 2 || die "-p requires a value";;
    *) die "unknown option: $1 (use -m/-t/-T/-d/-p)";;
  esac
done
case "$MODE" in
  full|callers|annotate) :;;
  *) die "invalid mode: $MODE (expected full|callers|annotate)";;
esac

# ---------- -T / -d validation (DepMap) ---------------------------------------
# -T: tumor (OncotreeLineage of Model.csv) or MBL (hematologic).
# -d: none | FILE1,FILE2,... (custom files in CRISPRGeneEffect format).
# The two are mutually exclusive; neither -> DepMap SKIPPED
# (there is no longer a "default" CNS).
if [ "$DEPMAP" = "default" ]; then
  die "-d default no longer exists: the DepMap annotation has no default. Use -T <tumor> (e.g. -T Lung, or -T MBL for hematologic) or -d FILE1,FILE2 for custom files. Without -T nor -d the DepMap is skipped."
fi
if [ -n "$TUMOR" ] && [ "$DEPMAP" != "none" ]; then
  die "-T $TUMOR and -d $DEPMAP are mutually exclusive (either the tumor type, or the custom files, not both)"
fi
if [ "$MODE" != "callers" ]; then
  if [ -n "$TUMOR" ]; then
    # With -T the two DepMap CSVs must be mounted (optional for users not using -T).
    [ -s "$T/Nanomatizer/Databsese/Model.csv" ] \
      || die "-T $TUMOR: Model.csv missing in $T/Nanomatizer/Databsese/ (download them with Nanomatizer/download_databases.sh and mount the folder with -v <your_db_folder>:/opt/nanomatizer/userdb)"
    [ -s "$T/Nanomatizer/Databsese/CRISPRGeneEffect.csv" ] \
      || die "-T $TUMOR: CRISPRGeneEffect.csv missing in $T/Nanomatizer/Databsese/ (download them with Nanomatizer/download_databases.sh and mount the folder with -v <your_db_folder>:/opt/nanomatizer/userdb)"
  elif [ "$DEPMAP" != "none" ]; then
    # custom files: they must already exist (in container or host), clear STOP
    IFS=',' read -ra _DFILES <<< "$DEPMAP"
    for _df in "${_DFILES[@]}"; do
      [ -s "$_df" ] || die "-d: custom file missing or empty: $_df"
    done
  fi
fi
# ---------- external databases (not shipped with the image) --------------------
# OncoKB, DECIPHER, TRRUST and DGV are not shipped: the user downloads them once
# with Nanomatizer/download_databases.sh and mounts the folder in $T/userdb
# (the Databsese entries are symlinks into it). Needed by STEP6 only.
if [ "$MODE" != "callers" ]; then
  for _db in OncoKB.tsv population_cnv_grch38.txt trrust_rawdata.human.tsv DGV.GS.hg38.gff3; do
    [ -s "$T/Nanomatizer/Databsese/$_db" ] \
      || die "$_db missing: it is not shipped with the image. Download it once with 'bash Nanomatizer/download_databases.sh' and mount the folder with -v <your_db_folder>:/opt/nanomatizer/userdb (see the tutorial, section 0)"
  done
fi
# RUN_SET: which callers to run.
#   annotate -> none; full/callers -> all (TOOLSEL=all) or the subset
#   derived from -t (token -> caller: sniffles->sniffles, cutesv->cutesv,
#   gladiator->nanogladiator, gasoline/intertrans->gasoline).
RUN_SET="all"
if [ "$MODE" = "annotate" ]; then
  RUN_SET="none"
elif [ "$TOOLSEL" != "all" ]; then
  RUN_SET=""
  IFS=',' read -ra _TOKS <<< "$TOOLSEL"
  for _tok in "${_TOKS[@]}"; do
    case "$_tok" in
      sniffles|sniffles2)   [ -z "${_RUNSET_sniffles:-}" ]  && RUN_SET="${RUN_SET:+$RUN_SET,}sniffles";   _RUNSET_sniffles=1;;
      cutesv)               [ -z "${_RUNSET_cutesv:-}" ]    && RUN_SET="${RUN_SET:+$RUN_SET,}cutesv";     _RUNSET_cutesv=1;;
      gladiator|ng100k)     [ -z "${_RUNSET_nangl:-}" ]     && RUN_SET="${RUN_SET:+$RUN_SET,}nanogladiator"; _RUNSET_nangl=1;;
      gasoline|intertrans)  [ -z "${_RUNSET_gas:-}" ]       && RUN_SET="${RUN_SET:+$RUN_SET,}gasoline";   _RUNSET_gas=1;;
      *) die "invalid tool in -t: $_tok (expected sniffles|cutesv|gasoline|intertrans|gladiator)";;
    esac
  done
  [ -n "$RUN_SET" ] || die "-t produced no caller: $TOOLSEL"
fi
in_set() { case "$1" in all) return 0;; none) return 1;; *) case ",$1," in *",$2,"*) return 0;; esac; return 1;; esac; }

# ---------- one instance at a time --------------------------------------------
exec 9>"$LOCK"
flock -n 9 || die "another pipeline instance is already running (lock: $LOCK)"

# ---------- global preconditions ----------------------------------------------
[ -n "$SAMPLE" ] || die "SAMPLE is required: pass -e SAMPLE=<name> (env SAMPLE) — the image is generic and has no default sample"
export SAMPLE
[ -x "$MM" ]  || die "micromamba missing: $MM"
# BAM/REF are only needed when running callers (full/callers); annotate does not use them.
# The indexes are NOT options: the pipeline checks them and, if missing, builds
# them itself (samtools index / samtools faidx). Idempotent.
ensure_bam_index() {
  local B=$1 I
  [ -f "$B" ] || return 0   # BAM presence is checked by the caller
  I="${B}.bai"
  if [ -s "$I" ]; then
    step "BAM index present: $I"
    return 0
  fi
  step "BAM index missing: building $I (samtools index)"
  # env with samtools guaranteed (sv-sniffles does NOT include it; nanogladiator does)
  if ! $MM run -p "$MAMBA_ROOT_PREFIX/envs/nanogladiator" samtools index "$B" \
       > "$LOGDIR/pipeline_bam_index.log" 2>&1; then
    die "could not build the BAM index $I (see pipeline_bam_index.log — the input file must be sorted and mounted read/write)"
  fi
  [ -s "$I" ] || die "BAM index $I not created (see pipeline_bam_index.log)"
  step "BAM index built: $I"
}
ensure_ref_index() {
  local R=$1 I
  [ -f "$R" ] || return 0
  I="${R}.fai"
  if [ -s "$I" ]; then
    step "REF index present: $I"
    return 0
  fi
  step "REF index missing: building $I (samtools faidx)"
  # env with samtools guaranteed (sv-sniffles does NOT include it; nanogladiator does)
  if ! $MM run -p "$MAMBA_ROOT_PREFIX/envs/nanogladiator" samtools faidx "$R" \
       > "$LOGDIR/pipeline_ref_index.log" 2>&1; then
    die "could not build the REF index $I (see pipeline_ref_index.log — the FASTA must be mounted read/write)"
  fi
  [ -s "$I" ] || die "REF index $I not created (see pipeline_ref_index.log)"
  step "REF index built: $I"
}
if [ "$RUN_SET" != "none" ]; then
  [ -f "$BAM" ] || die "BAM missing: $BAM"
  [ -f "$REF" ] || die "reference missing: $REF (it must be GRCh38/hg38 — the pipeline is GRCh38-only)"
  ensure_bam_index "$BAM"
  ensure_ref_index "$REF"
fi
# Cores to use: -p option (env NPROC) or all available.
if [ -n "$NPROC_OPT" ]; then
  case "$NPROC_OPT" in
    ''|*[!0-9]*) die "-p: invalid number of cores: $NPROC_OPT";;
  esac
  [ "$NPROC_OPT" -ge 1 ] || die "-p: invalid number of cores: $NPROC_OPT (minimum 1)"
  NPROC=$NPROC_OPT
else
  NPROC=$(nproc)
fi
step "PIPELINE START (core=$NPROC) — sample=$SAMPLE mode=$MODE tools=$TOOLSEL tumor=${TUMOR:-skip} depmap=$DEPMAP"

# =============================================================================
# STEP 1/5 — NanoGLADIATOR (run6, R server, 1 Mb windows)
# Marker:  .../R6_W1Mb/Results/<LABEL>/HSLMResults_<LABEL>.txt + 24/24 per-chrom
# Genome fingerprint: $NGOUT/.ng_target = "<fingerprint> <target>"
#   - standard GRCh38 no-alt (fingerprint == STANDARD_FP) -> bundled target
#     hg38w1000000 (never regenerated, BRIEF);
#   - non-standard genome -> 1 Mb windows built from the mounted genome
#     (build_ng_target.sh, cache $OUT/.ng_targets/<name>) + injection into the
#     tool via SYMLINK: a runtime artifact recreated on every run, no tool
#     file is written (BRIEF 1).
# =============================================================================
NGD=$T/NanoGLADIATOR_1.0
TARGET_STD=hg38w1000000
TARGET=hg38w1000000
WDIR=R6_W1Mb
NGL_LABEL="${SAMPLE}w1mb_germinal"
SRCW=$OUT/NanoGLADIATOR/$SAMPLE/W1Mb            # RC/RCNorm prepare (read-only)
NGOUT=$OUT/NanoGLADIATOR/$SAMPLE/$WDIR
NGMARKER=$NGOUT/Results/${NGL_LABEL}/HSLMResults_${NGL_LABEL}.txt
NGTGTMARKER=$NGOUT/.ng_target
# Genome fingerprint: md5 of the "name length" lines of the 24 primary
# chromosomes (chr1..chr22, chrX, chrY) of the .fai, sorted. A constant,
# verified value for the standard genome used by the pipeline (GRCh38 no-alt,
# 24 primaries with the chr prefix): identical computed from the .fai and from
# the @SQ of the BAM.
STANDARD_FP="1475bd4e081404f7aee6a7bc6e314c3a"
genome_fp() {
  awk '{print $1, $2}' "$1" | grep -E '^(chr([1-9]|1[0-9]|2[0-2])|chrX|chrY) [0-9]+$' | LC_ALL=C sort | md5sum | awk '{print $1}'
}

# NanoGLADIATOR 1 Mb windows target:
#   - standard genome (fingerprint == STANDARD_FP) -> bundled target
#     hg38w1000000 already in the image: NEVER regenerated (BRIEF).
#   - non-standard genome -> 1 Mb windows built from the mounted genome with
#     build_ng_target.sh (cache in $OUT/.ng_targets/<name>, idempotent) and
#     injected into the tool via SYMLINK: a runtime artifact recreated on
#     every run, no tool file is written or modified (BRIEF 1).
# The target is a function ONLY of the genome (the windows depend on the .fai,
# the GC on the FASTA, the mappability is constant 1.0 in the builder), so the
# cache is reused on every run: the symlink is the only recreated element.
FP=""
if [ -f "${REF}.fai" ]; then
  FP=$(genome_fp "${REF}.fai")
  if [ "$FP" = "$STANDARD_FP" ]; then
    TARGET=$TARGET_STD
  else
    TARGET="w1000000_${FP:0:8}"
    BUILT=$OUT/.ng_targets/$TARGET
    step "[1/5] NanoGLADIATOR: NON-standard genome (fp=$FP) — 1 Mb windows to build: $BUILT"
    $MM run -p "$MAMBA_ROOT_PREFIX/envs/nanogladiator" \
        bash "$T/Nanomatizer/Scripts/NewScripts/build_ng_target.sh" \
        --fasta "$REF" --fai "${REF}.fai" --name "$TARGET" --out "$OUT/.ng_targets" \
        > "$LOGDIR/pipeline_ng_target_build.log" 2>&1
    rc=$?
    [ $rc -eq 0 ] || die "NanoGLADIATOR target build rc=$rc (see pipeline_ng_target_build.log)"
    [ -s "$BUILT/TargetWindow.bed" ] || die "target build: TargetWindow.bed missing in $BUILT"
    [ -d "$NGD/data/targets/hg38" ] || die "$NGD/data/targets/hg38 missing (tool incomplete?)"
    # symlink injection: runtime artifact, recreated on every run
    ln -sfn "$BUILT" "$NGD/data/targets/hg38/$TARGET"
    [ -d "$NGD/data/targets/hg38/$TARGET" ] || die "target symlink not working: $NGD/data/targets/hg38/$TARGET"
    step "[1/5] NanoGLADIATOR: non-standard target ready ($TARGET, $(wc -l < "$BUILT/TargetWindow.bed") windows) — symlink $NGD/data/targets/hg38/$TARGET"
  fi
fi

if ! in_set "$RUN_SET" nanogladiator; then
  step "[1/5] NanoGLADIATOR SKIPPED: not requested (mode=$MODE)"
else
step "[1/5] NanoGLADIATOR: status check"
if [ -f "$NGMARKER" ]; then
    # genome coherence: the existing run must have been produced with the
    # current genome (fingerprint registered in $NGTGTMARKER at completion)
    if [ -n "$FP" ]; then
        [ -f "$NGTGTMARKER" ] \
            || die "NanoGLADIATOR: marker present but $NGTGTMARKER missing (a previous run before the pipeline with fingerprint) — inconsistent state, I remove nothing, a decision is required"
        _MFP=$(awk '{print $1; exit}' "$NGTGTMARKER")
        _MTGT=$(awk '{print $2; exit}' "$NGTGTMARKER")
        [ "$_MFP" = "$FP" ] && [ "$_MTGT" = "$TARGET" ] \
            || die "NanoGLADIATOR: marker present but genome INCONSISTENT (registered '$_MFP $_MTGT', current '$FP $TARGET') — I remove nothing, a decision is required"
    fi
    NCHR=$(find "$NGOUT" -path "*Results*" -name "chr*_FractionCallResults_*.txt" 2>/dev/null | wc -l)
    if [ "$NCHR" -eq 24 ]; then
        step "[1/5] NanoGLADIATOR SKIPPED: marker present, per-chrom 24/24, coherent genome (run already valid)"
    else
        die "NanoGLADIATOR in INCONSISTENT state: marker present but per-chrom $NCHR/24 (I remove nothing, a decision is required)"
    fi
else
    [ -f "$NGD/lib/F77/F4R.so" ] || die "NanoGLADIATOR: F4R.so missing in lib/F77"
    [ -f "$NGD/lib/F77/FastJointSLMLibrary.so" ] || die "NanoGLADIATOR: FastJointSLMLibrary.so missing in lib/F77"
    grep -q "TotalPredBreak <- JointSeg(DataSeq,eta,omega,muk,mi,smu,sepsilon)" \
        "$NGD/lib/R/NanoGLADIATORInferenceGenome.R" \
        || die "NanoGLADIATORInferenceGenome.R is NOT the server version (whole chromosome)"
    NGT=10   # processors for NanoGLADIATOR (<= nproc, same as the run6 launcher)
    [ "$NGT" -le "$NPROC" ] || NGT=$NPROC
    # RC/RCNorm are the products of the NanoGLADIATOR PREPARE in $SRCW (W1Mb).
    # Sample already prepared (e.g. validation): they are already present and we
    # copy them. NEW sample: missing -> we run the PREPARE here in $SRCW
    # (~35-40 min at 10 processors) and then continue with the normal flow
    # (copy + analysis).
    # Partial state (some RData but not 24/24) = inconsistent -> STOP, touch nothing.
    RC_N=$(ls "$SRCW/RC/"*.RData 2>/dev/null | wc -l)
    if [ "$RC_N" -gt 0 ] && [ "$RC_N" -lt 24 ]; then
        die "NanoGLADIATOR in INCONSISTENT state: RC in $SRCW/RC = $RC_N/24 (I remove nothing, a decision is required)"
    fi
    if [ "$RC_N" -lt 24 ] || [ ! -f "$SRCW/RCNorm/${NGL_LABEL}.NRC.RData" ]; then
        [ -f "$NGD/data/centromere/CentromerePosition_hg38.txt" ] \
            || die "NanoGLADIATOR: hg38 centromere file missing in $NGD/data/centromere/CentromerePosition_hg38.txt"
        mkdir -p "$SRCW"
        step "[1/5] NanoGLADIATOR PREPARE RUNNING (new sample, $TARGET, $NGT processors) — ~35-40 min"
        PREP="$SRCW/${SAMPLE}.w1mb.experimental.file.txt"
        printf '%s %s %s\n' "$BAM" "$SRCW" "$NGL_LABEL" > "$PREP"
        s=$(date +%s)
        ( cd "$NGD" && $MM run -p "$MAMBA_ROOT_PREFIX/envs/nanogladiator" \
            perl NanoGLADIATORPrepare.pl "$PREP" --processors $NGT \
            --target "$TARGET" --assembly hg38 ) \
            > "$LOGDIR/pipeline_nanogladiator_prepare.log" 2>&1
        rc=$?
        e=$(($(date +%s) - s))
        step "[1/5] NanoGLADIATOR prepare rc=$rc elapsed=$((e/60))min"
        [ $rc -eq 0 ] || die "NanoGLADIATOR prepare rc=$rc (see pipeline_nanogladiator_prepare.log)"
        RC_N=$(ls "$SRCW/RC/"*.RData 2>/dev/null | wc -l)
        [ "$RC_N" -eq 24 ] || die "NanoGLADIATOR prepare finished but RC in $SRCW/RC = $RC_N/24"
        step "[1/5] NanoGLADIATOR prepare COMPLETED (RC 24/24, $((e/60)) min)"
    fi
    [ -f "$SRCW/RCNorm/${NGL_LABEL}.NRC.RData" ] || die "RCNorm ${NGL_LABEL}.NRC.RData missing in $SRCW/RCNorm"
    mkdir -p "$NGOUT"
    [ -d "$NGOUT/RC" ]     || { mkdir -p "$NGOUT/RC";     cp -p "$SRCW/RC/"*.RData     "$NGOUT/RC/"; }
    [ -d "$NGOUT/RCNorm" ] || { mkdir -p "$NGOUT/RCNorm"; cp -p "$SRCW/RCNorm/"*.RData "$NGOUT/RCNorm/"; }
    NRC=$(ls "$NGOUT/RC"/*.RData 2>/dev/null | wc -l)
    [ "$NRC" -eq 24 ] || die "RC RData in $NGOUT/RC = $NRC/24"

    step "[1/5] NanoGLADIATOR RUNNING (R server, $TARGET, $NGT processors)"
    ANA="$NGOUT/${SAMPLE}.w1mb.experimental.file.analysis.txt"
    printf '%s %s %s\n' "T1" "$NGOUT" "$NGL_LABEL" > "$ANA"
    s=$(date +%s)
    # Env identical to the run6 launcher: HOME/XDG on dedicated folders + PATH
    # with the nanogladiator conda env (so the bare 'R' invoked by
    # NanoGLADIATORAnalysis.pl uses the conda R 4.5.3 with Hmisc/skellam,
    # in the container and on the host).
    ( cd "$NGD" && HOME="$T/fakehome" XDG_CACHE_HOME="$T/cache" MAMBA_ROOT_PREFIX="$MAMBA_ROOT_PREFIX" \
        PATH="$T/fakehome/bin:$MAMBA_ROOT_PREFIX/envs/nanogladiator/bin:$PATH" \
        perl NanoGLADIATORAnalysis.pl "$ANA" --processors $NGT \
        --target "$TARGET" --assembly hg38 --output "$NGOUT" --mode nocontrol ) \
        > "$LOGDIR/pipeline_nanogladiator.log" 2>&1
    rc=$?
    e=$(($(date +%s) - s))
    step "[1/5] NanoGLADIATOR rc=$rc elapsed=$((e/60))min"
    if grep -q "number of items to replace is not a multiple of replacement length" \
         "$LOGDIR/pipeline_nanogladiator.log" 2>/dev/null; then
        NCHR=$(find "$NGOUT" -name "chr*_FractionCallResults_*.txt" 2>/dev/null | wc -l)
        die "NanoGLADIATOR MatrixSeg crash — per-chrom $NCHR/24"
    fi
    [ $rc -eq 0 ] && [ -f "$NGMARKER" ] || die "NanoGLADIATOR rc=$rc, marker $([ -f "$NGMARKER" ] && echo present || echo MISSING)"
    NCHR=$(find "$NGOUT" -path "*Results*" -name "chr*_FractionCallResults_*.txt" 2>/dev/null | wc -l)
    [ "$NCHR" -eq 24 ] || die "NanoGLADIATOR completed but per-chrom $NCHR/24"
    # register the genome fingerprint + target used: needed by the coherence
    # check of later runs (marker present + changed genome = inconsistent)
    echo "$FP $TARGET" > "$NGTGTMARKER"
    step "[1/5] NanoGLADIATOR COMPLETED (per-chrom 24/24, $((e/60)) min) — target $TARGET fp=${FP:-(fai absent)}"
fi
fi

# =============================================================================
# STEP 2/5 — Sniffles2
# Marker: VCF present (the "signature" check on the SV count was removed: it is
#         sample-specific, no validation-sample thresholds apply)
# =============================================================================
SNDIR=$OUT/Sniffles2/$SAMPLE
SVCF=$SNDIR/${SAMPLE}_no_alt.srt.vcf
S_T=10; [ "$S_T" -le "$NPROC" ] || S_T=$NPROC

sv_count() { grep -vc '^#' "$1" 2>/dev/null || echo 0; }

if ! in_set "$RUN_SET" sniffles; then
  step "[2/5] Sniffles SKIPPED: not requested (mode=$MODE)"
else
step "[2/5] Sniffles: status check"
if [ -f "$SVCF" ] && [ -s "$SVCF" ]; then
    NSV=$(sv_count "$SVCF")
    step "[2/5] Sniffles SKIPPED: VCF present ($NSV SV) — skip"
else
    step "[2/5] Sniffles RUNNING ($S_T threads)"
    mkdir -p "$SNDIR"
    s=$(date +%s)
    $MM run -p "$MAMBA_ROOT_PREFIX/envs/sv-sniffles" sniffles \
        -i "$BAM" -v "$SVCF" --threads $S_T --reference "$REF" \
        > "$LOGDIR/pipeline_sniffles.log" 2>&1
    rc=$?
    e=$(($(date +%s) - s))
    [ $rc -eq 0 ] && [ -f "$SVCF" ] || die "Sniffles rc=$rc, VCF $([ -f "$SVCF" ] && echo present || echo MISSING) (see pipeline_sniffles.log)"
    NSV=$(sv_count "$SVCF")
    [ "$NSV" -gt 0 ] || die "Sniffles: VCF present but 0 SV (see pipeline_sniffles.log)"
    step "[2/5] Sniffles COMPLETED ($NSV SV, $((e/60)) min)"
fi
fi

# =============================================================================
# STEP 3/5 — cuteSV
# Marker: VCF present ; RAM >= 4 GB before (SV-count signature removed)
# =============================================================================
CSDIR=$OUT/CuteSV/$SAMPLE
CVCF=$CSDIR/${SAMPLE}cutesv_genotype.vcf
C_T=20; [ "$C_T" -le "$NPROC" ] || C_T=$NPROC

if ! in_set "$RUN_SET" cutesv; then
  step "[3/5] cuteSV SKIPPED: not requested (mode=$MODE)"
else
step "[3/5] cuteSV: status check"
if [ -f "$CVCF" ] && [ -s "$CVCF" ]; then
    NSV=$(sv_count "$CVCF")
    step "[3/5] cuteSV SKIPPED: VCF present ($NSV SV) — skip"
else
    freemem=$(free -g | awk 'NR==2{print $7}')
    step "free RAM before cuteSV: ${freemem} GB"
    [ "$freemem" -ge 4 ] || die "free RAM ${freemem} GB < 4 GB — not starting cuteSV"
    step "[3/5] cuteSV RUNNING ($C_T threads)"
    mkdir -p "$CSDIR"
    s=$(date +%s)
    $MM run -p "$MAMBA_ROOT_PREFIX/envs/sv-cutesv" cuteSV \
        --threads $C_T \
        --max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
        "$BAM" "$REF" "$CVCF" "$CSDIR/" \
        > "$LOGDIR/pipeline_cutesv.log" 2>&1
    rc=$?
    e=$(($(date +%s) - s))
    [ $rc -eq 0 ] && [ -f "$CVCF" ] || die "cuteSV rc=$rc, VCF $([ -f "$CVCF" ] && echo present || echo MISSING) (see pipeline_cutesv.log)"
    NSV=$(sv_count "$CVCF")
    [ "$NSV" -gt 0 ] || die "cuteSV: VCF present but 0 SV (see pipeline_cutesv.log)"
    step "[3/5] cuteSV COMPLETED ($NSV SV, $((e/60)) min)"
fi
fi

# =============================================================================
# STEP 4/5 — Gasoline 2.0 (Parse + Germline + copy)
# Marker: Parse/<SAMPLE>/<SAMPLE>.Results.tsv + Trans.tsv present and header
#         correct. The SV-count signatures of the validation sample have been removed
#         (an internal validation control of the reference sample, not
#         applicable to other BAMs): only the structure remains here.
# Note: the Gasoline/Germline/ folder is not used by the annotation (which
# reads Parse/); the copy is reproduced only for parity with the protocol.
# If the Parse dir exists but is NOT valid -> normally STOP (cleaning = the
# user's job), EXCEPT when the parse itself completed (ParsedSplit present):
# that is a resume point, and only the GermlineDetector is re-run.
# =============================================================================
GAS=$T/GASOLINE_2.0
GPARSE=$OUT/Gasoline/Parse/$SAMPLE
GGERM=$OUT/Gasoline/Germline/$SAMPLE
GREZ=$GPARSE/${SAMPLE}.Results.tsv
GTRA=$GPARSE/${SAMPLE}.Results.Trans.tsv
GANO=$GPARSE/${SAMPLE}.Results.Anno.tsv

gasoline_valid() {
    # $1 = Parse dir; returns 0 if the structure is valid, 1 otherwise
    # (only files + header: no signature on the counts, generic per sample)
    local D=$1 F
    F=$D/${SAMPLE}.Results.tsv
    [ -s "$F" ] || return 1
    [ -s "$D/${SAMPLE}.Results.Trans.tsv" ] || return 1
    # column 19 = "SV Type" (Gasoline 2.0 layout, sample-independent)
    head -1 "$F" 2>/dev/null | awk -F'\t' '{print $19}' | grep -qx "SV Type" || return 1
    return 0
}

if ! in_set "$RUN_SET" gasoline; then
  step "[4/5] Gasoline SKIPPED: not requested (mode=$MODE)"
else
step "[4/5] Gasoline: status check"
PARSPL=$GPARSE/AnalysisSupp/${SAMPLE}.ParsedSplit.txt
# NOTE: GASOLINEParse.pl only extracts the SV signatures into AnalysisSupp/
# (ParsedSplit + AlignedIns/Del files); Results.tsv is written LATER by the
# GermlineDetector. So:
#   - parse marker  = AnalysisSupp/<SAMPLE>.ParsedSplit.txt present and non-empty
#   - step marker   = <SAMPLE>.Results.tsv + Trans.tsv (gasoline_valid)
# A dir with a completed parse but no Results.tsv is a valid RESUME point
# (the detector runs straight on the ParsedSplit signatures).
if [ -e "$GPARSE" ]; then
    if gasoline_valid "$GPARSE" 2>/dev/null; then
        NROWS=$(($(wc -l < "$GREZ") - 1))
        step "[4/5] Gasoline SKIPPED: Parse valid ($NROWS rows in Results.tsv, signature OK)"
    elif [ -s "$PARSPL" ]; then
        step "[4/5] Gasoline RESUMING: Parse signatures present ($PARSPL), running GermlineDetector only (no re-parse)"
    else
        die "Gasoline Parse/$SAMPLE exists but is NOT valid/incomplete — STOP. Cleaning the dir is the user's job (BRIEF 2a): I remove nothing"
    fi
fi
if ! gasoline_valid "$GPARSE" 2>/dev/null; then
    [ -f "$GAS/GASOLINEParse.pl" ] || die "GASOLINE_2.0: GASOLINEParse.pl missing in $GAS"
    [ -f "$GAS/GASOLINEGermlineDetector.pl" ] || die "GASOLINE_2.0: GASOLINEGermlineDetector.pl missing in $GAS"
    GT=15   # -NT for the germline detector (same as the protocol)
    [ "$GT" -le "$NPROC" ] || GT=$NPROC
    # mkdir -p: no-op when the dirs already exist (resume case)
    mkdir -p "$GPARSE" "$GGERM/AnalysisSupp"
    s=$(date +%s)
    # GASOLINEParse calls samtools and GASOLINEGermlineDetector calls
    # `R --slave` (needs the igraph package): the container's base PATH does
    # NOT contain conda envs, so I build it explicitly. Order
    # annot:nanogladiator (verified in the image):
    #   R        -> nanomatizer-annot (has igraph; nanogladiator's R does NOT)
    #   samtools -> nanogladiator 1.24 (identical to sv-gasoline; annot lacks it)
    GPATH="$T/fakehome/bin:$MAMBA_ROOT_PREFIX/envs/nanomatizer-annot/bin:$MAMBA_ROOT_PREFIX/envs/nanogladiator/bin:$PATH"
    PATH="$GPATH" command -v samtools >/dev/null 2>&1 \
        || die "Gasoline: samtools not available ($MAMBA_ROOT_PREFIX/envs/nanogladiator/bin required) — STOP"
    PATH="$GPATH" R --vanilla -e 'requireNamespace("igraph", quietly=TRUE)' >/dev/null 2>&1 \
        || die "Gasoline: R without the igraph package ($MAMBA_ROOT_PREFIX/envs/nanomatizer-annot/bin required) — STOP"
    if [ -s "$PARSPL" ]; then
        step "[4/5] Gasoline: Parse SKIPPED (signatures already extracted)"
    else
        step "[4/5] Gasoline Parse RUNNING"
        ( cd "$GAS" && PATH="$GPATH" perl GASOLINEParse.pl -B "$BAM" -L "$SAMPLE" -O "$GPARSE" -R hg38 ) \
            > "$LOGDIR/pipeline_gasoline_parse.log" 2>&1
        rc=$?
        [ $rc -eq 0 ] && [ -s "$PARSPL" ] || die "Gasoline Parse rc=$rc, ParsedSplit $([ -s "$PARSPL" ] && echo present || echo MISSING) (see pipeline_gasoline_parse.log)"
    fi
    step "[4/5] Gasoline GermlineDetector RUNNING (-NT $GT)"
    ( cd "$GAS" && PATH="$GPATH" perl GASOLINEGermlineDetector.pl -B "$BAM" -L "$SAMPLE" -O "$GPARSE" -NT $GT ) \
        > "$LOGDIR/pipeline_gasoline_germline.log" 2>&1
    rc=$?
    [ $rc -eq 0 ] && [ -s "$GREZ" ] || die "Gasoline GermlineDetector rc=$rc, Results.tsv $([ -s "$GREZ" ] && echo present || echo MISSING) (see pipeline_gasoline_germline.log)"
    # NOTE: Gasoline 2.0 exits 0 even on partial failures -> the signature is decisive
    # NB: 'object TableGenCodeGASOLINEOverlap4Saving not found' = dirty input
    if grep -q "TableGenCodeGASOLINEOverlap4Saving" "$LOGDIR/pipeline_gasoline_germline.log" 2>/dev/null; then
        die "Gasoline: error 'object TableGenCodeGASOLINEOverlap4Saving not found' = dirty input"
    fi
    # copy Parse -> Germline (parity with the protocol; never overwrite valid data)
    if [ -e "$GGERM/${SAMPLE}.Results.tsv" ]; then
        step "[4/5] Gasoline: Germline/ already populated, skip the Parse->Germline copy"
    else
        cp -r "$GPARSE"/. "$GGERM"/ || die "Gasoline: Parse->Germline copy failed"
    fi
    e=$(($(date +%s) - s))
    gasoline_valid "$GPARSE" || die "Gasoline finished but signature NOT valid (see pipeline.log)"
    NROWS=$(($(wc -l < "$GREZ") - 1))
    step "[4/5] Gasoline COMPLETED ($NROWS rows in Results.tsv, $((e/60)) min)"
fi
fi

# =============================================================================
# STEP 5/5 — Annotations + aggregation (module v2)
# =============================================================================
# callers mode -> the callers are already done here; no annotation.
if [ "$MODE" = "callers" ]; then
  step "PIPELINE COMPLETED (mode=callers, tools=$TOOLSEL) — only callers run, no annotation/aggregation"
  echo "PIPELINE_EXIT_OK"
  exit 0
fi

MV2=$T/Nanomatizer/Scripts/NewScripts/step6/module_v2
# STEP6 mode and tool selection based on the pipeline mode:
#   full     -> all + TOOLSEL (TOOLSEL=all == original behavior, identical)
#   annotate -> all + TOOLSEL  (annotation+aggregation on the results already in OUT)
case "$MODE" in
  full)     STEP6_MODE=all; STEP6_TOOLS=$TOOLSEL;;
  annotate) STEP6_MODE=all; STEP6_TOOLS=$TOOLSEL;;
esac

step "[5/5] Annotation+aggregation (mode=$STEP6_MODE, tools=$STEP6_TOOLS, tumor=${TUMOR:-skip}, depmap=$DEPMAP, outdir=$OUT)"
bash "$MV2/run_module.sh" "$STEP6_MODE" "$OUT" "$STEP6_TOOLS" "$DEPMAP" "$TUMOR" >> "$PLLOG" 2>&1
rc=$?
[ $rc -eq 0 ] || die "STEP5 annotations rc=$rc (see $PLLOG and $OUT/step6_v2.log)"
[ -s "$OUT/Nanomatizer_${SAMPLE}_aggregated.xlsx" ] || die "STEP5: aggregate missing $OUT/Nanomatizer_${SAMPLE}_aggregated.xlsx"
step "PIPELINE COMPLETED — aggregate: $OUT/Nanomatizer_${SAMPLE}_aggregated.xlsx (mode=$MODE, tools=$STEP6_TOOLS, tumor=${TUMOR:-skip}, depmap=$DEPMAP)"
echo "PIPELINE_EXIT_OK"
exit 0
