#!/bin/bash
# =============================================================================
# run_module.sh — wrapper per il modulo annotazioni v2 (STEP 6)
# ----------------------------------------------------------------------------
# Esegue nanomatizer_annotation_v2.r sotto env nanomatizer-annot (stesso
# incapsulamento del wrapper originale step6/run_step6_annotations.sh).
#
# Uso:
#   bash run_module.sh [MODE] [OUTDIR] [TOOLS] [DEPMAP] [TUMOR]
#     MODE   : all | annotation | aggregation   (default: all)
#     OUTDIR : cartella output (default: $T/Nanomatizer/Results)
#     TOOLS  : all | elenco (sniffles|cutesv|gasoline|intertrans|gladiator)
#              (default: all). In "annotation" esegue solo i tool scelti;
#              in "aggregation" aggrega solo i tool scelti (i CSV dei tool
#              esclusi, se presenti, vengono azzerati; se assenti -> placeholder).
#     DEPMAP : none (DEFAULT: DepMap skipped — there is no longer a CNS default)
#              | FILE[,FILE..] (custom files, CRISPRGeneEffect format)
#              | default (technical back-door: 3 CNS/Pediatric/MBL groups,
#                not exposed by the pipeline CLI)
#     TUMOR  : tumor type (OncotreeLineage of Model.csv, e.g. "Lung", or
#              "MBL" for hematologic). +Pediatric suffix = only pediatric
#              models (AgeCategory 'Pediatric') of that lineage, e.g.
#              "CNS/Brain+Pediatric". Generates ONE DepMap group from the two
#              CSVs in $T/Nanomatizer/Databsese (Model.csv +
#              CRISPRGeneEffect.csv, the user's release — the files stay
#              optional). Mutually exclusive with the custom DEPMAP files.
#              Default: empty (skipped).
#
#     SAMPLE : REQUIRED via env — sample name used in the per-tool paths and
#              in the aggregate name. No default: the module is generic and
#              the sample is the user's.
#
# Examples (SAMPLE always via env):
#   SAMPLE=mionome bash run_module.sh                    # full run
#   SAMPLE=mionome bash run_module.sh aggregation $OUTDIR
#                                                       # aggregation only on already-annotated data
#   SAMPLE=mionome bash run_module.sh annotation $OUTDIR gladiator
#
# Rules: STOP at the first error, STATUS for every event. Removes nothing.
# The OUTDIR folder may also differ from the production one (use: validation
# on a scratch dir without touching the validated outputs).
# =============================================================================
# Layout generico: T = 6 livelli sopra questa cartella
# (…/Nanomatizer/Scripts/NewScripts/step6/module_v2).
# In immagine: /opt/nanomatizer. Nessun path utente hardcodato.
T="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../../../.." && pwd)"
MV2=$T/Nanomatizer/Scripts/NewScripts/step6/module_v2
MODE=${1:-all}
OUTDIR=${2:-$T/Nanomatizer/Results}
TOOLS=${3:-all}
DEPMAP=${4:-none}
TUMOR=${5:-}
SAMPLE=${SAMPLE:-}
LOG=$OUTDIR/step6_v2.log
STATUS=$T/Nanomatizer/Results/STATUS.txt
GFF3=$T/Nanomatizer/Databsese/EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3
ST(){ echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" >>$STATUS; }

mkdir -p "$OUTDIR"

# --- PRECONDITION: SAMPLE required (env) --------------------------------------
# The module is generic: there is no default sample. The pipeline already
# exports SAMPLE; here it stops if it is empty (no silent default).
if [ -z "$SAMPLE" ]; then
  MSG="STEP6-v2 STOP: SAMPLE is required (export SAMPLE=<name> before launching)"
  ST "STOP SAMPLE: $MSG"; exit 3
fi

# --- PRECONDITION: GFF3 v112 present (the 24h watcher is in the original wrapper) ---
if [ ! -s "$GFF3" ] || ! head -1 "$GFF3" | grep -qE "^#+gff-version"; then
  MSG="STEP6-v2 STOP: GFF3 v112 missing or invalid: $GFF3"
  ST "STOP GFF3 v2: $MSG"; exit 3
fi

# --- PRECONDITION: per-tool data (if mode != aggregation) --------------------
# Requires the inputs ONLY of the selected tools. In "annotation" with
# TOOLS=all (default) it is identical to the previous behavior (all 5 required).
# In "aggregation" no per-tool input is needed (the already-produced CSVs are read).
if [ "$MODE" = "all" ] || [ "$MODE" = "annotation" ]; then
  has_tool(){ case ",$1," in *",$2,"*) return 0;; *) return 1;; esac; }
  want_all=0; case "$TOOLS" in all) want_all=1;; esac
  req_list=()
  if [ "$want_all" -eq 1 ] || has_tool "$TOOLS" cutesv;     then req_list+=("$OUTDIR/CuteSV/$SAMPLE/"*.vcf); fi
  if [ "$want_all" -eq 1 ] || has_tool "$TOOLS" gasoline;   then req_list+=("$OUTDIR/Gasoline/Parse/$SAMPLE/$SAMPLE.Results.tsv"); fi
  if [ "$want_all" -eq 1 ] || has_tool "$TOOLS" intertrans; then req_list+=("$OUTDIR/Gasoline/Parse/$SAMPLE/$SAMPLE.Results.Trans.tsv"); fi
  if [ "$want_all" -eq 1 ] || has_tool "$TOOLS" sniffles;   then req_list+=("$OUTDIR/Sniffles2/$SAMPLE/"*.vcf); fi
  if [ "$want_all" -eq 1 ] || has_tool "$TOOLS" gladiator;  then req_list+=("$OUTDIR/NanoGLADIATOR/$SAMPLE/R6_W1Mb/Results/"); fi
  for REQ in "${req_list[@]}"; do
    ls $REQ >/dev/null 2>&1 || { ST "STOP missing input: $REQ"; exit 3; }
  done
fi

# --- R v2 DRIVER EXECUTION ---------------------------------------------------
mkdir -p "$T/Nanomatizer/Results/logs"   # STATUS log (fixed path)
ST "STEP6-v2 START mode=$MODE tools=$TOOLS depmap=$DEPMAP tumor=${TUMOR:-skip} outdir=$OUTDIR"
START=$(date +%s)
HOME=$T/fakehome XDG_CACHE_HOME=$T/cache MAMBA_ROOT_PREFIX=$T/mamba-root \
PATH="$T/fakehome/bin:$T/mamba-root/envs/nanomatizer-annot/bin:$PATH" \
"$T/mamba-root/envs/nanomatizer-annot/bin/Rscript" --vanilla "$MV2/nanomatizer_annotation_v2.r" \
  "sample=$SAMPLE" "outdir=$OUTDIR" "dbdir=$T/Nanomatizer/Databsese" \
  "mode=$MODE" "tools=$TOOLS" "depmap=$DEPMAP" "tumor=$TUMOR" >>"$LOG" 2>&1
RC=$?
ELAPSED=$(( ( $(date +%s) - START ) / 60 ))
ST "STEP6-v2 Rscript finished rc=$RC after ${ELAPSED} min"

if [ $RC -ne 0 ]; then
  ST "STEP6-v2 ERROR rc=$RC mode=$MODE (see $LOG)"
  exit $RC
fi

# --- AGGREGATE CHECK ---------------------------------------------------------
XLSX=$OUTDIR/Nanomatizer_${SAMPLE}_aggregated.xlsx
if [ "$MODE" = "all" ] || [ "$MODE" = "aggregation" ]; then
  if [ ! -s "$XLSX" ]; then
    ST "STEP6-v2 CHECK FAILED: aggregate missing: $XLSX"
    exit 4
  fi
fi

LASTLINE=$(grep "STEP6-AGGREGATO" "$LOG" | tail -1)
ST "STEP6-v2 COMPLETED in ${ELAPSED} min (mode=$MODE) — ${LASTLINE:-}"
echo "STEP6v2_EXIT_OK"
exit 0
