#!/bin/bash
# =============================================================================
# download_databases.sh — fetch the annotation databases that are NOT shipped
# with the nanoCASTLE image (their terms of use do not allow us to
# redistribute them, or we prefer you to get them from the original source).
#
#   REQUIRED by the annotation step (-m full / -m annotate):
#     OncoKB.tsv                 OncoKB Cancer Gene List   (https://www.oncokb.org)
#     population_cnv_grch38.txt  DECIPHER population CNVs  (https://www.deciphergenomics.org)
#     trrust_rawdata.human.tsv   TRRUST v2                 (https://www.grnpedia.org/trrust)
#     DGV.GS.hg38.gff3           DGV Gold Standard, hg38   (http://dgv.tcag.ca)
#   OPTIONAL (only for the DepMap annotation, option -T):
#     CRISPRGeneEffect.csv       DepMap Public (Broad Institute, CC BY 4.0)
#     Model.csv                  DepMap Public (Broad Institute, CC BY 4.0)
#
# By downloading these files you accept the terms of use of each provider.
#
# Usage:
#   bash download_databases.sh [-o OUTDIR] [--no-depmap] [--depmap-24q4] [--force]
#
#   -o OUTDIR    destination folder (default: /opt/nanomatizer/userdb when run
#                inside the container, ./nanocastle_db otherwise)
#   --no-depmap  skip the two DepMap CSVs (~430 MB)
#   --depmap-24q4  do not look for the latest DepMap release, use the pinned 24Q4
#   --force      re-download files that are already present
#
# The folder is then mounted in the container at /opt/nanomatizer/userdb:
#   -v ~/nanocastle_db:/opt/nanomatizer/userdb
#
# It works on the host (needs curl, gzip, md5sum) or inside the image:
#   docker run --rm --user $(id -u):$(id -g) \
#     -v ~/nanocastle_db:/opt/nanomatizer/userdb \
#     gmatteigmattei/nanocastle:1.0.0 bash Nanomatizer/download_databases.sh
#
# DepMap release: the script asks the DepMap portal for the file list and takes
# the most recent "DepMap Public" release. The portal may answer with a human
# verification page instead of the list: in that case the script falls back to
# 24Q4, the last release with a stable public URL (figshare
# 10.25452/figshare.plus.27993248). The release actually downloaded is written
# in OUTDIR/DEPMAP_RELEASE.txt. To use a specific release:
#   - export DEPMAP_CRISPR_URL=... DEPMAP_MODEL_URL=... [DEPMAP_CRISPR_MD5=... DEPMAP_MODEL_MD5=...]
#   - or download CRISPRGeneEffect.csv + Model.csv from
#     https://depmap.org/portal/data_page/?tab=allData and copy them into OUTDIR.
# =============================================================================
set -euo pipefail

ONCOKB_URL=${ONCOKB_URL:-https://www.oncokb.org/api/v1/utils/cancerGeneList.txt}
DECIPHER_URL=${DECIPHER_URL:-https://www.deciphergenomics.org/files/downloads/population_cnv_grch38.txt.gz}
TRRUST_URL=${TRRUST_URL:-https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv}
DGV_URL=${DGV_URL:-http://dgv.tcag.ca/dgv/docs/DGV.GS.hg38.gff3}
# md5 of the versions the pipeline was validated with (set to empty to skip the check)
TRRUST_MD5=${TRRUST_MD5-63bed867f32fa6519a4ee481bcee0785}
DGV_MD5=${DGV_MD5-c70d30e944d396e533c71e03b67442c1}

# DepMap: pinned fallback release (24Q4, figshare) and portal file list
DEPMAP_FALLBACK_RELEASE="DepMap Public 24Q4"
DEPMAP_FALLBACK_CRISPR_URL=https://ndownloader.figshare.com/files/51064667
DEPMAP_FALLBACK_CRISPR_MD5=6edf7ade09b9b34199210b559d4745d3
DEPMAP_FALLBACK_MODEL_URL=https://ndownloader.figshare.com/files/51065297
DEPMAP_FALLBACK_MODEL_MD5=675210d17675f3517b0ce39a3c274f16
DEPMAP_LIST_URL=${DEPMAP_LIST_URL:-https://depmap.org/portal/api/download/files}
# explicit URLs given by the user win over everything (md5 optional)
DEPMAP_USER_URLS=0
if [ -n "${DEPMAP_CRISPR_URL:-}" ] || [ -n "${DEPMAP_MODEL_URL:-}" ]; then
  [ -n "${DEPMAP_CRISPR_URL:-}" ] && [ -n "${DEPMAP_MODEL_URL:-}" ] \
    || { echo "ERROR: set both DEPMAP_CRISPR_URL and DEPMAP_MODEL_URL" >&2; exit 2; }
  DEPMAP_USER_URLS=1
  DEPMAP_RELEASE="user-provided URLs"
  DEPMAP_CRISPR_MD5=${DEPMAP_CRISPR_MD5:-}
  DEPMAP_MODEL_MD5=${DEPMAP_MODEL_MD5:-}
fi

OUT=""
DEPMAP=1
DEPMAP_LATEST=1
FORCE=0
while [ $# -gt 0 ]; do
  case "$1" in
    -o) OUT=${2:-}; [ -n "$OUT" ] || { echo "ERROR: -o requires a value" >&2; exit 2; }; shift 2;;
    --no-depmap) DEPMAP=0; shift;;
    --depmap-24q4) DEPMAP_LATEST=0; shift;;
    --force) FORCE=1; shift;;
    -h|--help) sed -n '2,45p' "$0"; exit 0;;
    *) echo "ERROR: unknown option: $1 (use -o OUTDIR, --no-depmap, --depmap-24q4, --force)" >&2; exit 2;;
  esac
done
if [ -z "$OUT" ]; then
  if [ -d /opt/nanomatizer/userdb ]; then OUT=/opt/nanomatizer/userdb; else OUT=./nanocastle_db; fi
fi
mkdir -p "$OUT"
[ -w "$OUT" ] || { echo "ERROR: $OUT is not writable" >&2; exit 1; }

for tool in curl gzip md5sum; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: '$tool' is required" >&2; exit 1; }
done

log() { echo "[$(date '+%F %T')] $*"; }

# fetch URL DEST [MD5] [gz]
fetch() {
  local url=$1 dest=$2 md5=${3:-} gz=${4:-} tmp
  if [ -s "$dest" ] && [ "$FORCE" -eq 0 ]; then
    log "present, skipped: $dest (use --force to re-download)"
    return 0
  fi
  tmp="$dest.part"
  log "downloading $(basename "$dest") <- $url"
  curl -fL --retry 3 --retry-delay 5 -o "$tmp" "$url" || { rm -f "$tmp"; echo "ERROR: download failed: $url" >&2; return 1; }
  if [ "$gz" = "gz" ]; then
    gzip -dc "$tmp" > "$tmp.txt" && mv "$tmp.txt" "$tmp"
  fi
  [ -s "$tmp" ] || { rm -f "$tmp"; echo "ERROR: empty download: $url" >&2; return 1; }
  if [ -n "$md5" ]; then
    local got
    got=$(md5sum "$tmp" | awk '{print $1}')
    [ "$got" = "$md5" ] || { rm -f "$tmp"; echo "ERROR: md5 mismatch for $(basename "$dest") (expected $md5, got $got)" >&2; return 1; }
  fi
  mv "$tmp" "$dest"
  log "ok: $dest ($(wc -c < "$dest") bytes)"
}

fetch "$ONCOKB_URL" "$OUT/OncoKB.tsv"
head -1 "$OUT/OncoKB.tsv" | grep -q "^Hugo Symbol" \
  || { echo "ERROR: $OUT/OncoKB.tsv does not look like the OncoKB cancer gene list (header changed?)" >&2; exit 1; }
[ "$(head -1 "$OUT/OncoKB.tsv" | awk -F'\t' '{print $7}')" = "Gene Type" ] \
  || { echo "ERROR: $OUT/OncoKB.tsv: column 7 is not 'Gene Type' (format changed)" >&2; exit 1; }

fetch "$DECIPHER_URL" "$OUT/population_cnv_grch38.txt" "" gz
head -1 "$OUT/population_cnv_grch38.txt" | grep -q "^#population_cnv_id" \
  || { echo "ERROR: $OUT/population_cnv_grch38.txt does not look like the DECIPHER population CNV file" >&2; exit 1; }

fetch "$TRRUST_URL" "$OUT/trrust_rawdata.human.tsv" "$TRRUST_MD5"
[ "$(head -1 "$OUT/trrust_rawdata.human.tsv" | awk -F'\t' '{print NF}')" = "4" ] \
  || { echo "ERROR: $OUT/trrust_rawdata.human.tsv does not look like the TRRUST table (4 tab-separated columns expected)" >&2; exit 1; }

fetch "$DGV_URL" "$OUT/DGV.GS.hg38.gff3" "$DGV_MD5"
head -1 "$OUT/DGV.GS.hg38.gff3" | grep -q "copy_number_variation_region" \
  || { echo "ERROR: $OUT/DGV.GS.hg38.gff3 does not look like the DGV Gold Standard gff3" >&2; exit 1; }

# Latest "DepMap Public" release from the portal file list
# (csv: release,release_date,filename,url,md5_hash). Sets DEPMAP_* and returns
# 0 only if both files of the same release are found.
depmap_latest() {
  local list=$OUT/.depmap_files.csv rel line
  curl -fsSL --retry 2 --max-time 120 -o "$list" "$DEPMAP_LIST_URL" 2>/dev/null || { rm -f "$list"; return 1; }
  head -1 "$list" | grep -q "^release,release_date,filename,url" || { rm -f "$list"; return 1; }
  rel=$(awk -F',' 'NR>1 && $1 ~ /^"?DepMap Public [0-9][0-9]Q[1-4]"?$/ && $3 ~ /^"?CRISPRGeneEffect\.csv"?$/ {gsub(/"/,"",$1); print $1}' "$list" \
        | sort -t' ' -k3,3 | tail -1)
  [ -n "$rel" ] || { rm -f "$list"; return 1; }
  line=$(awk -F',' -v r="$rel" '{gsub(/"/,"")} $1==r && $3=="CRISPRGeneEffect.csv" {print $4 "\t" $5; exit}' "$list")
  DEPMAP_CRISPR_URL=${line%%$'\t'*}; DEPMAP_CRISPR_MD5=${line##*$'\t'}
  line=$(awk -F',' -v r="$rel" '{gsub(/"/,"")} $1==r && $3=="Model.csv" {print $4 "\t" $5; exit}' "$list")
  DEPMAP_MODEL_URL=${line%%$'\t'*}; DEPMAP_MODEL_MD5=${line##*$'\t'}
  rm -f "$list"
  case "$DEPMAP_CRISPR_URL" in *://*) ;; *) return 1;; esac
  case "$DEPMAP_MODEL_URL"  in *://*) ;; *) return 1;; esac
  # md5 column is optional in the list: keep it only if it looks like one
  echo "$DEPMAP_CRISPR_MD5" | grep -qE '^[0-9a-f]{32}$' || DEPMAP_CRISPR_MD5=""
  echo "$DEPMAP_MODEL_MD5"  | grep -qE '^[0-9a-f]{32}$' || DEPMAP_MODEL_MD5=""
  DEPMAP_RELEASE=$rel
  return 0
}
depmap_fallback() {
  DEPMAP_RELEASE=$DEPMAP_FALLBACK_RELEASE
  DEPMAP_CRISPR_URL=$DEPMAP_FALLBACK_CRISPR_URL; DEPMAP_CRISPR_MD5=$DEPMAP_FALLBACK_CRISPR_MD5
  DEPMAP_MODEL_URL=$DEPMAP_FALLBACK_MODEL_URL;   DEPMAP_MODEL_MD5=$DEPMAP_FALLBACK_MODEL_MD5
}
depmap_check() {
  head -c 8 "$OUT/Model.csv" | grep -q "^ModelID," \
    || { echo "ERROR: $OUT/Model.csv does not look like the DepMap Model.csv" >&2; return 1; }
  head -c 200 "$OUT/CRISPRGeneEffect.csv" | grep -qE '^[^,]*,[A-Za-z0-9.-]+ \([0-9]+\),' \
    || { echo "ERROR: $OUT/CRISPRGeneEffect.csv does not look like the DepMap gene effect matrix" >&2; return 1; }
}

if [ "$DEPMAP" -eq 1 ]; then
  if [ -s "$OUT/Model.csv" ] && [ -s "$OUT/CRISPRGeneEffect.csv" ] && [ "$FORCE" -eq 0 ]; then
    log "present, skipped: DepMap ($(cat "$OUT/DEPMAP_RELEASE.txt" 2>/dev/null || echo "release not recorded"); use --force to re-download)"
  else
    if [ "$DEPMAP_USER_URLS" -eq 0 ]; then
      if [ "$DEPMAP_LATEST" -eq 1 ] && depmap_latest; then
        log "DepMap: latest release on the portal is '$DEPMAP_RELEASE'"
      else
        [ "$DEPMAP_LATEST" -eq 1 ] && log "DepMap: the portal file list is not reachable from a script (human verification); using $DEPMAP_FALLBACK_RELEASE. For a newer release see the header of this script."
        depmap_fallback
      fi
    fi
    FORCE=1
    if ! { fetch "$DEPMAP_MODEL_URL" "$OUT/Model.csv" "$DEPMAP_MODEL_MD5" \
           && fetch "$DEPMAP_CRISPR_URL" "$OUT/CRISPRGeneEffect.csv" "$DEPMAP_CRISPR_MD5" \
           && depmap_check; }; then
      [ "$DEPMAP_USER_URLS" -eq 0 ] && [ "$DEPMAP_RELEASE" != "$DEPMAP_FALLBACK_RELEASE" ] \
        || { echo "ERROR: DepMap download failed ($DEPMAP_RELEASE)" >&2; exit 1; }
      log "DepMap: download of '$DEPMAP_RELEASE' failed, falling back to $DEPMAP_FALLBACK_RELEASE"
      rm -f "$OUT/Model.csv" "$OUT/CRISPRGeneEffect.csv"
      depmap_fallback
      fetch "$DEPMAP_MODEL_URL"  "$OUT/Model.csv"            "$DEPMAP_MODEL_MD5"
      fetch "$DEPMAP_CRISPR_URL" "$OUT/CRISPRGeneEffect.csv" "$DEPMAP_CRISPR_MD5"
      depmap_check || exit 1
    fi
    echo "$DEPMAP_RELEASE" > "$OUT/DEPMAP_RELEASE.txt"
    log "DepMap release in $OUT: $DEPMAP_RELEASE"
  fi
else
  log "DepMap skipped (--no-depmap): the -T option will not be available"
fi

log "DONE. Mount the folder in the container with:  -v $(cd "$OUT" && pwd):/opt/nanomatizer/userdb"
