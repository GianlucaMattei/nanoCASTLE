#!/bin/bash
# =============================================================================
# build_ng_target.sh — target NanoGLADIATOR a finestre 1 Mb da un genoma montato
#
# Compagno di run_pipeline_sequential.sh (strategia fingerprint del genoma):
# quando il genoma montato NON e' lo standard GRCh38 no-alt della pipeline
# (fingerprint 24 primari diverso dalla costante STANDARD_FP), la pipeline
# genera qui le finestre 1 Mb per quel genoma e lo "inietta" nel tool con un
# SYMLINK (vedi la pipeline): NESSUN file dei tool e' scritto o modificato,
# il symlink e' un artefatto di runtime ricreato a ogni run.
#
# Uso:
#   build_ng_target.sh --fasta F --fai F.fai --name NAME --out DIR [--window 1000000]
#
# Layout prodotto in DIR/NAME/ (IDENTICO a un target bundled nel tool,
# es. NanoGLADIATOR_1.0/data/targets/hg38/hg38w1000000 — gli stessi file
# letti da BAM2RC.sh / LoadRC.R / LibraryGenomeRC.R):
#   NAME_chromosome.txt    i 24 cromosomi primari, una riga, spazi
#   TargetWindow.bed       finestre 1 Mb consecutive (chr start end aN)
#   NAME.RData             MyTarget  (matrice 4 colonne, colname "a")
#   GCC/GCC.<chr>.RData    GCContent (numeric, FRAZIONE 0-1, ordine finestre)
#   MAP/Map.<chr>.RData    MapMed    (numeric, costante 1.0 — see below)
#
# Regole e trade-off dichiarati:
#   - 24 primari CH&R (chr1..chr22, chrX, chrY): il tool li usa hardcoded
#     (LoadRC.R / NanoGLADIATORInferenceGenome.R). Assenti nel .fai -> STOP.
#   - le finestre partono da 10001 (stesso start del target stock hg38w100000)
#     e sono consecutive fino alla fine del cromosoma (ultima troncata).
#     Differenza dichiarata col target stock 1 Mb (hg38w1000000), che usa
#     maschera dei gap UCSC (2851 finestre); qui: 1 Mb consecutive senza
#     maschera (il gap-masking dipende dal genoma e non e' generalizzabile).
#   - GC: REALE per finestra, bedtools nuc sul FASTA montato. La scala del
#     tool e' la FRAZIONE (lo stock hg38w100000 ha max ~0.61;
#     NanoGLADIATORNormalizationGenome.R moltiplica per 100). Se bedtools
#     restituisce percentuali (max > 1.5) normalizziamo dividendo per 100.
#   - MAP: costante 1.0 per finestra. Nel percorso standard (GRCh38) si usa
#     il target stock con la mappability reale; questo builder serve le
#     varianti GRCh38 e, non avendo il .wig di mappability dell'utente,
#     assume mappability ~1 (riferimento no-alt). CorrectMAP/QuantileMAP a
#     MAP costante 100 (dopo *100) lasciano i RC invariati (un solo bucket).
#   - FRB/ NON generata: non e' letta in analisi whole-genome (verificato
#     nel tool: solo LibraryFractionCall.R la usa, non nel flusso prepare).
#
# Idempotenza: se il target e' gia' completo, esce 0 senza toccare nulla.
# Stato parziale (file a meta') = STOP: non sovrascrivo, decisione richiesta.
# Richiede (dall'ambiente chiamante): R, bedtools, awk, md5sum. La pipeline
# lo lancia dentro l'env conda nanogladiator (tutti presenti).
# =============================================================================
set -euo pipefail

FASTA=; FAI=; NAME=; OUT=; WINDOW=1000000
die() { echo "STOP: $*" >&2; exit 1; }
while [ $# -gt 0 ]; do
  case "$1" in
    --fasta)  FASTA=${2:-};  shift 2 || die "--fasta richiede un valore";;
    --fai)    FAI=${2:-};    shift 2 || die "--fai richiede un valore";;
    --name)   NAME=${2:-};   shift 2 || die "--name richiede un valore";;
    --out)    OUT=${2:-};    shift 2 || die "--out richiede un valore";;
    --window) WINDOW=${2:-}; shift 2 || die "--window richiede un valore";;
    *) die "opzione sconosciuta: $1";;
  esac
done
[ -n "$FASTA" ] && [ -n "$FAI" ] && [ -n "$NAME" ] && [ -n "$OUT" ] \
  || die "mancano argomenti (usa --fasta --fai --name --out [--window])"
case "$WINDOW" in ''|*[!0-9]*) die "--window non valido: $WINDOW";; esac
[ -f "$FASTA" ] || die "FASTA mancante: $FASTA"
[ -f "$FAI" ]   || die ".fai mancante: $FAI"
command -v R >/dev/null 2>&1 || die "R non disponibile (lancia il builder dentro l'env conda nanogladiator)"
command -v bedtools >/dev/null 2>&1 || die "bedtools non disponibile (lancia il builder dentro l'env conda nanogladiator)"

D=$OUT/$NAME

# ---------- idempotenza / anti-sovrascrittura (read-only, prima di scrivere) --
# NB: la dir di output e' un caso valido alla prima run (non esiste):
# "|| true" dentro la pipeline perche' con set -o pipefail il rc non zero
# di find su dir assente abortirebbe lo script (set -e).
N_GC=$( { find "$D/GCC" -maxdepth 1 -name "*.RData" 2>/dev/null || true; } | wc -l)
N_MAP=$( { find "$D/MAP" -maxdepth 1 -name "*.RData" 2>/dev/null || true; } | wc -l)
if [ -f "$D/TargetWindow.bed" ] && [ -f "$D/$NAME.RData" ] && [ -f "$D/${NAME}_chromosome.txt" ] \
   && [ "$N_GC" -eq 24 ] && [ "$N_MAP" -eq 24 ]; then
  echo "OK: target gia' presente e completo: $D"
  exit 0
fi
if [ -f "$D/TargetWindow.bed" ] || [ -f "$D/$NAME.RData" ] || [ "$N_GC" -gt 0 ] || [ "$N_MAP" -gt 0 ]; then
  die "target in stato parziale in $D (bed=$([ -f "$D/TargetWindow.bed" ] && echo si || echo no), GCC=$N_GC/24, MAP=$N_MAP/24) — non tocco nulla, decisione richiesta"
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT

# ---------- 24 primari dal .fai (PRIMA di creare l'output dir) ---------------
awk '$1 ~ /^(chr([1-9]|1[0-9]|2[0-2])|chrX|chrY)$/ {print $1"\t"$2}' "$FAI" > "$WORK/prim"
N_PRIM=$(wc -l < "$WORK/prim")
[ "$N_PRIM" -eq 24 ] || die "fai: attesi 24 cromosomi primari (chr1..chr22, chrX, chrY), trovate $N_PRIM — il tool lavora solo sui 24 primari con prefisso chr"

mkdir -p "$D/GCC" "$D/MAP"

# ---------- finestre: bed + MyTarget + chromosome.txt (R) --------------------
# invocazione identica a quella del tool (R --slave --args <csv> < script)
R --slave --args "$WORK/prim,$WINDOW,$NAME,$D" <<'REOF' >/dev/null
# stesso idiom degli script R del tool (vedi MakeReadCount.R)
vars <- commandArgs()[length(commandArgs())]
sv <- unlist(strsplit(vars, ","))
PrimIn <- sv[1]; Window <- as.numeric(sv[2]); Name <- sv[3]; DirOut <- sv[4]
prim <- read.table(PrimIn, header = FALSE)
chr.order <- c(paste0("chr", 1:22), "chrX", "chrY")
stopifnot(setequal(as.character(prim[[1]]), chr.order))
prim <- prim[match(chr.order, prim[[1]]), , drop = FALSE]
chrs <- as.character(prim[[1]]); lens <- as.numeric(prim[[2]])
# finestre consecutive da 10001 (convenzione tool: salta i primi 10 kb)
rows <- list()
for (i in seq_along(chrs)) {
  s <- 10001
  while (s <= lens[i]) {
    e <- min(s + Window - 1, lens[i])
    rows[[length(rows) + 1]] <- c(chrs[i], as.character(s), as.character(e))
    s <- e + 1
  }
}
n <- length(rows)
stopifnot(n > 0)
TotalMatrixCoord <- do.call(rbind, rows)   # n x 3, senza colnames (come CreateBedGenome.R)
a <- paste("a", seq_len(n), sep = "")
MyTarget <- cbind(TotalMatrixCoord, a)     # 4a colonna = id, colname "a" (verbatim del tool)
MyChr <- unique(TotalMatrixCoord[, 1])
write.table(MyTarget, file.path(DirOut, "TargetWindow.bed"),
            col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
save(MyTarget, file = file.path(DirOut, paste0(Name, ".RData")))
write.table(data.frame(rbind(MyChr)), file.path(DirOut, paste0(Name, "_chromosome.txt")),
            col.names = FALSE, row.names = FALSE, quote = FALSE)
cat("OK: finestre generate: ", n, "finestre 1 Mb su 24 cromosomi\n", sep = "")
REOF

# ---------- GCC reale (bedtools nuc) + MAP costante 1.0 (R) ------------------
R --slave --args "$WORK/prim,$FASTA,$NAME,$D" <<'REOF' >/dev/null
# stesso idiom degli script R del tool (vedi MakeReadCount.R)
vars <- commandArgs()[length(commandArgs())]
sv <- unlist(strsplit(vars, ","))
PrimIn <- sv[1]; Fasta <- sv[2]; Name <- sv[3]; DirOut <- sv[4]
prim <- read.table(PrimIn, header = FALSE)
chr.order <- c(paste0("chr", 1:22), "chrX", "chrY")
prim <- prim[match(chr.order, prim[[1]]), , drop = FALSE]
chrs <- as.character(prim[[1]])
win <- read.table(file.path(DirOut, "TargetWindow.bed"), header = FALSE)
for (ch in chrs) {
  idx <- which(win[[1]] == ch)
  nwin <- length(idx)
  if (nwin == 0) stop("nessuna finestra per ", ch)
  # bed per-chr in 0-based (start-1, come fa TargetCreate.sh del tool)
  tb <- file.path(tempdir(), paste0("bed.", ch, ".bed"))
  write.table(cbind(win[idx, 1], as.integer(win[idx, 2]) - 1,
                    as.integer(win[idx, 3]), ch),
              file = tb, col.names = FALSE, row.names = FALSE, sep = "\t", quote = FALSE)
  # bedtools nuc: prendo la colonna pct_gc PARSENDOLA DALL'HEADER (l'ordine
  # delle colonne cambia tra versioni: in v2.31.1 e' $5 per un bed a 3
  # colonne; il tool stock usava cut -f6 col suo bedtools vecchio, dove quella
  # era la GC). Robusto a qualsiasi versione.
  out.nuc <- system(paste("bedtools nuc -fi", shQuote(Fasta), "-bed", shQuote(tb)),
                    intern = TRUE)
  if (length(out.nuc) < 2 || !grepl("^#", out.nuc[1]))
    stop("bedtools nuc: output inatteso per ", ch)
  # NB: l'header bedtools numera le colonne (es. "6_pct_gc" per un bed a
  # 3 colonne, "5_pct_gc" per un bed a 4): il match e' sul suffisso pct_gc,
  # robusto a versione e al numero di colonne del bed di input.
  hcols <- unlist(strsplit(sub("^#", "", out.nuc[1]), "\t"))
  idx.gc <- which(grepl("^[_0-9]*pct_gc$", hcols))
  length(idx.gc) == 1 || stop("bedtools nuc: colonna pct_gc non trovata (",
                              length(idx.gc), " match in [",
                              paste(hcols, collapse=","), "]) per ", ch)
  cols.gc <- vapply(out.nuc[-1], function(line)
                     strsplit(line, "\t")[[1]][idx.gc], character(1))
  GCContent <- suppressWarnings(as.numeric(cols.gc))
  if (length(GCContent) != nwin) stop("bedtools nuc: attese ", nwin,
                                      " righe per ", ch, ", ottenute ", length(GCContent))
  if (anyNA(GCContent)) stop("GC non numerico per ", ch)
  # scala del tool = FRAZIONE 0-1 (stock max ~0.61; NormalizationGenome.R
  # moltiplica per 100). pct_gc recente gia' 0-1; se una versione vecchia
  # restituisse percentuali 0-100, normalizziamo qui.
  if (max(GCContent, na.rm = TRUE) > 1.5) GCContent <- GCContent / 100
  save(GCContent, file = file.path(DirOut, "GCC", paste0("GCC.", ch, ".RData")))
  # MAP costante 1.0 (mappability ~1 per riferimento no-alt, trade-off dichiarato)
  MapMed <- rep(1.0, nwin)
  save(MapMed, file = file.path(DirOut, "MAP", paste0("Map.", ch, ".RData")))
  cat("OK:", ch, nwin, "finestre (GC mediano", format(median(GCContent)), ")\n")
}
REOF

# ---------- verifica finale --------------------------------------------------
N_GC=$(find "$D/GCC" -maxdepth 1 -name "*.RData" 2>/dev/null | wc -l)
N_MAP=$(find "$D/MAP" -maxdepth 1 -name "*.RData" 2>/dev/null | wc -l)
[ "$N_GC" -eq 24 ] && [ "$N_MAP" -eq 24 ] \
  || die "verifica finale fallita: GCC=$N_GC/24 MAP=$N_MAP/24"
N_WIN=$(wc -l < "$D/TargetWindow.bed")
echo "OK: target 1 Mb completo in $D ($N_WIN finestre, 24 cromosomi)"
