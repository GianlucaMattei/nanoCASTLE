#!/usr/bin/env Rscript
# =============================================================================
# build_chm13_dbs.r — one-off conversion of the coordinate-based annotation
# databases from GRCh38 to T2T-CHM13v2.0 (UCSC liftOver, official T2T chain
# grch38-chm13v2). Files already present in OUTDIR are NOT rebuilt (run it as
# many times as you like: only the missing ones are converted).
#
# Same rule as GASOLINE_2.0_T2T/DB/build_DB_chm13.r: two-point liftover (start
# and end lifted separately, unique mappings only); a record is kept only if
# both ends map to the same chromosome and the size stays compatible
# (difference < 1 kb, or new/old size ratio 0.5-2). chrM is copied unchanged
# (rCRS in both references).
#
# Converted (output keeps the SAME file name and format as the GRCh38 file, so
# the annotation step reads them unchanged from the chm13 folder):
#   DGV.GS.hg38.gff3                      (user database, cols 1,4,5)
#   population_cnv_grch38.txt             (DECIPHER, user database, chr/start/end)
#   GRCh38.nr_{deletions,duplications,insertions}*.tsv   (ClinVar/dbVar, outermost_start/stop)
#   gnomad.v4.1.cnv.all.vcf               (POS + INFO END)
#   gnomad.v4.1.sv.sites.selected.columns.bed  (start/end; POS2/END2 for inter-chromosomal)
#   ENCODE.blacklist.CNV.bed
#   EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3
#   EnsDb.Hsapiens.v86 genes + CDS + exons -> EnsDb.Hsapiens.v86_chm13.rds (list(genes, cds, exons))
#   (a gene is kept only if the gene body and ALL its CDS/exon pieces are converted)
# Not converted (no coordinates): OncoKB, TRRUST, STRING, ClinGen, DepMap, MSigDB.
# A report with the converted fraction of every file is written to OUTDIR/LIFTOVER_REPORT.tsv.
#
# Usage:
#   Rscript build_chm13_dbs.r <DBDIR_GRCh38> <OUTDIR_chm13> [CHAIN] [LIFTOVER_BIN]
#   env ONLY=file1,file2  convert only these files (e.g. the user databases
#                          DGV.GS.hg38.gff3,population_cnv_grch38.txt at run time)
# =============================================================================
suppressPackageStartupMessages(library(data.table))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) < 2) stop("Usage: build_chm13_dbs.r <DBDIR_GRCh38> <OUTDIR_chm13> [CHAIN] [LIFTOVER_BIN]")
DB    <- a[1]; OUT <- a[2]
CHAIN <- if (length(a) >= 3) a[3] else Sys.getenv("CHM13_CHAIN", "/opt/nanomatizer/reference/grch38-chm13v2.chain")
LO    <- if (length(a) >= 4) a[4] else Sys.getenv("LIFTOVER", "liftOver")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(CHAIN))
REPORT <- file.path(OUT, "LIFTOVER_REPORT.tsv")
report <- function(f, n_in, n_out) {
  line <- data.table(file = f, records_GRCh38 = n_in, records_chm13 = n_out,
                     converted_pct = round(100 * n_out / max(n_in, 1), 2), date = format(Sys.Date()))
  fwrite(line, REPORT, sep = "\t", append = file.exists(REPORT))
  cat(sprintf("%-60s %10d -> %10d (%.2f%%)\n", f, n_in, n_out, 100 * n_out / max(n_in, 1)))
}
ONLY <- strsplit(Sys.getenv("ONLY", ""), ",")[[1]]
todo <- function(f) {
  if (length(ONLY) && !f %in% ONLY) return(FALSE)
  src <- file.path(DB, f); dst <- file.path(OUT, f)
  if (file.exists(dst)) { cat("present, skipped:", f, "\n"); return(FALSE) }
  if (!file.exists(src)) { cat("source missing, skipped:", f, "\n"); return(FALSE) }
  TRUE
}

## ---- liftOver of single 1-based positions; NA if not uniquely mapped ----
lift_points <- function(chr, pos) {
  i <- tempfile(fileext = ".bed"); o <- tempfile(fileext = ".bed")
  ok <- !is.na(chr) & !is.na(pos)
  fwrite(data.table(chr[ok], as.integer(pos[ok]) - 1L, as.integer(pos[ok]), which(ok)), i, sep = "\t", col.names = FALSE)
  system2(LO, c(i, CHAIN, o, "/dev/null"), stdout = FALSE, stderr = FALSE)
  r <- data.table(chr = rep(NA_character_, length(pos)), pos = rep(NA_real_, length(pos)))
  if (file.exists(o) && file.size(o) > 0) {
    l <- fread(o, header = FALSE, col.names = c("chr", "s", "e", "idx"))
    r[l$idx, `:=`(chr = l$chr, pos = as.numeric(l$e))]
  }
  unlink(c(i, o)); r
}
## two-point liftover of intervals (chr with 'chr' prefix): list(chr, start, end, keep)
lift_iv <- function(chr, start, end) {
  start <- as.numeric(start); end <- as.numeric(end)
  isM <- chr %in% c("chrM", "chrMT")
  s <- lift_points(chr, start); e <- lift_points(chr, end)
  ns <- pmin(s$pos, e$pos); ne <- pmax(s$pos, e$pos)
  old <- end - start; new <- ne - ns
  keep <- !is.na(s$chr) & !is.na(e$chr) & s$chr == e$chr &
          (abs(new - old) < 1000 | (new >= 0.5 * old & new <= 2 * old))
  keep[is.na(keep)] <- FALSE
  nc <- s$chr; nc[isM] <- "chrM"; ns[isM] <- start[isM]; ne[isM] <- end[isM]; keep[isM] <- TRUE
  list(chr = nc, start = ns, end = ne, keep = keep)
}
addchr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
rmchr  <- function(x) sub("^chr", "", x)

## ---- 1. DGV gold standard gff3 (no header; chr-prefixed; cols 1,4,5) ----
f <- "DGV.GS.hg38.gff3"
if (todo(f)) {
  x <- fread(file.path(DB, f), header = FALSE, sep = "\t", quote = "")
  l <- lift_iv(addchr(x$V1), x$V4, x$V5)
  x[, `:=`(V1 = l$chr, V4 = as.integer(l$start), V5 = as.integer(l$end))]
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", col.names = FALSE, quote = FALSE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 2. DECIPHER population CNV (header '#population_cnv_id'; chr without prefix) ----
f <- "population_cnv_grch38.txt"
if (todo(f)) {
  x <- fread(file.path(DB, f), sep = "\t")
  l <- lift_iv(addchr(x$chr), x$start, x$end)
  x[, `:=`(chr = rmchr(l$chr), start = as.integer(l$start), end = as.integer(l$end))]
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", quote = FALSE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 3. ClinVar/dbVar non-redundant SVs (line 1 '#NR_SVs GRCh38', header '#chr'; chr without prefix) ----
for (f in c("GRCh38.nr_deletions.pathogenic.tsv", "GRCh38.nr_duplications.somatic.tsv", "GRCh38.nr_insertions.somatic.tsv")) {
  if (!todo(f)) next
  first <- readLines(file.path(DB, f), n = 1)
  x <- fread(file.path(DB, f), sep = "\t", skip = 1, colClasses = "character")
  l <- lift_iv(addchr(x[["#chr"]]), x$outermost_start, x$outermost_stop)
  x[, `:=`(`#chr` = rmchr(l$chr), outermost_start = as.character(as.integer(l$start)),
           outermost_stop = as.character(as.integer(l$end)))]
  con <- file(file.path(OUT, f), "w"); writeLines(sub("GRCh38", "T2T-CHM13v2.0 (liftOver from GRCh38)", first), con); close(con)
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", quote = FALSE, append = TRUE, col.names = TRUE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 4. gnomAD v4.1 CNV vcf (POS + INFO END) ----
f <- "gnomad.v4.1.cnv.all.vcf"
if (todo(f)) {
  hdr <- grep("^#", readLines(file.path(DB, f), n = 5000), value = TRUE)
  x <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(DB, f))), header = FALSE, sep = "\t", quote = "")
  end <- as.numeric(sub(";.*", "", sub(".*?(^|;)END=", "", x$V8, perl = TRUE)))
  l <- lift_iv(addchr(x$V1), x$V2, end)
  x[, `:=`(V1 = l$chr, V2 = as.integer(l$start))]
  x[, V8 := mapply(function(info, e) sub("(^|;)END=[0-9]+", paste0("\\1END=", e), info), V8, as.integer(l$end))]
  hdr <- c(hdr[-length(hdr)], "##liftover=T2T-CHM13v2.0 from GRCh38 (UCSC liftOver, chain grch38-chm13v2; POS and END only)", hdr[length(hdr)])
  writeLines(hdr, file.path(OUT, f))
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 5. gnomAD v4.1 SV sites (header; chrom/start/end; CHR2 POS2/END2) ----
f <- "gnomad.v4.1.sv.sites.selected.columns.bed"
if (todo(f)) {
  x <- fread(file.path(DB, f), sep = "\t", colClasses = list(character = c("CHR2", "END2", "POS2")))
  l <- lift_iv(addchr(x$chrom), x$start, x$end)
  x[, `:=`(chrom = l$chr, start = as.integer(l$start), end = as.integer(l$end))]
  for (cc in c("POS2", "END2")) {
    v <- suppressWarnings(as.numeric(x[[cc]])); k <- !is.na(v) & !is.na(x$CHR2) & x$CHR2 != "NA"
    if (any(k)) { p <- lift_points(addchr(x$CHR2[k]), v[k])
      nv <- x[[cc]]; nv[k] <- ifelse(is.na(p$pos), "NA", as.character(as.integer(p$pos))); set(x, j = cc, value = nv) }
  }
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", quote = FALSE, na = "NA")
  report(f, length(l$keep), sum(l$keep))
}

## ---- 6. ENCODE blacklist (no header; chr start end type) ----
f <- "ENCODE.blacklist.CNV.bed"
if (todo(f)) {
  x <- fread(file.path(DB, f), header = FALSE, sep = "\t")
  l <- lift_iv(x$V1, x$V2 + 1, x$V3)          # bed 0-based start
  x[, `:=`(V1 = l$chr, V2 = as.integer(l$start) - 1L, V3 = as.integer(l$end))]
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", col.names = FALSE, quote = FALSE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 7. Ensembl Regulatory Build v112 gff3 (chr without prefix; cols 4,5) ----
f <- "EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3"
if (todo(f)) {
  dir.create(file.path(OUT, "EnsemblRegulatory"), showWarnings = FALSE)
  hdr <- grep("^#", readLines(file.path(DB, f), n = 200), value = TRUE)   # '##gff-version 3' is checked by run_module.sh
  x <- fread(cmd = paste("grep -v '^#'", shQuote(file.path(DB, f))), header = FALSE, sep = "\t", quote = "")
  l <- lift_iv(addchr(x$V1), x$V4, x$V5)
  x[, `:=`(V1 = rmchr(l$chr), V4 = as.integer(l$start), V5 = as.integer(l$end))]
  writeLines(c(hdr, "# liftover: T2T-CHM13v2.0 from GRCh38 (UCSC liftOver, chain grch38-chm13v2)"), file.path(OUT, f))
  fwrite(x[l$keep], file.path(OUT, f), sep = "\t", col.names = FALSE, quote = FALSE, append = TRUE)
  report(f, length(l$keep), sum(l$keep))
}

## ---- 8. EnsDb.Hsapiens.v86 genes + CDS (used for gene/exon annotation of every SV/CNV) ----
f <- "EnsDb.Hsapiens.v86_chm13.rds"
if (!file.exists(file.path(OUT, f)) && (!length(ONLY) || f %in% ONLY)) {
  suppressPackageStartupMessages({library(EnsDb.Hsapiens.v86); library(GenomicRanges); library(GenomicFeatures)})
  g <- genes(EnsDb.Hsapiens.v86)
  cds <- cdsBy(EnsDb.Hsapiens.v86, by = "gene")
  exs <- exonsBy(EnsDb.Hsapiens.v86, by = "gene")
  lg <- lift_iv(addchr(as.character(seqnames(g))), start(g), end(g))
  cu <- unlist(cds, use.names = FALSE); cid <- rep(names(cds), lengths(cds))
  lc <- lift_iv(addchr(as.character(seqnames(cu))), start(cu), end(cu))
  eu <- unlist(exs, use.names = FALSE); eid <- rep(names(exs), lengths(exs))
  le <- lift_iv(addchr(as.character(seqnames(eu))), start(eu), end(eu))
  gchr <- function(id) lg$chr[match(id, names(g))]
  bad <- unique(c(cid[!lc$keep], cid[lc$keep & lc$chr != gchr(cid)], eid[!le$keep], eid[le$keep & le$chr != gchr(eid)],
                  names(g)[!lg$keep]))
  kg <- lg$keep & !(names(g) %in% bad)
  g2 <- GRanges(rmchr(lg$chr[kg]), IRanges(lg$start[kg], lg$end[kg]), strand = strand(g)[kg])
  mcols(g2) <- mcols(g)[kg, ]; names(g2) <- names(g)[kg]
  kc <- lc$keep & !(cid %in% bad)
  c2 <- GRanges(rmchr(lc$chr[kc]), IRanges(lc$start[kc], lc$end[kc]), strand = strand(cu)[kc])
  mcols(c2) <- mcols(cu)[kc, , drop = FALSE]
  cds2 <- split(c2, factor(cid[kc], levels = unique(cid[kc])))
  ke <- le$keep & !(eid %in% bad)
  e2 <- GRanges(rmchr(le$chr[ke]), IRanges(le$start[ke], le$end[ke]), strand = strand(eu)[ke])
  mcols(e2) <- mcols(eu)[ke, , drop = FALSE]
  exs2 <- split(e2, factor(eid[ke], levels = unique(eid[ke])))
  saveRDS(list(genes = g2, cds = cds2, exons = exs2, note = "EnsDb.Hsapiens.v86 lifted to T2T-CHM13v2.0 (seqnames without chr, as EnsDb)"),
          file.path(OUT, f))
  report("EnsDb.Hsapiens.v86_genes", length(g), length(g2))
  report("EnsDb.Hsapiens.v86_CDS_genes", length(cds), length(cds2))
  report("EnsDb.Hsapiens.v86_exon_genes", length(exs), length(exs2))
} else if (file.exists(file.path(OUT, f))) cat("present, skipped:", f, "\n")
cat("DONE ->", OUT, "\n")
