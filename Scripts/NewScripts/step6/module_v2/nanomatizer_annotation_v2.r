# ============================================================================
# nanomatizer_annotation_v2.r — STEP 6 annotazioni + aggregazione (driver v2)
# ----------------------------------------------------------------------------
# Scopo: riorganizzazione DELIZIATA del driver step6_annotation.r,
# comportamento IDENTICO (stesse funzioni, stessi argomenti, stessa catena).
#
# Architettura:
#   1. source delle 4 componenti (logiche VERBATIM, non riscritte):
#        - nanomatizer.annotation.local.r   (annotazioni per-tool, step6/)
#        - nanomatizer_core_functions.r     (46 funzioni core, step6/)
#        - nanomatizer_driver_helpers.r     (exonsEnsDb, remove.ins.in.introns,
#                                            build_sv_list, create.helping.col)
#        - nanomatizer_dbs_v2.r             (load_annotation_dbs, parametrica)
#   2. tre funzioni di orchestrazione pulite:
#        - run_tool_annotations(...)   : le 5 chiamate per-tool (argomenti uguali)
#        - run_aggregation(...)        : carica i 5 CSV per-tool + catena completa
#                                        centrality -> ... -> aggregato xlsx
#        - main(...)                   : parsing argomenti + timing
#
# Logica delle 46 funzioni core e del blocco DB: INVARIATA (source).
#
# Uso (tramite run_module.sh o Rscript diretto):
#   nanomatizer_annotation_v2.r [sample=X] [outdir=X] [dbdir=X]
#     sample = nome campione (obbligatorio, nessun default;
#              i path per-tool sono tutti derivati da 'sample', quindi il
#              driver e' generico per qualunque campione hg38)
#                               [mode=all|annotation|aggregation]
# Defaults:
#   sample=(obbligatorio)
#   outdir=/opt/nanomatizer/Nanomatizer/Results
#   dbdir =/opt/nanomatizer/Nanomatizer/Databsese
#   mode  =all
# ============================================================================

## ------------------------- percorso dello script --------------------------
args.all   <- commandArgs(trailingOnly = FALSE)
file.arg   <- grep("^--file=", args.all, value = TRUE)
SCRIPT_DIR <- if (length(file.arg) == 1) {
    normalizePath(dirname(sub("^--file=", "", file.arg[[1]])), winslash = "/")
} else {
    "/opt/nanomatizer/Nanomatizer/Scripts/NewScripts/step6/module_v2"   # fallback (normale: container)
}
P6         <- file.path(SCRIPT_DIR, "..")          # cartella step6/ (componenti condivise)

## ----------------------------- librerie -----------------------------------
suppressPackageStartupMessages({
    library(gtools)
    library(GenomicRanges)
    library(data.table)
    library(IRanges)
    library(igraph)
})

## --------------------------- componenti (verbatim) ------------------------
# ORDINE DI source = ORDINE TOP-LEVEL dell'originale step6_annotation.r:
#   1. lib locale  -> crea exonsEnsDb LOCALE (= CDS) usato dalle annotazioni per-tool
#   2. core        -> le 46 funzioni (solo definizioni, niente top-level)
#   3. run_tool_annotations()   [come gli loop per-tool dell'originale, AVANTI
#                                rispetto al block driver-helper]
#   4. load_annotation_dbs()    [come il blocco DB top-level dell'originale]
#   5. driver-helper (source in run_aggregation) -> sovrascrive exonsEnsDb con
#                                exonsBy(..., by='gene') PRIMA di build_sv_list
# NB: nanomatizer_driver_helpers.r NON va sourced all'avvio: la sua riga top-level
# `exonsEnsDb <- exonsBy(..., by='gene')` deve avvenire DOPO le annotazioni per-tool
# (che leggono l'exonsEnsDb-CDS della lib locale), esattamente come nell'originale.
source(file.path(P6, "nanomatizer.annotation.local.r"))
source(file.path(P6, "nanomatizer_core_functions.r"))
source(file.path(SCRIPT_DIR, "nanomatizer_dbs_v2.r"))
# Opzione -depmap (generalizzazione catena DepMap). Definizioni autonome
# (solo data.table/base): sicure da source all'avvio.
source(file.path(SCRIPT_DIR, "nanomatizer_depmap_v2.r"))

## ============================================================================
## OVERWRITES DI ROBUSTEZZA (solo per slot a 0 righe, subset di tool)
## ----------------------------------------------------------------------------
## Le 46 funzioni core sono logiche VERBATIM del tool (non si modifica il
## sorgente). Due di esse usano l'idioma `1:length(...)` / `by = 1:nrow(...)`,
## che in R con length/nrow == 0 degenera (1:0 = c(1,0)) e fallisce con
## "subscript out of bounds". Questo si verifica SOLO quando uno slot e' a 0
## righe, cioe' nel caso di AGGREGAZIONE SU SUBSET DI TOOL (tool esclusi ->
## placeholder 0-righe). Con tutti i tool pieni (default) il comportamento e'
## IDENTICO (seq_along(x) == 1:length(x) e seq_len(n) == 1:n per n>0), quindi
## l'anti-rotture byte-identico resta garantito. Ridefiniamo qui (nel driver,
## non nel sorgente) solo le 4 funzioni collegate, lasciando intatto il tool.
## ============================================================================
annotate.centrality.v2 <- function(input.list, db){
    for(i in c(1,3,4,5)){
        universe <- unique(input.list[[i]]$symbol)
        ind.list <- sapply(universe, function(x){
            grep(paste0("^", x, "$"), gsub(" \\(X\\)", "", input.list[[i]]$symbol))
        })
        ind.string <- sapply(universe, function(x){
            match(x, db$symbol)
        })
        ind.rm <- unlist(lapply(ind.string, function(z) {!is.na(z)}))
        ind.list <- ind.list[ind.rm]
        ind.string <- ind.string[ind.rm]
        input.list[[i]][, `:=`(counts = NA_integer_, tot.score = NA_integer_)]
        for(j in seq_along(ind.list)){           # was: 1:length(ind.list)
            values.to.add <- db[ind.string[j], ]
            input.list[[i]][ind.list[[j]], `:=`(counts = values.to.add$counts, tot.score = values.to.add$tot.score)]
        }
    }
    input.list
}
annotate.depmap.cns <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max")
        if (nrow(input.list[[i]]) == 0) next      # slot vuoto: nessuna righe da annotare
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.cns[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = seq_len(nrow(input.list[[i]])), .SDcols = "symbol"]  # was: 1:nrow(...)
    }
    input.list
}
annotate.depmap.cns.pediatric <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max")
        if (nrow(input.list[[i]]) == 0) next
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.cns.pediatric[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = seq_len(nrow(input.list[[i]])), .SDcols = "symbol"]
    }
    input.list
}
annotate.depmap.cns.mbl <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max")
        if (nrow(input.list[[i]]) == 0) next
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.mbl[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = seq_len(nrow(input.list[[i]])), .SDcols = "symbol"]
    }
    input.list
}
## --- auto.annotate: 0-row-safe ---------------------------------------------
## auto.annotate(x, ...) e' la funzione che OGNI ramo di annotate.results6.list.vGene
## (e di annotate.results6.list) chiama per ogni (db, sv.type). Con input a 0 righe
## (placeholder di un tool ESCLUSO nel subset) la selezione per sv.type e' sempre
## vuota, quindi la funzione ritorna NULL per OGNI sv.type; cosi'
## rbindlist(list(11 x NULL), fill=TRUE) produce una tabella 0-righe 0-COLONNE e
## l'aritmetica a valle del ramo NG100K (`length := end - start`, L1038) fallisce
## con "non-numeric argument to binary operator". I rami Gasoline/CuteSV/Sniffles
## sopravvivono (non hanno aritmetica post-rbindlist) ma restano a 0 colonne.
## Soluzione: se nrow(x)==0 restituiamo copy(x) (0 righe, COLONNE invariate, con
## start/end), cosi' il rbindlist conserva le colonne e ogni operazione a valle
## diventa "0 righe" (sicura). Con nrow(x)>0 delega VERBATIM all'implementazione
## del core (auto.annotate.og) -> identico comportamento originario.
## NB: nel path default (tutti i tool pieni) NESSUN input a auto.annotate e' a 0
## righe, quindi la branch nrow==0 non viene MAI eseguita -> anti-rotture byte-identico.
auto.annotate.og <- auto.annotate
auto.annotate <- function(x, sv.label, sv.type, db, db.tag, seqnames.col, start.col, end.col){
    if (nrow(x) == 0) return(copy(x))   # placeholder vuoto: conserva le colonne, 0 righe
    auto.annotate.og(x, sv.label, sv.type, db, db.tag, seqnames.col, start.col, end.col)
}

## ----------------- PLACEHOLDER 0-righe per tool esclusi (subset) -----------
## Quando si aggrega un SUBSET di tool, gli slot dei tool ESCLUSI devono restare
## nella lista (posizione fissa 1..5) ma a 0 righe, cosi' la catena posizionale
## (build_sv_list, annotate.*, collapse) funziona invariata. Se il CSV per-tool
## e' ASSENTE, costruiamo una data.table vuota collo SCHEMA ESATTO del CSV reale
## (header catturato dai CSV validati di nuovirisultati). I nomi file e gli
## header seguono il naming dell'originale.
SLOT.TOOLNAME <- c(gasoline = "gasoline", intertrans = "intertrans", cutesv = "cutesv",
                   sniffles2 = "sniffles2", ng100k = "ng100k")
TOOL.CSV.PATH <- function(sample, gladiator.window_tag) {
    list(
        gasoline   = paste0("Gasoline.",                sample, ".germline.SVs.exons.csv"),
        intertrans = paste0("GasolineInterTrans.",      sample, ".germline.SVs.exons.csv"),
        cutesv     = paste0("CuteSV.",                  sample, ".germline.SVs.exons.csv"),
        sniffles2  = paste0("Sniffles.",                sample, ".germline.SVs.exons.csv"),
        ng100k     = paste0("NanoGladiator.", gladiator.window_tag, ".", sample, ".germline.SVs.exons.csv")
    )
}
# Header reali (riga 1) dei CSV per-tool validati. Usati SOLO per i placeholder
# 0-righe quando il CSV e' assente: la struttura (colonne, anche duplicate)
# deve coincidere con quella che fread produrrebbe dal file reale.
TOOL.CSV.HEADER <- list(
    gasoline   = "EnsID,seqnames,start,end,Size,Start.SD,End.SD,Coesion,Clustered.Signatures,Gap.Del.Signatures,Gap.Ins.Signatures,Split.Del.Signatures,Split.Ins.Signatures,Split.Dup.Signatures,Split.Inv.Signatures,Del.Prob,Ins.Prob,Dup.Prob,Inv.Prob,SV.Type,Cluster.Category,RONorm,Total.Reads,Alt.Reads,Geno0.Prob,Geno1.Prob,Geno2.Prob,strand,gene_id,symbol,gene_name,gene_biotype,SYMBOL",
    intertrans = "seqnames.1,seqnames.2,start.1,start.2,Reads,Alt.Reads,Prob.Homo.WT,Prob.Ete,Prom.Homo,Prob.Gen,gene1,regulatory1,regulatory.type1,gene2,regulatory2,regulatory.type2,nearest.gene1,nearest.gene1.distance,nearest.gene2,nearest.gene2.distance",
    cutesv     = "EnsID,seqnames,start,end,sv,reads,length,seqnames.1,start.1,end.1,width,strand,gene_id,symbol,gene_name,gene_biotype,SYMBOL",
    sniffles2  = "EnsID,seqnames,start,end,length,svtype,filter,info,format,sample,seqnames.1,start.1,end.1,width,strand,gene_id,symbol,gene_name,gene_biotype,SYMBOL",
    ng100k     = "seqnames,start,end,boh1,CNV.Type,Copies,Freq,Freq.Error,Shift,start,end,symbol,gene_biotype"
)
# Carica il CSV per-tool; se il file e' assente e il tool e' ESCLUSO, restituisce
# un placeholder 0-righe con lo schema reale. Se assente e tool SELEZIONATO -> STOP.
load.tool.csv <- function(tool, sample, output.dir, gladiator.window_tag, selected) {
    fn   <- TOOL.CSV.PATH(sample, gladiator.window_tag)[[tool]]
    path <- file.path(output.dir, fn)
    if (file.exists(path)) {
        args <- if (tool == "intertrans") list(sep = ",", fill = TRUE) else list()
        do.call(fread, c(list(path), args))
    } else if (selected) {
        stop("per-tool CSV missing for SELECTED tool (", tool, "): ", path,
             "\n  (to aggregate a tool its CSV is required; to exclude it just do not select it)")
    } else {
        log.step(paste0("AGGREG: ", tool, " escluso e CSV assente -> placeholder 0-righe (", fn, ")"))
        fread(text = paste0(TOOL.CSV.HEADER[[tool]], "\n"))
    }
}

## ------------------------------- helper -----------------------------------
log.step <- function(msg) {
    cat("[" , format(Sys.time(), "%F %T"), "] ", msg, "\n", sep = "")
    flush.console()
}

# ----------------------------------------------------------------------------
# run_tool_annotations
#   Esegue le 5 annotazioni per-tool, con esattamente gli stessi argomenti
#   del driver originale (step6_annotation.r L11-46):
#     CuteSV:     annotateCuteSV(sample.file, f, output.dir)
#     Gasoline:   annotateGasoline2(sample.file, f, output.dir)
#     InterTrans: annotateInterTrans2(sample.file, f, output.dir)
#     Sniffles:   annotateSniffles(sample.file, f, output.dir)
#     Gladiator:  annotateGladiator(sample.file, f, "W1Mb", output.dir,
#                                   shift.thr = 0.35, freq.error.thr = 0.35)
#   I CSV per-tool vengono scritti in output.dir (root), come nell'originale.
#   Parametro tools (default "all"): esegue SOLO le annotazioni dei tool
#   selezionati (abilita "un solo caller + annotazione" e l'annotazione su
#   sottoinsieme). "all" = le 5 chiamate, ordine e argomenti invariati.
# ----------------------------------------------------------------------------
run_tool_annotations <- function(sample, output.dir, tools = "all",
                                 gladiator.window.dir = c("R6_W1Mb"),
                                 gladiator.window_tag = "W1Mb",
                                 gladiator.shift.thr = 0.35,
                                 gladiator.freq.error.thr = 0.35) {

    keep.slots <- parse.tools(tools)
    on <- function(slot) slot %in% keep.slots
    log.step(paste0("ANNOT: tools=", tools, " (slot ", paste(keep.slots, collapse = ","), ")"))

    ## ------------------- CuteSV (slot 3) -------------------
    if (on(TOOLS.SLOT["cutesv"])) {
        log.step("ANNOT: CuteSV")
        def <- paste0(output.dir, "/CuteSV/")
        sample.file <- list.files(paste0(def, sample, "/"), pattern = paste0(sample, ".*vcf"), full.names = TRUE)
        annotateCuteSV(sample.file, sample, output.dir)
    }

    ## ------------------- Gasoline (slot 1) -------------------
    if (on(TOOLS.SLOT["gasoline"])) {
        def <- paste0(output.dir, "/Gasoline/Parse/")
        sample.file <- list.files(paste0(def, sample, "/"), pattern = paste0(sample, ".Results.tsv"), full.names = TRUE)
        annotateGasoline2(sample.file, sample, output.dir)
    }

    ## ------------------- Gasoline InterTrans (slot 2) -------------------
    if (on(TOOLS.SLOT["intertrans"])) {
        def <- paste0(output.dir, "/Gasoline/Parse/")
        sample.file <- list.files(paste0(def, sample, "/"), pattern = paste0(sample, ".Results.Trans.tsv"), full.names = TRUE)
        annotateInterTrans2(sample.file, sample, output.dir)
    }

    ## ------------------- Sniffles (slot 4) -------------------
    if (on(TOOLS.SLOT["sniffles"])) {
        def <- paste0(output.dir, "/Sniffles2/")
        sample.file <- list.files(paste0(def, sample, "/"), pattern = "vcf", full.names = TRUE)
        if (length(sample.file) > 0) {
            log.step("ANNOT: Sniffles")
            annotateSniffles(sample.file, sample, output.dir)
        } else {
            stop("Sniffles: no VCF file in ", def, sample, "/ — STOP")
        }
    }

    ## ------------------- NanoGladiator (run6, W1Mb) (slot 5) -------------------
    # RUN6: finestra R6_W1Mb (finestre 1 Mb, target hg38w1000000; risultati
    # con lo script R del server); il tag di output e' "W1Mb" (CSV
    # NanoGladiator.W1Mb.<sample>.germline.SVs.exons.csv).
    # ATTENZIONE: il pattern e' letteralmente "^FractionCallResults_" (come l'originale):
    # nel subdir Results/<LABEL>/ matcha SOLO il file consolidato
    # "FractionCallResults_<LABEL>.txt" (i file chrN_FractionCallResults_* non matchano).
    if (on(TOOLS.SLOT["gladiator"])) {
        def <- paste0(output.dir, "/NanoGLADIATOR/")
        for (w in gladiator.window.dir) {
            ops <- list.files(paste0(def, sample, "/", w, "/Results/"))
            sample.file <- list.files(paste0(def, sample, "/", w, "/Results/", ops),
                                      pattern = "^FractionCallResults_", full.names = TRUE)
            log.step(paste0("ANNOT: NanoGladiator (", w, ", ", length(sample.file), " file FractionCallResults)"))
            if (length(sample.file) == 0)
                stop("NanoGladiator: no FractionCallResults_ file in ", def, sample, "/", w, "/ — STOP")
            annotateGladiator(sample.file, sample, gladiator.window_tag, output.dir,
                              shift.thr = gladiator.shift.thr, freq.error.thr = gladiator.freq.error.thr)
        }
    }

    invisible(NULL)
}

# ------------------------- mapping -tools -> slot --------------------------
# Dopo annotate.results6.list.vGene la lista è nominata e POSITIONALE:
#   1=GASOLINE  2=INTERTRANS  3=CUTESV  4=SNIFFLES2  5=NG100K(GladiatoR)
# NB: lo slot 2 (INTERTRANS) NON alimenta collapse.genes (che usa solo 1,3,4,5),
# quindi includerlo/escluderlo non cambia l'output (comunque onorato).
TOOLS.SLOT <- c(gasoline = 1, intertrans = 2, cutesv = 3,
                sniffles = 4, sniffles2 = 4, gladiator = 5, ng100k = 5)
parse.tools <- function(spec) {
    spec <- tolower(trimws(spec))
    if (spec == "all") return(1:5)
    toks <- trimws(strsplit(spec, ",")[[1]]); toks <- tolower(toks[nzchar(toks)])
    bad  <- setdiff(toks, names(TOOLS.SLOT))
    if (length(bad) > 0)
        stop("tools: invalid value(s): ", paste(bad, collapse = ", "),
             " (expected all|sniffles|cutesv|gasoline|intertrans|gladiator)")
    slots <- unique(TOOLS.SLOT[toks])
    if (length(slots) == 0) stop("tools: no tool selected")
    slots
}

# ----------------------------------------------------------------------------
# run_aggregation
#   Carica i 5 CSV per-tool prodotti da run_tool_annotations e riproduce
#   ESATTAMENTE la catena del driver originale (L2945-3157), con due opzioni:
#     -tools  : "all" (default, identico) oppure elenco di tool da INCLUDERE
#               (sniffles|cutesv|gasoline|intertrans|gladiator). I tool esclusi
#               vengono azzerati (0 righe) DOPO l'annotazione e PRIMA del
#               collapse: collapse.genes scarta gli slot vuoti (anti-rotture).
#     -depmap : "none" (NUOVO DEFAULT della pipeline: nessun gruppo DepMap —
#               il default CNS e' stato rimosso), oppure "FILE[,FILE..]"
#               (0..N file custom, formato CRISPRGeneEffect). "default" resta
#               come back-door tecnica (3 gruppi CNS/Pediatric/MBL verbatim,
#               NON esposta dalla CLI) per l'anti-rotture.
#     -tumor  : <OncotreeLineage> di Model.csv (es. "Lung") oppure "MBL"
#               (ematologici); con suffisso "+Pediatric" (es.
#               "CNS/Brain+Pediatric") seleziona il solo sottoinsieme
#               pediatrico (AgeCategory 'Pediatric') di quel lineage.
#               Genera UN gruppo DepMap dal Model.csv + CRISPRGeneEffect.csv
#               montati (release dell'utente, file opzionali). Mutuamente
#               esclusivo coi file custom.
#               Per none/custom/tumor si usano le funzioni .v2 dinamiche
#               (collapse/helping/final.cleaning).
#   CATENA (default, tutti gli argomenti invariati):
#     build_sv_list -> annotate.centrality.v2 -> rename.duplicated.col ->
#     annotate.clingen -> [depmap] -> annotate.results6.list.vGene ->
#     annotate_encode_blacklist -> [tools-zeroing] -> collapse.genes(.v2) ->
#     clean -> replace.long.cells -> rm.duplicated.info -> flatten ->
#     sort.chrom -> create.helping.col(.v2) ->
#     annotate.macrocategories(annotateODB(annotateTF)) -> final.cleaning(.v2)
#     -> filtro Shift -> riordino colonne -> xlsx foglio S.<sample>
# ----------------------------------------------------------------------------
run_aggregation <- function(sample, output.dir,
                            gladiator.window_tag = "W1Mb",
                            encode.bed = "/opt/nanomatizer/Nanomatizer/Databsese/ENCODE.blacklist.CNV.bed",
                            shift.thr = 0.3,
                            tools = "all",
                            depmap = "none",
                            tumor = "",
                            dbdir = NULL) {

    ## --- parsing opzioni ---
    keep.slots <- parse.tools(tools)
    dm <- build.depmap.groups(depmap, dbdir = dbdir, tumor = tumor)
    log.step(paste0("AGGREG: tools=", tools, " (slot ", paste(keep.slots, collapse = ","),
                    ") | depmap=", dm$mode))
    if (length(dm$skipped) > 0) {
        log.step(paste0("AGGREG: -depmap option: file skipped (invalid): ",
                        paste(dm$skipped, collapse = ", ")))
    }

    ## --- caricamento 5 CSV per-tool (nomi identici all'originale) ---
    ## Tool SELEZIONATI: il CSV e' obbligatorio (assente -> STOP, come prima).
    ## Tool ESCLUSI: se il CSV e' presente lo carichiamo (per poi azzerarlo,
    ## come nel comportamento precedente); se assente costruiamo un placeholder
    ## 0-righe con lo schema reale. Questo abilita l'aggregazione su subset di
    ## tool anche quando i CSV dei tool esclusi non esistono affatto.
    log.step("AGGREG: caricamento CSV per-tool")
    selected.set <- unique(SLOT.TOOLNAME[keep.slots])
    env <- new.env()
    for (slot in 1:5) {
        tool <- names(SLOT.TOOLNAME)[slot]
        selected <- slot %in% keep.slots
        dt <- load.tool.csv(tool, sample, output.dir, gladiator.window_tag, selected)
        env[[if (tool == "intertrans")
                 paste0("gasoline.", sample, ".intertrans")
             else paste0(tool, ".", sample, ".exons")]] <- dt
    }

    ## --- costruzione lista SV (funzione verbatim) ---
    result_list <- build_sv_list(sample, env = env)
    result_list <- annotate.centrality.v2(result_list, string.abs.high)
    result_list <- rename.duplicated.col(result_list)
    result_list <- annotate.clingen(result_list)

    ## --- DEPMAP: default (verbatim, back-door) oppure none/custom/tumor
    ##     (funzione .v2 dinamica; tumor = un solo gruppo dall'opzione -T) ---
    if (dm$mode == "default") {
        result_list <- annotate.depmap.cns(result_list)
        result_list <- annotate.depmap.cns.pediatric(result_list)
        result_list <- annotate.depmap.cns.mbl(result_list)
    } else if (dm$mode == "custom" || dm$mode == "tumor") {
        result_list <- annotate.depmap.custom(result_list, dm$groups)
    }
    # mode == "none": nessuna colonna DepMap aggiunta (default della pipeline).

    result_list <- annotate.results6.list.vGene(result_list)
    result_list <- lapply(result_list, function(x) annotate_encode_blacklist(x, encode.bed))

    ## --- TOOLS: azzeramento (0 righe) degli slot esclusi, PRIMA del collapse ---
    ## NB: il collapse (verbatim e .v2) usa SOLO gli slot 1/3/4/5 (gasoline,
    ## cutesv, sniffles, ng100k). Lo slot 2 (intertrans) non viene mai usato,
    ## e dopo annotate_encode_blacklist puo' essere una STRINGA (il caso
    ## "skipping one annotation" ritorna invisibile, senza data.table): per
    ## questo azzeriamo [FALSE, ] SOLO gli slot che sono davvero data.table.
    excluded <- setdiff(1:5, keep.slots)
    if (length(excluded) > 0) {
        for (p in excluded) {
            if (is.data.table(result_list[[p]])) {
                result_list[[p]] <- result_list[[p]][FALSE, ]
            }
        }
        log.step(paste0("AGGREG: tool esclusi (slot azzerati): ",
                        paste(excluded, collapse = ",")))
    }

    overlap.names <- c("decipher.reciprocal.overlap", "gAD.CNV.reciprocal.overlap",
                       "gAD.reciprocal.overlap", "clinvar.reciprocal.overlap")
    ## Scegliamo il collapse:
    ##  - verbatim  collapse.genes       -> SOLO quando tools="all" e depmap="default"
    ##    (garantisce byte-identical rispetto al vecchio output)
    ##  - robusto   collapse.genes.v2    -> quando e' escluso un tool oppure depmap!=default
    ##    (stesso schema; per depmap="default" riceve le 15 colonne standard, quindi
    ##     riproduce lo schema verbatim + mette il guard su 'tag'/'seqnames' assenti)
    DEPMAP.15.COL <- c(
        "DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max",
        "DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max",
        "DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max")
    use.robust <- (length(excluded) > 0) || (dm$mode != "default")
    if (!use.robust) {
        result_list <- collapse.genes(result_list, overlap.names)
    } else {
        dm.cols.for.collapse <- if (dm$mode == "default") DEPMAP.15.COL else dm$depmap.cols
        result_list <- collapse.genes.v2(result_list, overlap.names, dm.cols.for.collapse)
    }

    cols.to.process <- names(result_list)
    result_list[, (cols.to.process) := lapply(.SD, clean), .SDcols = cols.to.process]

    result_list <- tryCatch({
        replace.long.cells(result_list)
    }, error = function(e) {
        cat("ERROR on result_list: object not modified. (", conditionMessage(e), ")\n", sep = "")
        result_list
    })

    cols.to.process <- names(result_list)
    sheet <- sort.chrom(flatten.data.table(
        rm.duplicated.info(result_list[, (cols.to.process) := lapply(.SD, clean), .SDcols = cols.to.process])
    ))

    ## --- colonne di aiuto + annotazione DB finali ---
    ## Quando use.robust (tool esclusi o depmap!=default) usiamo le versioni .v2,
    ## che sono robuste a colonne mancanti; per depmap="default" ricevono le 15
    ## colonne DepMap standard e le 4 colonne sum-score standard (schema identico
    ## al verbatim). Altrimenti (all+default) usiamo il verbatim byte-identico.
    SUMSCORE.4.COL <- c("CNS.Sum.Score", "CNS.Sum.Score.match",
                        "Pediatric.Sum.Score", "Pediatric.Sum.Score.match")
    if (use.robust) {
        dm.cols.post    <- if (dm$mode == "default") DEPMAP.15.COL else dm$depmap.cols
        sumscore.cols.post <- if (dm$mode == "default") SUMSCORE.4.COL else dm$sumscore.cols
        sheet <- create.helping.col.v2(sheet, sumscore.cols.post)
    } else {
        sheet <- create.helping.col(sheet)
    }
    sheet <- annotate.macrocategories(annotateODB(annotateTF(sheet, tf)))
    if (use.robust) {
        sheet <- final.cleaning.v2(sheet, dm.cols.post, sumscore.cols.post)
    } else {
        sheet <- final.cleaning(sheet)
    }

    ## --- filtro Shift + riordino colonne + scrittura xlsx ---
    sheet.shift <- sheet[as.numeric(sheet$Shift) < -shift.thr | as.numeric(sheet$Shift) > shift.thr | is.na(as.numeric(sheet$Shift)), ]
    sheet.shift <- setcolorder(sheet.shift, c(2:5, 1, 6:ncol(sheet.shift)))

    out.xlsx <- paste0(output.dir, "/Nanomatizer_", sample, "_aggregated.xlsx")
    writexl::write_xlsx(setNames(list(sort.chrom(sheet.shift)), paste0("S.", sample)), out.xlsx)
    log.step(paste0("STEP6-AGGREGATO: file scritto ", out.xlsx,
                    " righe = ", nrow(sheet.shift), " colonne = ", ncol(sheet.shift)))
    return(out.xlsx)
}

# ----------------------------------------------------------------------------
# main
# ----------------------------------------------------------------------------
main <- function() {
    cfg <- list(
        sample = Sys.getenv("SAMPLE", ""),   # generico: nessuna sample preconfigurata
        outdir = "/opt/nanomatizer/Nanomatizer/Results",
        dbdir  = "/opt/nanomatizer/Nanomatizer/Databsese",
        mode   = "all",                      # all | annotation | aggregation
        tools  = "all",                      # all | elenco (sniffles|cutesv|gasoline|intertrans|gladiator)
        depmap = "none",                     # none (DEFAULT: DepMap skipped) | FILE[,FILE..] | "default" (back-door tecnica, non esposta dalla CLI)
        tumor  = ""                          # -T: OncotreeLineage di Model.csv | MBL | "" (skipped)
    )
    for (a in commandArgs(trailingOnly = TRUE)) {
        kv <- strsplit(a, "=", fixed = TRUE)[[1]]
        if (length(kv) == 2 && kv[1] %in% names(cfg)) cfg[[kv[1]]] <- kv[2]
    }
    if (!cfg$mode %in% c("all", "annotation", "aggregation"))
        stop("invalid mode: ", cfg$mode, " (expected all|annotation|aggregation)")
    # Spec effective per build.depmap.groups:
    #   tumor non vuoto -> "tumor" (UN gruppo dai CSV in dbdir, release utente)
    #   altrimenti      -> il valore -depmap (none | FILE... | default-back-door)
    cfg$depmap_spec <- if (trimws(cfg$tumor) != "") "tumor" else cfg$depmap

    t_start <- Sys.time()
    log.step(paste0("STEP6 v2 AVVIO — sample=", cfg$sample, " mode=", cfg$mode,
                    " outdir=", cfg$outdir))

    ## Ordine ESATTO del top-level originale:
    ##   (a) loop per-tool          -> exonsEnsDb = quello CDS della lib locale
    ##   (b) blocco DB              -> load_annotation_dbs (oggetti in globalenv)
    ##   (c) block driver-helper    -> source: sovrascrive exonsEnsDb (exons by gene)
    ##   (d) freads 5 CSV + catena  -> run_aggregation
    if (cfg$mode %in% c("all", "annotation")) {
        run_tool_annotations(cfg$sample, cfg$outdir, tools = cfg$tools)
    }

    if (cfg$mode %in% c("all", "aggregation")) {
        load_annotation_dbs(cfg$dbdir, depmap = cfg$depmap_spec)
        source(file.path(P6, "nanomatizer_driver_helpers.r"))
        run_aggregation(cfg$sample, cfg$outdir,
                        tools     = cfg$tools,
                        encode.bed = file.path(cfg$dbdir, "ENCODE.blacklist.CNV.bed"),
                        depmap = cfg$depmap_spec,
                        tumor  = cfg$tumor,
                        dbdir  = cfg$dbdir)
    }

    log.step(paste0("STEP6 v2 COMPLETATO in ",
                    round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1), " min"))
    invisible(NULL)
}

# main() viene eseguito SOLO quando lo script e' lanciato come programma
# (Rscript .../nanomatizer_annotation_v2.r). Se viene SOLO sourced (es. dal
# harness di test che imposta SOURCE_ONLY=TRUE), le funzioni restano disponibili
# ma la pipeline NON parte -> nessun sovrascrittura accidentale di output.
if (!isTRUE(get0("SOURCE_ONLY", ifnotfound = FALSE))) main()
