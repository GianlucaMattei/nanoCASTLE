# ============================================================================
# nanomatizer_dbs_v2.r  —  caricamento database di annotazione (STEP 6)
# ----------------------------------------------------------------------------
# Funzione unica: load_annotation_dbs(dir)
#
# Il blocco di caricamento VERBATIM (dal driver originale step6_annotation.r
# L2738-2941) vive in:  nanomatizer_dbs_block_v2.r  (stessa cartella).
# Unica modifica rispetto all'originale: i path assoluti "Databsese/..." sono
# diventati "paste0(dir, '/...')"; nessun'altra riga del blocco e' toccata.
#
# Scopo del wrapper: nel driver originale il blocco girava in TOP LEVEL, quindi
# gli oggetti DB finivano nel globalenv (dove le funzioni core li cercano per
# scope: annotateODB, annotateTF, annotate.clingen, annotate.depmap.*, ...).
# Qui il blocco viene valutato in un ambiente temporaneo (con `dir` definito)
# e TUTTI gli oggetti creati vengono copiati nel globalenv: comportamento
# identico, ma la cartella DB resta un parametro.
#
# Uso:  load_annotation_dbs(dir = "/home/.../Databsese", depmap = "none")
#
# `depmap` (spec effettiva: none | "tumor" | lista file | "default" back-door)
# e' il valore che arriva dalla pipeline: il blocco DB carica
# CRISPRGeneEffect.csv/Model.csv (oggetti depmap, depmap.models, depmap.cns,
# ...) SOLO quando depmap == "default" (back-door tecnica dei 3 gruppi
# CNS/Pediatric/MBL, non piu' esposta dalla CLI). Per "none" (nuovo default:
# DepMap skipped) e "custom" quei file non servono affatto, quindi non devono
# esistere (es. immagine Docker senza i due CSV). Per "tumor" i due CSV sono
# letti DIRETTAMENTE da build.depmap.groups (nanomatizer_depmap_v2.r), non da
# qui — il gruppo tumore e' costruito al volo dal Model.csv + CRISPRGeneEffect.csv
# dell'utente (release arbitraria).
# ============================================================================

load_annotation_dbs <- function(dir = "/opt/nanomatizer/Nanomatizer/Databsese",
                                depmap = "none") {

    ## --- localizza il blocco (stessa cartella di questo file) ---
    args.all <- commandArgs(trailingOnly = FALSE)
    file.arg <- grep("^--file=", args.all, value = TRUE)
    here <- if (length(file.arg) == 1) {
        normalizePath(dirname(sub("^--file=", "", file.arg[[1]])), winslash = "/")
    } else {
        "/opt/nanomatizer/Nanomatizer/Scripts/NewScripts/step6/module_v2"
    }
    blockfile <- file.path(here, "nanomatizer_dbs_block_v2.r")

    if (!dir.exists(dir)) stop("Cartella DB mancante: ", dir, " (STEP6) — STOP")
    dir <- normalizePath(dir, winslash = "/")
    if (!file.exists(blockfile)) stop("Blocco DB mancante: ", blockfile, " — STOP")

    cat("[STEP6] DB: caricamento da ", dir, "\n", sep = "")
    t0 <- Sys.time()

    ## --- valuta il blocco in ambiente temporaneo, parent = globalenv ---
    e <- new.env(parent = globalenv())
    assign("dir", dir, envir = e)
    assign("depmap.mode", depmap, envir = e)
    eval(parse(blockfile), envir = e)

    ## --- copia tutti gli oggetti DB nel globalenv (come nel top-level originale) ---
    n <- 0
    for (nm in ls(e, all.names = TRUE)) {
        if (nm == "dir" || nm == "depmap.mode") next
        assign(nm, get(nm, envir = e), envir = globalenv())
        n <- n + 1
    }

    cat("[STEP6] DB: caricamento completato in ",
        round(as.numeric(difftime(Sys.time(), t0, units = "secs"))),
        "s — ", n, " oggetti DB in globalenv\n", sep = "")
    invisible(NULL)
}
