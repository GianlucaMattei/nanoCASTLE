# ============================================================================
# nanomatizer_depmap_v2.r — opzioni -depmap / -tumor (generalizzazione catena
# DepMap)
# ----------------------------------------------------------------------------
# Scopo: permettere di scegliere QUALI gruppi DEPMAP usare nell'aggregazione.
# L'opzione -depmap (spec) accetta:
#   - "none"         : nessun gruppo DepMap (NUOVO DEFAULT della pipeline —
#                      senza -T e senza file custom il DepMap viene SKIPPED;
#                      non esiste piu' un default CNS). Le colonne DepMap sono
#                      assenti, ma la catena NON si rompe (funzioni .v2
#                      dinamiche).
#   - "FILE[,FILE..]" : 0..N file custom, stesso formato di CRISPRGeneEffect.csv
#                      (colonna 1 = ModelID, colonne 2..N = geni). Ogni file
#                      genera un gruppo di colonne DepMap.<TAG>.{Min,Q1,Median,Q3,Max}
#                      dove <TAG> = nome del file (senza estensione).
#   - "default"      : back-door tecnica (NON piu' esposta dalla CLI della
#                      pipeline): i 3 gruppi standard CNS/Pediatric/MBL -> le 3
#                      funzioni verbatim annotate.depmap.cns/.pediatric/.mbl
#                      (comportamento IDENTICO alla baseline storica, usato
#                      solo per l'anti-rotture).
# L'opzione -tumor (tumor, argomento separato, mutuamente esclusivo coi file
# custom) accetta:
#   - <OncotreeLineage> esatto di Model.csv (es. "Lung",
#     "Ovary/Fallopian Tube", "CNS/Brain") oppure
#   - "<OncotreeLineage>+Pediatric" : sottoinsieme PEDIATRICO di quel lineage
#     (OncotreeLineage == X E AgeCategory == 'Pediatric'), es.
#     "CNS/Brain+Pediatric". E' l'equivalente CLI del gruppo "Pediatric"
#     della back-door default (tag = "CNS_Brain_Pediatric").
#   - "MBL" : tutti i modelli ematologici (DepmapModelType == "MBL").
# Genera UN gruppo con tag = valore -T (sanitizzato): le due fonti vengono
# lette da dbdir (Model.csv + CRISPRGeneEffect.csv, release dell'utente —
# i file restano OPZIONALI: servono solo quando -T e' specificato).
#
# Perche' servono le funzioni .v2:
#   collapse.genes / create.helping.col / final.cleaning (core verbatim)
#   hardcodano le 15 colonne DepMap (CNS/Pediatric/MBL). Per supportare gruppi
#   custom (o nessuno) servono versioni DINAMICHE che selezionano le colonne
#   DepMap in base ai gruppi reali presenti.
#
# INVARIA anti-rotture:
#   * Per il caso DEFAULT il driver non usa MAI le funzioni .v2: chiama le
#     funzioni verbatim -> output bit-identico alla baseline.
#   * Le funzioni .v2 replicano ESATTO il blocco non-DepMap delle versioni
#     originali (stessi nomi, stesso ordine, stesse rename) e sostituiscono
#     solo il blocco DepMap con le colonne dinamiche.
#
# Anti-rotture caricamento file custom: ogni file e' processato con tryCatch;
# un file invalido viene SCIPPATO (avviso a terminale) e gli altri procedono.
# ============================================================================

## ------------------- sanitizzazione tag da nome file ----------------------
depmap.sanitize.tag <- function(fname) {
    tag <- basename(fname)
    tag <- sub("\\.[^.]+$", "", tag)                 # togli estensione
    tag <- gsub("[^A-Za-z0-9_.]", "_", tag)          # solo alfanumerici/_/.
    tag <- gsub("^_+", "", sub("_+$", "", tag))
    if (nzchar(tag)) tag else "custom"
}

## ------------------- carico UN file custom -------------------------------
# Ritorna list(tag=, effects=<data.table geni x modelli>) oppure lancia error.
load.one.depmap.file <- function(f, dbdir = NULL) {
    cand <- f
    if (!is.null(dbdir)) cand <- c(cand, file.path(dbdir, f))
    cand <- cand[file.exists(cand)]
    if (length(cand) == 0) stop("file not found: ", f)
    f.real <- cand[1]

    dt <- data.table::fread(f.real, header = TRUE, check.names = FALSE)
    if (ncol(dt) < 2) stop("invalid format (>=2 columns required: ModelID + genes)")
    if (nrow(dt) < 1)  stop("empty file (no data rows)")

    # colonna 1 = ModelID (scartata); colonne 2..N = geni
    gene.raw  <- colnames(dt)[-1]
    gene.names <- gsub(" .*", "", gene.raw)                       # "A1BG (1)" -> "A1BG"
    gene.names <- gsub("[^A-Za-z0-9_.\\-]", "_", gene.names)
    effects <- as.data.table(dt[, -1, with = FALSE])
    colnames(effects) <- gene.names

    # rimuovi geni duplicati (tieni il primo)
    dup.idx <- duplicated(colnames(effects))
    if (any(dup.idx)) effects <- effects[, !dup.idx, with = FALSE]
    if (ncol(effects) < 1) stop("no valid gene column")

    list(tag = depmap.sanitize.tag(f.real), effects = effects, nmodels = nrow(effects))
}

## ------------------- build.depmap.groups ---------------------------------
# spec:  "none" | "default" (back-door) | "FILE1[,FILE2,...]"
# tumor: "" | <OncotreeLineage> | "MBL"   (mutuamente esclusivo coi file custom)
# Ritorna list(mode, groups, depmap.cols, sumscore.cols, skipped, files).
#   mode        : "default" | "none" | "custom" | "tumor"
#   groups      : named list tag -> data.table (per custom/tumor; default/none: vuoto)
#   depmap.cols : vettore (ordinato) di colonne DepMap.<tag>.{Min..Max}
#   sumscore.cols: vettore (ordinato) di colonne <TAG>.Sum.Score/.match
#   skipped     : file custom non validi (scippati)
#   files       : file custom richiesti
build.depmap.groups <- function(spec, dbdir = NULL, tumor = "") {
    spec  <- trimws(spec)
    tumor <- trimws(tumor)
    if (tolower(spec) == "default") {
        return(list(mode = "default", groups = list(),
                    depmap.cols = NULL, sumscore.cols = NULL,
                    skipped = character(0), files = character(0)))
    }
    if (tolower(spec) == "none") {
        if (tumor != "") {
            stop("tumor=", tumor, " but spec=none: -T is used only with spec 'tumor'")
        }
        return(list(mode = "none", groups = list(),
                    depmap.cols = character(0),
                    sumscore.cols = c("CNS.Sum.Score", "CNS.Sum.Score.match",
                                      "Pediatric.Sum.Score", "Pediatric.Sum.Score.match"),
                    skipped = character(0), files = character(0)))
    }
    # ---- tumor: UN gruppo da Model.csv + CRISPRGeneEffect.csv in dbdir ----
    # La spec per il mode tumor e' "tumor" (fissato da main); qui si costruisce
    # il gruppo leggendo i due CSV dell'utente (release arbitraria).
    # Mutuamente esclusivo coi file custom: -T con -d FILE... e' un errore.
    if (tumor != "") {
        if (tolower(spec) != "tumor")
            stop("-tumor '", tumor, "' and the custom file list are mutually exclusive (spec=", spec, ")")
        is.mbl    <- tolower(tumor) == "mbl"
        # Suffix "+Pediatric": sottoinsieme pediatrico di un lineage
        # (OncotreeLineage == X E AgeCategory == 'Pediatric'), es.
        # -T "CNS/Brain+Pediatric". E' l'equivalente CLI del gruppo
        # "Pediatric" della back-door default.
        is.ped    <- !is.mbl &&
                     grepl("\\+Pediatric$", tumor, ignore.case = TRUE)
        lineage   <- if (is.ped) sub("\\+Pediatric$", "", tumor,
                                     ignore.case = TRUE) else tumor
        # NB: NO basename() — i lineage contengono "/" (es. "Ovary/Fallopian Tube")
        # che basename tratterebbe come separatore di path. Sanitizzazione diretta.
        tag       <- tumor
        tag       <- gsub("[^A-Za-z0-9_.]", "_", tag)
        tag       <- gsub("^_+", "", sub("_+$", "", tag))
        if (!nzchar(tag)) tag <- "custom"
        models.f  <- if (!is.null(dbdir)) file.path(dbdir, "Model.csv") else "Model.csv"
        crisper.f <- if (!is.null(dbdir)) file.path(dbdir, "CRISPRGeneEffect.csv") else "CRISPRGeneEffect.csv"
        if (!file.exists(models.f))  stop("Model.csv not found: ", models.f)
        if (!file.exists(crisper.f)) stop("CRISPRGeneEffect.csv not found: ", crisper.f)
        models <- data.table::fread(models.f, header = TRUE, check.names = FALSE)
        for (col in c("ModelID", "OncotreeLineage", "DepmapModelType", "AgeCategory")) {
            if (!col %in% colnames(models))
                stop("Model.csv: column '", col, "' missing (unrecognized DepMap format in ", models.f, ")")
        }
        if (is.mbl) {
            sel <- subset(models, DepmapModelType == "MBL")
        } else if (is.ped) {
            # NB: '&' (non '&&'): subset.data.table valuta la condizione
            # in modo vettoriale.
            sel <- subset(models, OncotreeLineage == lineage &
                          AgeCategory == "Pediatric")
        } else {
            sel <- subset(models, OncotreeLineage == lineage)
        }
        if (nrow(sel) == 0) {
            desc <- if (is.mbl) "MBL"
                    else if (is.ped) paste0("OncotreeLineage '", lineage,
                                            "' pediatric")
                    else paste0("OncotreeLineage '", tumor, "'")
            stop(paste0("no DepMap models for ", desc,
                       " in Model.csv. Check the exact name (e.g. 'Lung', 'Ovary/Fallopian Tube', 'CNS/Brain+Pediatric') or use MBL for the hematologic models."))
        }
        # CRISPRGeneEffect.csv: colonna 1 = ModelID, colonne 2..N = geni.
        # Stessa pulizia dei nomi gene di load.one.depmap.file (replicata qui
        # perche' serve conservare la colonna ModelID per il subset).
        raw <- data.table::fread(crisper.f, header = TRUE, check.names = FALSE)
        if (ncol(raw) < 2) stop("CRISPRGeneEffect.csv invalid format (<= 1 column): ", crisper.f)
        modelids <- as.character(raw[[1]])                     # colonna 1 = ModelID (per posizione)
        raw <- raw[modelids %in% sel$ModelID, ]
        if (nrow(raw) == 0) stop("no model of tumor '", tumor, "' present in CRISPRGeneEffect.csv")
        gene.names <- gsub(" .*", "", colnames(raw)[-1])
        gene.names <- gsub("[^A-Za-z0-9_.\\-]", "_", gene.names)
        effects <- as.data.table(raw[, -1, with = FALSE])
        colnames(effects) <- gene.names
        dup.idx <- duplicated(colnames(effects))
        if (any(dup.idx)) effects <- effects[, !dup.idx, with = FALSE]
        if (ncol(effects) < 1) stop("no valid gene column in CRISPRGeneEffect.csv")
        groups <- list()
        groups[[tag]] <- effects
        cat("[DEPMAP] tumor group '", tag, "' (", nrow(effects), " models, ",
            ncol(effects), " genes) from ", crisper.f, " x ", models.f, "\n", sep = "")
        depmap.cols <- paste0("DepMap.", tag, c(".Min", ".Q1", ".Median", ".Q3", ".Max"))
        sumscore.cols <- c(paste0(tag, ".Sum.Score"), paste0(tag, ".Sum.Score.match"))
        return(list(mode = "tumor", groups = groups, depmap.cols = depmap.cols,
                    sumscore.cols = sumscore.cols,
                    skipped = character(0), files = c(models.f, crisper.f)))
    }
    # ---- custom: elenco di file separati da virgola ----
    files  <- trimws(strsplit(spec, ",")[[1]])
    files  <- files[nzchar(files)]
    groups <- list()
    skipped <- character(0)
    for (f in files) {
        res <- tryCatch(load.one.depmap.file(f, dbdir),
                        error = function(e) {
                            cat("[DEPMAP] ERROR loading '", f, "': ",
                                conditionMessage(e), "\n", sep = "")
                            NULL
                        })
        if (is.null(res)) { skipped <- c(skipped, f); next }
        groups[[res$tag]] <- res$effects
        cat("[DEPMAP] loaded group '", res$tag, "' (", res$nmodels,
            " models, ", ncol(res$effects), " genes) from ", f, "\n", sep = "")
    }
    # colonne DepMap ordinate in ordine di inserimento dei gruppi
    depmap.cols <- unlist(lapply(names(groups), function(t)
        paste0("DepMap.", t, c(".Min", ".Q1", ".Median", ".Q3", ".Max"))))
    # sum-score: sempre i 4 standard (CNS/Pediatric) per stabilita' di schema;
    # + una coppia per ogni tag custom (esclusi CNS/Pediatric/MBL, gia' standard)
    sumscore.cols <- c("CNS.Sum.Score", "CNS.Sum.Score.match",
                       "Pediatric.Sum.Score", "Pediatric.Sum.Score.match")
    custom.tags <- setdiff(names(groups), c("CNS", "Pediatric", "MBL"))
    if (length(custom.tags) > 0) {
        sumscore.cols <- c(sumscore.cols, unlist(lapply(custom.tags, function(t)
            c(paste0(t, ".Sum.Score"), paste0(t, ".Sum.Score.match")))))
    }
    list(mode = "custom", groups = groups, depmap.cols = depmap.cols,
         sumscore.cols = sumscore.cols, skipped = skipped, files = files)
}

## ------------------- annotate.depmap.custom ------------------------------
# Stesso meccanismo delle 3 verbatim (quantile per riga per slot c(1,3,4,5))
# ma con il gruppo/tag parametrizzato. Espressione quantile IDENTICA al verbatim.
annotate.depmap.custom <- function(input.list, groups) {
    tags <- names(groups)
    if (length(tags) == 0) return(input.list)
    for (i in c(1, 3, 4, 5)) {
        if (is.null(input.list[[i]]) || nrow(input.list[[i]]) == 0) next
        for (t in tags) {
            new_cols <- paste0("DepMap.", t, c(".Min", ".Q1", ".Median", ".Q3", ".Max"))
            eff <- groups[[t]]
            input.list[[i]][, (new_cols) := {
                current_symbol <- .SD$symbol
                v <- if (is.null(current_symbol) || length(current_symbol) == 0) {
                    NULL
                } else {
                    eff[[current_symbol]]
                }
                if (is.null(v) || length(v) == 0 || all(is.na(v))) {
                    list(NA_real_, NA_real_, NA_real_, NA_real_, NA_real_)
                } else {
                    as.list(quantile(v, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
                }
            }, by = 1:nrow(input.list[[i]]), .SDcols = "symbol"]
        }
    }
    input.list
}

## ------------------- collapse.genes.v2 -----------------------------------
# Identico a collapse.genes (core) tranne che def.col e' costruito
# dinamicamente: il blocco delle 15 colonne DepMap standard e' sostituito
# da `depmap.cols` (le colonne DepMap reali presenti).
collapse.genes.v2 <- function(x, overlap.names, depmap.cols) {
    def.col.orig <- c(
        "counts", "tot.score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT",
        "DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max",
        "DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max",
        "DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max",
        "decipher.#population_cnv_id", "decipher.deletion_observations", "decipher.deletion_frequency",
        "decipher.duplication_observations", "decipher.duplication_frequency", "decipher.observations",
        "decipher.frequency", "decipher.type", "decipher.sample_size", "decipher.study", "decipher.reciprocal.overlap",
        "gAD.CNV.sv.CNV", "gAD.CNV.symbol.CNV", "gAD.CNV.SF", "gAD.CNV.SF_nfe", "gAD.CNV.reciprocal.overlap",
        "gnomAD.svtype", "gnomAD.AF", "gnomAD.FREQ_HOMREF", "gnomAD.FREQ_HET", "gnomAD.FREQ_HOMALT",
        "gnomAD.AF_nfe", "gnomAD.FREQ_HOMREF_nfe", "gnomAD.FREQ_HET_nfe", "gAD.reciprocal.overlap",
        "clinvar.dbVAR.SV", "clinvar.dbVAR.id.variant", "clinvar.dbVAR.clinical", "clinvar.dbVAR.id.clinical",
        "clinvar.reciprocal.overlap", "tool", "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category"
    )
    dm.idx  <- grep("^DepMap\\.", def.col.orig)
    def.col <- c(def.col.orig[1:(dm.idx[1] - 1)], depmap.cols,
                 def.col.orig[(dm.idx[length(dm.idx)] + 1):length(def.col.orig)])

    # Robustezza subset: su uno slot a 0 righe (tool escluso -> placeholder),
    # annotate.results6.list.vGene riduce lo slot a una tabella senza colonna
    # 'symbol' (auto.annotate restituisce NULL per ogni svtype -> rbindlist
    # svuota). Il raggruppamento by=symbol fallirebbe su colonna assente.
    # Guard identico al comportamento originale: per ogni slot REALE (che ha
    # 'symbol', come nel run default) la chiamata e' verbatim invariata; per
    # uno slot vuoto si restituisce una data.table vuota (-> poi NULL nel blocco
    # di formattazione).
    .group.sym <- function(dt) {
        if ("symbol" %in% colnames(dt) && nrow(dt) > 0)
            dt[, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
        else
            data.table()
    }
    gasoline.out  <- .group.sym(x[[1]])
    cutesv.out    <- .group.sym(x[[3]])
    sniffles2.out <- .group.sym(x[[4]])
    ng100k.out    <- .group.sym(x[[5]])

    if (nrow(data.frame(gasoline.out)) > 0) {
        colnames.gasoline.out <- c("tag", "seqnames", "start", "end", "Size", "SV.Type", "Clustered.Signatures", "Total.Reads", "symbol", def.col)
        gasoline.out.formatted <- gasoline.out[, ..colnames.gasoline.out]
        colnames(gasoline.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "Total.Reads", "symbol", def.col)
    } else { gasoline.out.formatted <- NULL }

    if (nrow(data.frame(cutesv.out)) > 0) {
        colnames.cutesv.out <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
        cutesv.out.formatted <- cutesv.out[, ..colnames.cutesv.out]
        colnames(cutesv.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
    } else { cutesv.out.formatted <- NULL }

    if (nrow(data.frame(sniffles2.out)) > 0) {
        colnames.sniffles2.out <- c("tag", "seqnames", "start", "end", "length", "svtype", "SUPPORT", "COVERAGE", "VAF", "hap", "symbol", def.col)
        sniffles2.out.formatted <- sniffles2.out[, ..colnames.sniffles2.out]
        colnames(sniffles2.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "coverage", "VAF", "hap", "symbol", def.col)
    } else { sniffles2.out.formatted <- NULL }

    if (nrow(data.frame(ng100k.out)) > 0) {
        colnames.ng100k.out <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
        ng100k.out.formatted <- ng100k.out[, ..colnames.ng100k.out]
        colnames(ng100k.out.formatted) <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    } else { ng100k.out.formatted <- NULL }

    combined.dt <- rbindlist(list(gasoline.out.formatted, cutesv.out.formatted, sniffles2.out.formatted, ng100k.out.formatted), fill = TRUE)
    combined.dt <- combined.dt[, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    # NB: quando un tool e' totalmente escluso (0 righe) alcune colonne
    # (es. 'tag' che esiste solo in gasoline/cutesv/sniffles, ma non in ng100k)
    # possono mancare. Il verbatim assume che siano presenti; qui mettiamo il
    # guard per non 'rompere' l'annotazione (stesso risultato, colonne assenti
    # semplicemente non vengono normalizzate).
    if ("seqnames" %in% colnames(combined.dt)) combined.dt[, seqnames := gsub("chrchr", "chr", seqnames)]
    if ("tag" %in% colnames(combined.dt))       combined.dt[, tag := gsub("chrchr", "chr", tag)]
    combined.dt
}

## ------------------- create.helping.col.v2 --------------------------------
# Costruisce <TAG>.Sum.Score / <TAG>.Sum.Score.match per ogni tag presente in
# sumscore.cols (CNS/Pediatric sempre, + tag custom). Se un tag non ha colonne
# DepMap, la somma e' 0 (stesso comportamento della verbatim su gruppo assente).
create.helping.col.v2 <- function(input, sumscore.cols) {
    tags <- unique(sub("\\.Sum\\.Score(\\.match)?$", "", sumscore.cols))
    for (t in tags) {
        cur <- input[, grep(paste0("DepMap\\.", t), colnames(input)), with = FALSE]
        score <- tryCatch({
            if (ncol(cur) == 0) rep(0, nrow(input))
            else apply(cur, 1, function(x) sum(as.numeric(x), na.rm = TRUE))
        }, error = function(e) rep(0, nrow(input)))
        score.col <- paste0(t, ".Sum.Score")
        match.col <- paste0(t, ".Sum.Score.match")
        input[[score.col]] <- score
        # anti-rotture: quando un tool e' escluso, 'sv'/'CNV.Type' possono mancare
        sv    <- if (!is.null(input$sv))       input$sv       else rep(NA_character_, nrow(input))
        cnvt  <- if (!is.null(input$CNV.Type)) input$CNV.Type else rep(NA_character_, nrow(input))
        m <- rep(NA, nrow(input))
        dupmask <- (grepl("DUP", sv) | grepl("TDup", sv) | grepl("uplicat", cnvt))
        dupmask[is.na(dupmask)] <- FALSE
        m[score > 0 & dupmask] <- 1
        m[score < 0 & dupmask] <- -1
        input[[match.col]] <- m
    }
    input
}

## ------------------- final.cleaning.v2 ------------------------------------
# Identico a final.cleaning (core) tranne che il blocco DepMap e il blocco
# sum-score sono sostituiti da `depmap.cols` e `sumscore.cols` (dinamici).
# Il blocco non-DepMap e' IDENTICO (stessi nomi, stesso ordine, stesse rename).
final.cleaning.v2 <- function(x, depmap.cols, sumscore.cols) {
    sel.base <- c(
        "symbol", "seqnames", "start", "end", "length", "sv", "CNV.Type", "Shift", "Copies", "reads", "hap",
        "counts", "tot.score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT",
        "DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max",
        "DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max",
        "DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max",
        "decipher.#population_cnv_id", "decipher.deletion_frequency", "decipher.duplication_frequency", "decipher.frequency",
        "decipher.reciprocal.overlap", "gAD.CNV.sv.CNV", "gAD.CNV.symbol.CNV", "gAD.CNV.SF", "gAD.CNV.SF_nfe", "gAD.CNV.reciprocal.overlap",
        "gnomAD.svtype", "gnomAD.AF", "gnomAD.FREQ_HOMREF", "gnomAD.FREQ_HET", "gnomAD.FREQ_HOMALT",
        "gnomAD.AF_nfe", "gnomAD.FREQ_HOMREF_nfe", "gnomAD.FREQ_HET_nfe", "gAD.reciprocal.overlap", "clinvar.dbVAR.SV",
        "clinvar.dbVAR.id.variant", "clinvar.dbVAR.clinical", "clinvar.dbVAR.id.clinical", "clinvar.reciprocal.overlap",
        "tool", "CNS.Sum.Score", "CNS.Sum.Score.match", "Pediatric.Sum.Score", "Pediatric.Sum.Score.match",
        "i.TF_string", "i.n_TF", "onkoDB.Gene Type", "macrocategories", "hallmarks",
        "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category", "coverage"
    )
    ren.base <- c(
        "symbol", "seqnames", "start", "end", "length", "sv", "CNV.Type", "Shift", "Copies", "reads", "hap",
        "String.Interactions", "String.Interactions.Score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT",
        "DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max",
        "DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max",
        "DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max",
        "decipher.population_cnv_id", "decipher.deletion_frequency", "decipher.duplication_frequency", "decipher.frequency",
        "decipher.reciprocal.overlap", "gAD.CNV.sv.CNV", "gAD.CNV.symbol.CNV", "gAD.CNV.SF", "gAD.CNV.SF_nfe", "gAD.CNV.reciprocal.overlap",
        "gnomAD.svtype", "gnomAD.AF", "gnomAD.FREQ_HOMREF", "gnomAD.FREQ_HET", "gnomAD.FREQ_HOMALT",
        "gnomAD.AF_nfe", "gnomAD.FREQ_HOMREF_nfe", "gnomAD.FREQ_HET_nfe", "gAD.reciprocal.overlap", "clinvar.dbVAR.SV",
        "clinvar.dbVAR.id.variant", "clinvar.dbVAR.clinical", "clinvar.dbVAR.id.clinical", "clinvar.reciprocal.overlap",
        "tool", "CNS.Sum.Score", "CNS.Sum.Score.match", "Pediatric.Sum.Score", "Pediatric.Sum.Score.match",
        "TF.targets", "TF.targets.number", "onkoDB.Gene Type", "macrocategories", "hallmarks",
        "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category", "coverage"
    )
    # sostituzione blocco DepMap
    dm.idx <- grep("^DepMap\\.", sel.base)
    sel <- c(sel.base[1:(dm.idx[1] - 1)], depmap.cols, sel.base[(dm.idx[length(dm.idx)] + 1):length(sel.base)])
    ren <- c(ren.base[1:(dm.idx[1] - 1)], depmap.cols, ren.base[(dm.idx[length(dm.idx)] + 1):length(ren.base)])
    # sostituzione blocco sum-score (ricerca DOPO il shift del blocco DepMap)
    ss.std <- c("CNS.Sum.Score", "CNS.Sum.Score.match", "Pediatric.Sum.Score", "Pediatric.Sum.Score.match")
    ss.idx <- which(sel %in% ss.std)
    if (length(ss.idx) > 0) {
        sel <- c(sel[1:(ss.idx[1] - 1)], sumscore.cols, sel[(ss.idx[length(ss.idx)] + 1):length(sel)])
        ren <- c(ren[1:(ss.idx[1] - 1)], sumscore.cols, ren[(ss.idx[length(ss.idx)] + 1):length(ren)])
    }

    # anti-rotture: quando un tool e' escluso alcune colonne di `sel` possono
    # non esistere (es. 'reads'/'hap'/'coverage' che provengono da gasoline/
    # cutesv/sniffles e non da ng100k). Le aggiungiamo come NA prima della
    # selezione, cosi' lo schema resta completo e l'annotazione non si rompe.
    missing.cols <- setdiff(sel, colnames(x))
    if (length(missing.cols) > 0) {
        x[, (missing.cols) := lapply(missing.cols, function(cc) rep(NA_character_, nrow(x)))]
    }
    x.2 <- cbind(x[, ..sel])
    to.add <- c(grep("specific", colnames(x)), grep("common", colnames(x)))
    if (length(to.add) > 0) {
        x.2 <- cbind(x.2, x[, ..to.add])
        to.add.names <- c(grep("specific", colnames(x), value = T), grep("common", colnames(x), value = T))
        colnames(x.2) <- c(ren, to.add.names)
    } else {
        colnames(x.2) <- ren
    }

    ind <- grep("DepMap", colnames(x.2))
    if (length(ind) > 0) {
        x.2[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]
    }
    x.2
}
