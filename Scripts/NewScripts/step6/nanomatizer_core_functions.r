final.cleaning <- function(x){
    x.2 <- cbind(x[,c("symbol", "seqnames", "start", "end", "length", "sv", "CNV.Type", "Shift", "Copies", "reads", "hap", "counts", "tot.score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT", 
        "DepMap.CNS.Min", "DepMap.CNS.Q1",  "DepMap.CNS.Median",  "DepMap.CNS.Q3",  "DepMap.CNS.Max",  "DepMap.Pediatric.Min",  "DepMap.Pediatric.Q1",  "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3",  "DepMap.Pediatric.Max",  "DepMap.MBL.Min",  "DepMap.MBL.Q1",  "DepMap.MBL.Median",  "DepMap.MBL.Q3",  "DepMap.MBL.Max",  "decipher.#population_cnv_id", "decipher.deletion_frequency", "decipher.duplication_frequency", "decipher.frequency",
        "decipher.reciprocal.overlap",  "gAD.CNV.sv.CNV",  "gAD.CNV.symbol.CNV",  "gAD.CNV.SF",  "gAD.CNV.SF_nfe",  "gAD.CNV.reciprocal.overlap",  "gnomAD.svtype",  "gnomAD.AF",  "gnomAD.FREQ_HOMREF",  "gnomAD.FREQ_HET",  "gnomAD.FREQ_HOMALT",  "gnomAD.AF_nfe",  "gnomAD.FREQ_HOMREF_nfe",  "gnomAD.FREQ_HET_nfe",  "gAD.reciprocal.overlap",  "clinvar.dbVAR.SV",
        "clinvar.dbVAR.id.variant",  "clinvar.dbVAR.clinical",  "clinvar.dbVAR.id.clinical",  "clinvar.reciprocal.overlap",  "tool",   "CNS.Sum.Score",  "CNS.Sum.Score.match",  "Pediatric.Sum.Score",  "Pediatric.Sum.Score.match",  "i.TF_string", "i.n_TF", "onkoDB.Gene Type", "macrocategories", "hallmarks", "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category", "coverage"
    )])
    to.add <- c(grep("specific", colnames(x)), grep("common", colnames(x)))
    if(length(to.add)>0) {
        x.2 <- cbind(x.2, x[,..to.add])
        to.add.names <- c(grep("specific", colnames(x), value = T), grep("common", colnames(x), value = T))
        colnames(x.2) <- c("symbol", "seqnames", "start", "end", "length", "sv", "CNV.Type", "Shift", "Copies", "reads", "hap", "String.Interactions", "String.Interactions.Score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT", 
            "DepMap.CNS.Min", "DepMap.CNS.Q1",  "DepMap.CNS.Median",  "DepMap.CNS.Q3",  "DepMap.CNS.Max",  "DepMap.Pediatric.Min",  "DepMap.Pediatric.Q1",  "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3",  "DepMap.Pediatric.Max",  "DepMap.MBL.Min",  "DepMap.MBL.Q1",  "DepMap.MBL.Median",  "DepMap.MBL.Q3",  "DepMap.MBL.Max",  "decipher.population_cnv_id", "decipher.deletion_frequency", "decipher.duplication_frequency", "decipher.frequency",
            "decipher.reciprocal.overlap",  "gAD.CNV.sv.CNV",  "gAD.CNV.symbol.CNV",  "gAD.CNV.SF",  "gAD.CNV.SF_nfe",  "gAD.CNV.reciprocal.overlap",  "gnomAD.svtype",  "gnomAD.AF",  "gnomAD.FREQ_HOMREF",  "gnomAD.FREQ_HET",  "gnomAD.FREQ_HOMALT",  "gnomAD.AF_nfe",  "gnomAD.FREQ_HOMREF_nfe",  "gnomAD.FREQ_HET_nfe",  "gAD.reciprocal.overlap",  "clinvar.dbVAR.SV",
            "clinvar.dbVAR.id.variant",  "clinvar.dbVAR.clinical",  "clinvar.dbVAR.id.clinical",  "clinvar.reciprocal.overlap",  "tool", "CNS.Sum.Score",  "CNS.Sum.Score.match",  "Pediatric.Sum.Score",  "Pediatric.Sum.Score.match",  "TF.targets", "TF.targets.number", "onkoDB.Gene Type", "macrocategories", "hallmarks", "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category", "coverage", to.add.names
        )        
    } else {
        colnames(x.2) <- c("symbol", "seqnames", "start", "end", "length", "sv", "CNV.Type", "Shift", "Copies", "reads", "hap", "String.Interactions", "String.Interactions.Score", "HAPLOINSUFFICIENCY", "TRIPLOSENSITIVITY", "ONLINE REPORT", 
            "DepMap.CNS.Min", "DepMap.CNS.Q1",  "DepMap.CNS.Median",  "DepMap.CNS.Q3",  "DepMap.CNS.Max",  "DepMap.Pediatric.Min",  "DepMap.Pediatric.Q1",  "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3",  "DepMap.Pediatric.Max",  "DepMap.MBL.Min",  "DepMap.MBL.Q1",  "DepMap.MBL.Median",  "DepMap.MBL.Q3",  "DepMap.MBL.Max",  "decipher.population_cnv_id", "decipher.deletion_frequency", "decipher.duplication_frequency", "decipher.frequency",
            "decipher.reciprocal.overlap",  "gAD.CNV.sv.CNV",  "gAD.CNV.symbol.CNV",  "gAD.CNV.SF",  "gAD.CNV.SF_nfe",  "gAD.CNV.reciprocal.overlap",  "gnomAD.svtype",  "gnomAD.AF",  "gnomAD.FREQ_HOMREF",  "gnomAD.FREQ_HET",  "gnomAD.FREQ_HOMALT",  "gnomAD.AF_nfe",  "gnomAD.FREQ_HOMREF_nfe",  "gnomAD.FREQ_HET_nfe",  "gAD.reciprocal.overlap",  "clinvar.dbVAR.SV",
            "clinvar.dbVAR.id.variant",  "clinvar.dbVAR.clinical",  "clinvar.dbVAR.id.clinical",  "clinvar.reciprocal.overlap",  "tool",  "CNS.Sum.Score",  "CNS.Sum.Score.match",  "Pediatric.Sum.Score",  "Pediatric.Sum.Score.match",  "TF.targets", "TF.targets.number", "onkoDB.Gene Type", "macrocategories", "hallmarks", "ENCODE.blacklist.overlap.pct", "ENCODE.blacklist.category", "coverage"
        )
    }

    # # qua devo splittare per chr e ordinare 
    # ind <- order(sapply(x$tag, function(y){
    #     new.y <- gsub("chr", "", y)
    #     new.y2 <- min(unlist(strsplit(new.y, ', ')))
    #     new.y2
    # }))
    # x.2 <- x.2[ind,]

    ind <- grep("DepMap", colnames(x.2))
    x.2[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]
    x.2
}

annotateODB <- function(x){
    new.merged <- merge(onkoDB[,c(1,7)], x, by.x = 'Hugo Symbol', by.y = "symb", all.y = T)
    new.merged <- cbind(new.merged[,c(-1,-2)], onkoDB = new.merged[,2])
    new.merged
}

annotateTF <- function(x, tf){

    # Assicurati che x sia un data.table
    setDT(x)

    # Creo la colonna symb pulita se non ce l'hai già
    x[, symb := gsub(" \\(X\\)", "", symbol)]

    # Filtro e processo i TF
    tf_filtered <- tf[V2 %in% x$symb]
    tf_filtered[, symbol := fcase(
    V3 == "Activation", "+",
    V3 == "Repression", "-",
    V3 == "Unknown", "?",
    default = "?"
    )]
    tf_filtered[, TF_symbol := paste0(V1, symbol)]

    # Creo la tabella con le stringhe E il conteggio
    tf_summary <- tf_filtered[, .(
    TF_string = paste(TF_symbol, collapse = ", "),
    n_TF = .N
    ), by = V2]

    setnames(tf_summary, "V2", "symb")

    # Faccio il merge con x
    x <- merge(x, tf_summary, by = "symb", all.x = TRUE)

    # Per i geni senza match, metto 0 invece di NA nel conteggio
    x[is.na(n_TF), n_TF := 0]
    x <- tf_summary[x, on = "symb"]
    return(x)
}

annotate.macrocategories <- function(dt) {
    genes <- gsub(" \\(X\\)", "", dt$symbol)
    to.add1 <- unlist(sapply(genes, function(x){
        out <- apply(macrocateg[,c(-4:-1)], 1, function(y){
            match(x, y)
        })
        out <- which(!is.na(out))
        paste0(na.omit(unique(c(macrocateg[out, 4]))), collapse = "; ")
    }, simplify = FALSE))

    to.add2 <- unlist(sapply(genes, function(x){
        out <- apply(macrocateg[,c(-5:-1)], 1, function(y){
            match(x, y)
        })
        out <- which(!is.na(out))
        paste0(na.omit(unique((macrocateg[out, 1]))), collapse = "; ")
    }, simplify = FALSE))

    dt$macrocategories <- to.add1
    dt$hallmarks <- to.add2
    dt
}

annotate.macrocat <- function(x){
    unlist(
    sapply(c(genes[1:3], "MYC"), function(x){
        out <- apply(macrocateg[,c(-5:-1)], 1, function(y){
            match(x, y)
        })
        out <- which(!is.na(out))
        macrocateg[out, c(3,4)]
    }, simplify = FALSE))


  # Step 1: Crea mapping gene -> macrocategorie
  gene_macro_map <- list()

  for(i in 1:nrow(macrocateg)) {
    # Estrai tutti i geni da questo pathway (colonne V3 in poi)
    genes <- unlist(macrocateg[i, 3:ncol(macrocateg)])
    genes <- genes[!is.na(genes) & genes != ""]

    # Macrocategorie per questo pathway
    macro1 <- macrocateg[i, "macrocategory_1"]
    macro2 <- macrocateg[i, "macrocategory_2"]

    # Aggiungi macrocategorie per ogni gene
    for(gene in genes) {
      if(!gene %in% names(gene_macro_map)) {
        gene_macro_map[[gene]] <- character()
      }
      gene_macro_map[[gene]] <- c(gene_macro_map[[gene]], macro1)
      if(!is.na(macro2)) {
        gene_macro_map[[gene]] <- c(gene_macro_map[[gene]], macro2)
      }
    }
  }
  
  # Fai unique per ogni gene (massimo 2 categorie diverse)
  gene_macro_map <- lapply(gene_macro_map, unique)
  
  # Step 2: Aggiungi colonne al dt
  dt[, Macrocategory_1 := sapply(symbol, function(s) {
    if(s %in% names(gene_macro_map) && length(gene_macro_map[[s]]) >= 1) {
      gene_macro_map[[s]][1]
    } else {
      NA_character_
    }
  })]
  
  dt[, Macrocategory_2 := sapply(symbol, function(s) {
    if(s %in% names(gene_macro_map) && length(gene_macro_map[[s]]) >= 2) {
      gene_macro_map[[s]][2]
    } else {
      NA_character_
    }
  })]
  
  return(dt)
}

get.value.from.vector <- function(values, operation) {
  if (is.null(values) || length(values) == 0 || all(is.na(values))) {
    return(NA_real_)
  }
  num.values <- suppressWarnings(as.numeric(values))
  num.values <- num.values[!is.na(num.values)]
  if (length(num.values) == 0) {
    return(NA_real_)
  }
  if (operation == 'max') {
    return(max(num.values))
  } else if (operation == 'min') {
    return(min(num.values))
  } else {
    stop("L'argomento 'operation' deve essere 'max' o 'min'.")
  }
}

flatten.data.table <- function(input.dt) {
  dt.flat <- copy(input.dt)
  cols.to.flatten <- names(which(sapply(dt.flat, is.list)))
  
  if (length(cols.to.flatten) > 0) {
    dt.flat[, (cols.to.flatten) := lapply(.SD, function(col) {
      unlist(lapply(col, paste, collapse = ", "))
    }), .SDcols = cols.to.flatten]
  }
  
  return(dt.flat)
}

sort.chrom <- function(x){

    dt <- x
    dt[, min.start.for.sort := sapply(start, get.value.from.vector, operation = "min")]
    chromosome.levels <- mixedsort(unique(dt$seqnames))
    dt[, seqnames.for.sort := factor(seqnames, levels = chromosome.levels)]
    x.sorted <- dt[order(seqnames.for.sort, min.start.for.sort)]

    # Rimuovi le colonne temporanee dal nuovo data.table ordinato
    x.sorted[, c("min.start.for.sort", "seqnames.for.sort") := NULL]
    return(x.sorted)
}

clean <- function(column) {
    sapply(column, function(cell.value) {
        # 1. NUOVO CONTROLLO: Gestisce celle che sono NULL o vettori di lunghezza zero.
        #    Questo è il controllo che risolve l'errore "argument is of length zero".
        if (is.null(cell.value) || length(cell.value) == 0) {
            return(NA_character_)
        }

        # 2. CONTROLLO MIGLIORATO PER NA: Gestisce il caso in cui una cella
        #    possa contenere uno o più NA (es. list(NA, NA)).
        if (all(is.na(cell.value))) {
            return(NA_character_)
        }
        
        # Rimuoviamo gli NA espliciti se la cella è un vettore misto (es. c("A", NA))
        cell.value <- cell.value[!is.na(cell.value)]
        
        # Se dopo aver tolto gli NA non resta nulla, trattiamolo come vuoto/NA
        if (length(cell.value) == 0) {
            return(NA_character_)
        }

        # 3. Convertiamo la cella in carattere per sicurezza
        cell.string <- as.character(cell.value)

        # 4. Dividiamo la stringa in un vettore di elementi
        #    Il paste iniziale gestisce il caso in cui cell.value sia un vettore (es. c("A", "B"))
        elements <- unlist(strsplit(paste(cell.string, collapse=", "), ",[ ]*"))

        # 5. Troviamo gli elementi unici e non vuoti
        unique.elements <- unique(elements[elements != ""])

        # 6. Riuniamo gli elementi unici in una singola stringa
        return(paste(unique.elements, collapse = ", "))

    }, USE.NAMES = FALSE)
}

sign.specific.genes.list <- function(input, genes){
    if(!is.null(nrow(input[[1]]))) input[[1]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[2]]))) input[[2]][gene1 %in% genes, gene1 := paste0(gene1, " (X)")]
    if(!is.null(nrow(input[[2]]))) input[[2]][gene2 %in% genes, gene2 := paste0(gene2, " (X)")]
    if(!is.null(nrow(input[[3]]))) input[[3]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[4]]))) input[[4]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[5]]))) input[[5]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    return(input)
}

sign.specific.genes.dt <- function(input, genes){
    if(!is.null(nrow(input))) input[symbol %in% genes, symbol := paste0(symbol, " (X)")]
    return(input)
}

annotate.centrality.v1 <- function(input.list, input.gene){
    for(i in c(1,3,4,5)){
        input.list[[i]] <- input.list[[i]][, !duplicated(names(input.list[[i]])), with = FALSE] #! QUI ELIMINO DELLE COLONNE (HANNO NOMI DUPOLICATI E NON MI SERVONO A NIENTE (...GIUSTO??)
        input.list[[i]][string.abs.high, on = "symbol", tot.score := i.tot.score]
        input.list[[i]][counts, on = "symbol", string.count := i.N]
        input.list[[i]][is.na(string.count), string.count := 0] # sostituisco na con 0
    }
    input.list
}

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
        for(j in 1:length(ind.list)){
            values.to.add <- db[ind.string[j], ]
            input.list[[i]][ind.list[[j]], `:=`(counts = values.to.add$counts, tot.score = values.to.add$tot.score)]
        }
    }
    input.list
}

rename.duplicated.col <- function(input){
    for(i in 1:length(input)){
        ind <- which(duplicated(colnames(input[[i]])))
        if(length(ind) > 0){
            colnames(input[[i]])[ind] <- paste0(colnames(input[[i]])[ind], ".1")
        }
    }
    input
}

annotate.clingen <- function(input.list){
    for(i in c(1,3,4,5)){
        input.list[[i]] <- merge(input.list[[i]], clingen, all.x = TRUE, by ='symbol',  allow.cartesian = TRUE)
    }
    input.list
}

annotate.depmap.cns <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.CNS.Min", "DepMap.CNS.Q1", "DepMap.CNS.Median", "DepMap.CNS.Q3", "DepMap.CNS.Max")
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.cns[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = 1:nrow(input.list[[i]]), .SDcols = "symbol"]
    }
    input.list
}

annotate.depmap.cns.pediatric <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.Pediatric.Min", "DepMap.Pediatric.Q1", "DepMap.Pediatric.Median", "DepMap.Pediatric.Q3", "DepMap.Pediatric.Max")
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.cns.pediatric[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = 1:nrow(input.list[[i]]), .SDcols = "symbol"]
    }
    input.list
}

annotate.depmap.cns.mbl <- function(input.list){
    for(i in c(1,3,4,5)){
        new_cols <- c("DepMap.MBL.Min", "DepMap.MBL.Q1", "DepMap.MBL.Median", "DepMap.MBL.Q3", "DepMap.MBL.Max")
        input.list[[i]][, (new_cols) := {
            current_symbol <- .SD$symbol
            as.list(quantile(depmap.mbl[[current_symbol]], probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE))
        }, by = 1:nrow(input.list[[i]]), .SDcols = "symbol"]
    }
    input.list
}

auto.annotate <- function(x, sv.label, sv.type, db, db.tag, seqnames.col, start.col, end.col){
    sv.label.pos <- match(sv.label, colnames(x))

    x.sub <- x[x[[sv.label.pos]] == sv.type]

    if(sum(x[[sv.label.pos]] == sv.type) == 0) { 
        print(paste0("non ci sono ", sv.type," da annotare --> INTERROMPO"))
        out <- NULL
        return(out)
    } else {
        db.colnames <- colnames(db)
        templato <- db[0,]
        x.sub <- x[x[[sv.label.pos]] == sv.type]
        x.sub[, (db.colnames) := templato]         # cosi ogni colonna e della classe che mi serve, senno mi sbrocca
        x.sub[,(paste0(db.tag, ".reciprocal.overlap")) := NA_real_] # NA per valori numerici
        col.pos <- c(match(seqnames.col, colnames(x)), match(start.col, colnames(x)), match(end.col, colnames(x)))
        x.sub[, tag := do.call(paste, c(.SD, sep = ".")), .SDcols = col.pos]

        hits <- findOverlaps(makeGRangesFromDataFrame(x.sub), makeGRangesFromDataFrame(data.frame(db)))
        if(length(hits) > 0 ){
            gr1 <- makeGRangesFromDataFrame(x.sub[queryHits(hits),])
            gr2 <- makeGRangesFromDataFrame(db[subjectHits(hits)])
            intersection <- pintersect(gr1, gr2)
            length.intersection <- width(intersection)
            width.gr1 <- width(gr1)
            width.gr2 <- width(gr2)
            reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
            reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
            reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2, 3)

            x.sub.to.add <- x.sub[setdiff(seq_len(nrow(x.sub)), queryHits(hits)),]
            x.sub.annotated <- x.sub[queryHits(hits),]
            x.sub.annotated[ , (paste0(db.tag, ".reciprocal.overlap")) := reciprocal.overlap.mean]
            x.sub.annotated[ , (db.colnames) := db[subjectHits(hits),]]
            x.sub.to.add[ , (paste0(db.tag, ".reciprocal.overlap")) := NA]
            x.sub.to.add[ , (db.colnames) := NA]
            out <- rbind(x.sub.annotated, x.sub.to.add)
        } else {
            print("non mi trova hits in gasoline del gAD")
            out <- x.sub # se non lo crea allora mettilo come x.del.gdv (x.del che potrebbe essere come x.del se anche prima non ha trovato nulla) (mi serve per gAD)
        }

        # qui aggrego (prendo solo la prima, sticazzi delle varie annoatzioni che non ci servono) cio che ha un reciprocal overlap basso, senno ho troppe righe
        out.to.add <- out[out[[(paste0(db.tag, ".reciprocal.overlap"))]] < 0.7]
        col.overlap <- paste0(db.tag, ".reciprocal.overlap")
        out.to.add <- out.to.add[, .SD[which.max(get(col.overlap))], by = .(tag, symbol)]
        out <- rbind(out[out[[(paste0(db.tag, ".reciprocal.overlap"))]] >= 0.7], out.to.add, out[is.na(out[[(paste0(db.tag, ".reciprocal.overlap"))]])])

        return(out)
    }
}

collapse.genes <- function(x, overlap.names){

    def.col <- c(
        "counts",
        "tot.score",
        "HAPLOINSUFFICIENCY",
        "TRIPLOSENSITIVITY",
        "ONLINE REPORT",
        "DepMap.CNS.Min",
        "DepMap.CNS.Q1",
        "DepMap.CNS.Median",
        "DepMap.CNS.Q3",
        "DepMap.CNS.Max",
        "DepMap.Pediatric.Min",
        "DepMap.Pediatric.Q1",
        "DepMap.Pediatric.Median",
        "DepMap.Pediatric.Q3",
        "DepMap.Pediatric.Max",
        "DepMap.MBL.Min",
        "DepMap.MBL.Q1",
        "DepMap.MBL.Median",
        "DepMap.MBL.Q3",
        "DepMap.MBL.Max",
        "decipher.#population_cnv_id",
        "decipher.deletion_observations",
        "decipher.deletion_frequency",
        "decipher.duplication_observations",
        "decipher.duplication_frequency",
        "decipher.observations",
        "decipher.frequency",
        "decipher.type",
        "decipher.sample_size",
        "decipher.study",
        "decipher.reciprocal.overlap",
        "gAD.CNV.sv.CNV",
        "gAD.CNV.symbol.CNV",
        "gAD.CNV.SF",
        "gAD.CNV.SF_nfe",
        "gAD.CNV.reciprocal.overlap",
        "gnomAD.svtype",
        "gnomAD.AF",
        "gnomAD.FREQ_HOMREF",
        "gnomAD.FREQ_HET",
        "gnomAD.FREQ_HOMALT",
        "gnomAD.AF_nfe",
        "gnomAD.FREQ_HOMREF_nfe",
        "gnomAD.FREQ_HET_nfe",
        "gAD.reciprocal.overlap",
        "clinvar.dbVAR.SV",
        "clinvar.dbVAR.id.variant",
        "clinvar.dbVAR.clinical",
        "clinvar.dbVAR.id.clinical",
        "clinvar.reciprocal.overlap",
        "tool",
        "ENCODE.blacklist.overlap.pct",
        "ENCODE.blacklist.category"
    )

    gasoline.out <- x[[1]][, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    cutesv.out <- x[[3]][, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    sniffles2.out <- x[[4]][, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    ng100k.out <- x[[5]][, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    # ng10k.out <- x[[6]][, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]

    if(nrow(data.frame(gasoline.out)) > 0){
        colnames.gasoline.out <- c("tag", "seqnames", "start", "end", "Size", "SV.Type", "Clustered.Signatures", "Total.Reads", "symbol", def.col)
        gasoline.out.formatted <- gasoline.out[, ..colnames.gasoline.out]
        colnames(gasoline.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "Total.Reads", "symbol", def.col)
    } else {
        gasoline.out.formatted <- NULL
    }

    if(nrow(data.frame(cutesv.out)) > 0){
        colnames.cutesv.out <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
        cutesv.out.formatted <- cutesv.out[, ..colnames.cutesv.out]
        colnames(cutesv.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
    } else {
        cutesv.out.formatted <- NULL
    }

    if(nrow(data.frame(sniffles2.out)) > 0){
        colnames.sniffles2.out <- c("tag", "seqnames", "start", "end", "length", "svtype", "SUPPORT", "COVERAGE", "VAF", "hap", "symbol", def.col)
        sniffles2.out.formatted <- sniffles2.out[, ..colnames.sniffles2.out]
        colnames(sniffles2.out.formatted) <-  c("tag", "seqnames", "start", "end", "length", "sv", "reads", "coverage", "VAF", "hap", "symbol", def.col)
    } else {
        sniffles2.out.formatted <- NULL
    }

    if(nrow(data.frame(ng100k.out)) > 0){
        colnames.ng100k.out <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
        ng100k.out.formatted <- ng100k.out[, ..colnames.ng100k.out]
        colnames(ng100k.out.formatted) <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    } else {
        ng100k.out.formatted <- NULL
    }

    # if(nrow(data.frame(ng10k.out)) > 0){
    #     colnames.ng10k.out <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    #     ng10k.out.formatted <- ng10k.out[, ..colnames.ng10k.out]
    #     colnames(ng10k.out.formatted) <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    # } else {
    #     ng10k.out.formatted <- NULL
    # }

    print("ultimo")
    combined.dt <- rbindlist(list(gasoline.out.formatted , cutesv.out.formatted, sniffles2.out.formatted, ng100k.out.formatted), fill = TRUE)
    # combined.dt[, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]
    combined.dt <- combined.dt[, lapply(.SD, function(x) paste(x, collapse = ",")), by = symbol]

    combined.dt[, seqnames := gsub("chrchr", "chr", seqnames)]
    combined.dt[, tag := gsub("chrchr", "chr", tag)]
    combined.dt
}

sign.specific.genes <- function(input, genes){
    if(!is.null(nrow(input[[1]]))) input[[1]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[2]]))) input[[2]][gene1 %in% genes, gene1 := paste0(gene1, " (X)")]
    if(!is.null(nrow(input[[2]]))) input[[2]][gene2 %in% genes, gene2 := paste0(gene2, " (X)")]
    if(!is.null(nrow(input[[3]]))) input[[3]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[4]]))) input[[4]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    if(!is.null(nrow(input[[5]]))) input[[5]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    # if(!is.null(nrow(input[[6]]))) input[[6]][symbol %in% genes, symbol := paste0(symbol, " (X)")]
    return(input)
}

convert.list.to.df <- function(input, rm.int){
    out <- data.table(
        seqnames = c(
            input[[1]]$Chr, 
            apply(cbind(data.frame(input[[2]])[,1], data.frame(input[[2]][,2])), 1, function(x) {paste0(x[1], "-", x[2])}),
            input[[3]]$seqnames,
            input[[4]]$seqnames,
            input[[5]]$seqnames
            # input[[6]]$seqnames
        ),
        start = c(
            input[[1]]$start, 
            apply(cbind(data.frame(input[[2]])[,3], data.frame(input[[2]][,3])), 1, function(x) {paste0(x[1], "-", x[2])}),
            input[[3]]$start,
            input[[4]]$start,
            input[[5]]$start
            # input[[6]]$start
        ),
        end = c(
            input[[1]]$end, 
            apply(cbind(data.frame(input[[2]])[,3], data.frame(input[[2]][,3])), 1, function(x) {paste0(x[1], "-", x[2])}),
            input[[3]]$end,
            input[[4]]$end,
            input[[5]]$end
            # input[[6]]$end
        ),
        sv = c(
            input[[1]]$SVLabelMostProb,
            rep("Intertrans", nrow(input[[2]])),
            input[[3]]$sv,
            input[[4]]$svtype,
            rep("NG100K", nrow(input[[5]]))
            # rep("NG10K", nrow(input[[6]]))
        ),
        sheet = c(
            rep(1, nrow(input[[1]])),
            rep(2, nrow(input[[2]])),
            rep(3, nrow(input[[3]])),
            rep(4, nrow(input[[4]])),
            rep(5, nrow(input[[5]]))
            # rep(6, nrow(input[[6]]))
        ),
        tag = c(
            input[[1]]$tag,
            input[[2]]$tag,
            input[[3]]$tag,
            input[[4]]$tag,
            input[[5]]$tag
            # input[[6]]$tag
        )
    )
    out <- out[grep("chr", seqnames, invert = T), seqnames := paste0("chr", seqnames)]
    if(rm.int) out <- subset(out, !is.na(as.numeric(out$start)))
    return(out)
}

filtra.per.overlap.reciproco.vold1 <- function(lista.granges, soglia.overlap = 0.9) {
    #' @title Trova Range Partner Congiunti basato su Reciprocal Overlap
    #' @description Data una lista di N set di range, questa funzione identifica,
    #'   per ogni set, i range che hanno un "partner forte" (con reciprocal
    #'   overlap > soglia) in OGNUNO degli altri set.
    #'
    #' @param lista.granges Un oggetto GRangesList contenente 2 o più set di range.
    #' @param soglia.overlap Soglia (0-1) per il reciprocal overlap mean. Default: 0.7.
    #'
    #' @return Un oggetto GRangesList con la stessa struttura dell'input. Ciascun
    #'   elemento contiene i range originali che hanno soddisfatto il criterio
    #'   di avere un partner forte in tutti gli altri set.
    #'

    # --- 1. Validazione e Setup Iniziale ---
    if (!is(lista.granges, "GRangesList") || length(lista.granges) < 2) {
        stop("L'input deve essere un oggetto GRangesList contenente almeno 2 elementi.")
    }
    if (is.null(names(lista.granges))) {
        names(lista.granges) <- paste0("set_", seq_along(lista.granges))
    }
    set_names <- names(lista.granges)

    # --- 2. Calcolo Pairwise degli Overlap Forti ---
    # Trova, per ogni coppia di set, quali range hanno un partner forte.

    # Struttura per salvare gli indici: pairwise_indices$setA$setB = indici di A forti contro B
    pairwise_indices <- stats::setNames(lapply(set_names, function(x) list()), set_names)
    
    # Genera tutte le coppie uniche di set
    pairs <- utils::combn(set_names, 2, simplify = FALSE)

    cat(sprintf("INFO: Analisi di %d coppie di set...\n", length(pairs)))

    for (pair in pairs) {
        name1 <- pair[1]
        name2 <- pair[2]
        set1 <- lista.granges[[name1]]
        set2 <- lista.granges[[name2]]

        hits <- findOverlaps(set1, set2)
        if (length(hits) == 0) {
            # Se non c'è overlap, gli indici forti sono un set vuoto
            pairwise_indices[[name1]][[name2]] <- integer(0)
            pairwise_indices[[name2]][[name1]] <- integer(0)
            next
        }

        # Calcolo del reciprocal overlap per la coppia
        overlapping_1 <- set1[queryHits(hits)]
        overlapping_2 <- set2[subjectHits(hits)]
        intersection_width <- width(pintersect(overlapping_1, overlapping_2))
        overlap1 <- intersection_width / width(overlapping_1)
        overlap2 <- intersection_width / width(overlapping_2)
        overlap1[is.nan(overlap1)] <- 0
        overlap2[is.nan(overlap2)] <- 0
        reciprocal_mean <- (overlap1 + overlap2) / 2
        
        strong_hits <- hits[reciprocal_mean > soglia.overlap]
        
        # Salva gli indici dei range che hanno superato la soglia per questa coppia
        pairwise_indices[[name1]][[name2]] <- unique(queryHits(strong_hits))
        pairwise_indices[[name2]][[name1]] <- unique(subjectHits(strong_hits))
    }

    # --- 3. Intersezione degli Indici per Trovare i Partner Congiunti ---
    # Per ogni set, trova i range che hanno un partner forte in TUTTI gli altri set.

    risultati_filtrati <- lapply(set_names, function(current_set_name) {
        original_ranges <- lista.granges[[current_set_name]]
        
        # Se c'è un solo altro set, la logica è semplice
        if (length(set_names) == 2) {
            other_name <- set_names[set_names != current_set_name]
            final_indices_to_keep <- pairwise_indices[[current_set_name]][[other_name]]
            return(original_ranges[final_indices_to_keep])
        }
        
        # Nomi degli altri set con cui fare il confronto
        other_sets <- set_names[set_names != current_set_name]
        
        # Raccogli le liste di indici forti per il set corrente contro tutti gli altri
        list_of_strong_indices <- lapply(other_sets, function(other_name) {
            pairwise_indices[[current_set_name]][[other_name]]
        })
        
        # Calcola l'intersezione: un range deve essere forte contro TUTTI gli altri
        final_indices_to_keep <- Reduce(intersect, list_of_strong_indices)
        
        return(original_ranges[final_indices_to_keep])
    })

    # --- 4. Restituzione dei Set Filtrati ---
    names(risultati_filtrati) <- set_names
    return(GRangesList(risultati_filtrati))
}

filtra.per.overlap.reciproco <- function(lista.granges, soglia.overlap = 0.9) {
    #' @title Trova Range Partner Congiunti basato su Reciprocal Overlap
    
    # --- 1. Validazione e Setup Iniziale ---
    if (!is(lista.granges, "GRangesList") || length(lista.granges) < 2) {
        if(any(lengths(lista.granges) == 0)) return(lista.granges) 
        stop("L'input deve essere un oggetto GRangesList con almeno 2 elementi.")
    }
    if (is.null(names(lista.granges))) {
        names(lista.granges) <- paste0("set_", seq_along(lista.granges))
    }
    set_names <- names(lista.granges)

    # --- 2. Calcolo Pairwise degli Overlap Forti ---
    pairwise_indices <- stats::setNames(lapply(set_names, function(x) list()), set_names)
    pairs <- utils::combn(set_names, 2, simplify = FALSE)

    for (pair in pairs) {
        name1 <- pair[1]; name2 <- pair[2]
        set1 <- lista.granges[[name1]]; set2 <- lista.granges[[name2]]
        
        hits <- findOverlaps(set1, set2)
        if (length(hits) == 0) {
            pairwise_indices[[name1]][[name2]] <- integer(0)
            pairwise_indices[[name2]][[name1]] <- integer(0)
            next
        }

        overlapping_1 <- set1[queryHits(hits)]
        overlapping_2 <- set2[subjectHits(hits)]
        intersection_width <- width(pintersect(overlapping_1, overlapping_2))
        overlap1 <- intersection_width / width(overlapping_1)
        overlap2 <- intersection_width / width(overlapping_2)
        overlap1[is.nan(overlap1)] <- 0; overlap2[is.nan(overlap2)] <- 0
        reciprocal_mean <- (overlap1 + overlap2) / 2
        
        strong_hits <- hits[reciprocal_mean > soglia.overlap]
        
        pairwise_indices[[name1]][[name2]] <- unique(queryHits(strong_hits))
        pairwise_indices[[name2]][[name1]] <- unique(subjectHits(strong_hits))
    }

    # --- 3. Intersezione degli Indici per Trovare i Partner Congiunti ---
    risultati_filtrati <- lapply(set_names, function(current_set_name) {
        original_ranges <- lista.granges[[current_set_name]]
        
        if (length(set_names) == 2) {
            other_name <- set_names[set_names != current_set_name]
            final_indices_to_keep <- pairwise_indices[[current_set_name]][[other_name]]
            return(original_ranges[final_indices_to_keep])
        }
        
        other_sets <- set_names[set_names != current_set_name]
        list_of_strong_indices <- lapply(other_sets, function(other_name) {
            pairwise_indices[[current_set_name]][[other_name]]
        })
        
        final_indices_to_keep <- Reduce(intersect, list_of_strong_indices)
        return(original_ranges[final_indices_to_keep])
    })
    
    names(risultati_filtrati) <- set_names
    return(GRangesList(risultati_filtrati))
}

analisi.per.blocchi <- function(lista.granges.completa, funzione.analisi, block.size = 1000000) {
    #' @title Esegue un'analisi genomica per blocchi per controllare l'uso della memoria.
    
    # --- ECCO LA RIGA CORRETTA ---
    # Questo metodo è robusto e più leggibile per ottenere i nomi unici dei cromosomi
    all.chromosomes <- as.character(unique(seqnames(unlist(lista.granges.completa))))
    
    final.results.list <- list()

    cat(sprintf("INFO: Inizio analisi per %d cromosomi, usando blocchi da %.1f Mb...\\n", 
                length(all.chromosomes), block.size / 1e6))

    for (chr in all.chromosomes) {
        cat(sprintf("  -> Processando chr: %s...\\n", chr))
        
        grl.chr <- endoapply(lista.granges.completa, function(gr) gr[seqnames(gr) == chr])
        
        max.coord <- max(unlist(lapply(grl.chr, function(gr) if(length(gr) > 0) end(gr) else 0)))
        if (is.infinite(max.coord) || max.coord == 0) next
        
        for (start.pos in seq(1, max.coord, by = block.size)) {
            end.pos <- start.pos + block.size - 1
            current.block <- GRanges(seqnames = chr, ranges = IRanges(start = start.pos, end = end.pos))
            grl.block <- endoapply(grl.chr, function(gr) subsetByOverlaps(gr, current.block))
            
            if(any(lengths(grl.block) == 0)) next
            
            result.block <- tryCatch({
                # funzione.analisi(grl.block, soglia.overlap = 0.7)
                funzione.analisi(grl.block)
            }, error = function(e) {
                cat(sprintf("      - Errore nel blocco %s:%d-%d: %s\\n", chr, start.pos, end.pos, e$message))
                return(NULL)
            })
            
            if (!is.null(result.block) && sum(lengths(result.block)) > 0) {
                final.results.list <- append(final.results.list, list(result.block))
            }
        }
    }

    if (length(final.results.list) == 0) {
      cat("ATTENZIONE: Nessun risultato trovato.\\n")
      return(GRangesList())
    }

    cat("INFO: Unione dei risultati di tutti i blocchi...\\n")
    set.names <- names(final.results.list[[1]])
    
    final.granges.by.set <- lapply(set.names, function(name) {
        granges.list.for.name <- lapply(final.results.list, function(res) res[[name]])
        return(unlist(GRangesList(granges.list.for.name)))
    })
    names(final.granges.by.set) <- set.names
    
    final.grl <- GRangesList(final.granges.by.set)
    final.grl.unique <- endoapply(final.grl, unique)
    
    cat("INFO: Analisi completata!\\n")
    return(final.grl.unique)
}

revert.df.to.list <- function(input, s, original.list){
    original.list.revamped <- list(
        gasoline = c(),
        intertrans = c(),
        cutesv = c(),
        sniffles2 = c(),
        NG100k = c()
        # NG10k = c()
    )
    for(i in unique(input[[s]]$sheet)){
        cur.input <- subset(input[[s]], input[[s]]$sheet == i)
        original.list.revamped[[i]] <- subset(original.list[[i]], original.list[[i]]$tag %in% cur.input$tag)
    }
    original.list.revamped[setdiff(1:5, unique(input[[s]]$sheet))] <- ''
    original.list.revamped[[2]] <- original.list[[2]]
    original.list.revamped
}

collapse.ranges.in.list.6 <- function(x){
    cur1 <- x[[1]]
    if(!is.null(nrow(cur1))){
        cur1[, tag := paste0(tag, "-", symbol)]
        cur1 <- cur1[, .SD[which.max(NumReadsClustered)], by = tag]
    }
    x[[1]] <- cur1

    cur2 <- x[[2]]
    if(!is.null(nrow(cur2))){
        cur2[, tag := paste0(tag, "-", gene1, "-", gene2)]
        cur2 <- cur2[, .SD[which.max(Reads)], by = tag]
    }
    x[[2]] <- cur2

    cur3 <- x[[3]]
    if(!is.null(nrow(cur3))){
        cur3[, tag := paste0(tag, "-", symbol)]
        cur3 <- cur3[, .SD[which.max(reads)], by = tag]
    }
    x[[3]] <- cur3

    cur4 <- x[[4]]
    if(!is.null(nrow(cur4))){
        cur4[, tag := paste0(tag, "-", symbol)]
        cur4 <- cur4[, .SD[which.max(SUPPORT)], by = tag]
    }
    x[[4]] <- cur4

    cur5 <- x[[5]]
    if(!is.null(nrow(cur5))){
        cur5[, tag := paste0(tag, "-", symbol)]
        cur5 <- cur5[, .SD[which.max(abs(Shift))], by = tag]
    }
    x[[5]] <- cur5
    
    # cur6 <- x[[6]]
    # if(!is.null(nrow(cur6))){
    #     cur6[, tag := paste0(tag, "-", symbol)]
    #     cur6 <- cur6[, .SD[which.max(abs(Shift))], by = tag]
    # }
    # x[[6]] <- cur6

    x

}

annotate.results6.list.vGene <- function(input){

    ####! ------------------------------------------ GASOLINE
    x <- input[[1]]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Del", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Del", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Del", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Del", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")


    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "SV.Type", sv.type = "TDup", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "SV.Type", sv.type = "TDup", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "SV.Type", sv.type = "TDup", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "SV.Type", sv.type = "TDup", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Inv", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Ins", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "SV.Type", sv.type = "Ins", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)

    gasoline.annotated.all.collapsed <- copy(x.annotated.all)
    gasoline.annotated.all.collapsed$tool <- "Gasoline"


    ####! ------------------------------------------ INTERTRANS
    x <- input[[2]]
    x[, Chr := paste0(seqnames.1, "-", seqnames.2)]
    x[, start := paste0(start.1, "-", start.2)]
    x[, gene := paste0(gene1, "-", gene2)]
    x.intertrans.aggregated <- x
    x.intertrans.aggregated$tools <- "Intertrans"


    ####! ------------------------------------------ CUTESV

    x <- input[[3]]
    if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DUP
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "sv", sv.type = "INV", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "sv", sv.type = "INS", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "sv", sv.type = "INS", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)


    cutesv.annotated.all.collapsed <- copy(x.annotated.all)
    cutesv.annotated.all.collapsed$tool <- "CuteSV"


    ####! ------------------------------------------ SNIFFLES2

    x <- input[[4]]
    colnames(x)[14:18] <- paste0("ens.", colnames(x)[14:18])
    if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DUP
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "svtype", sv.type = "INV", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "svtype", sv.type = "INS", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "svtype", sv.type = "INS", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
   
    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)

    x.annotated.all[, SUPPORT := gsub(";.*", "", gsub(".*SUPPORT=", "", x.annotated.all$info))]
    x.annotated.all[, COVERAGE := gsub(";.*", "", gsub(".*COVERAGE=", "", x.annotated.all$info))]
    x.annotated.all[, VAF := gsub(";.*", "", gsub(".*AF=", "", x.annotated.all$info))]
    x.annotated.all[, hap := gsub(":.*", "", x.annotated.all$sample)]

    sniffles.annotated.all.collapsed <- copy(x.annotated.all)
    sniffles.annotated.all.collapsed$tool <- "Sniffles2"


    ####! ------------------------------------------ NG100K
    if(length(input[[5]]) == 1 && input[[5]] == ''){
        ng100k.annotated.all.collapsed <- NULL
    } else {
        x <- input[[5]]
        if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]
        # aggregalo, non po esse che ogni gene ha una riga differente anche se viene dalla stessa SV
        x[, tag := paste(seqnames, start, end, sep = '.')]

        # x <- x[, .(
        #     seqnames = first(seqnames),
        #     start = first(start),
        #     end = first(end),
        #     boh1 = first(boh1),
        #     CNV.Type = first(CNV.Type),
        #     Copies = first(Copies),
        #     Freq = first(Freq),
        #     Freq.Error = first(Freq.Error),
        #     Shift = first(Shift),
        #     start.1 = first(start.1),
        #     end.1 = first(end.1),
        #     gene_biotype = first(gene_biotype),
        #     symbol = paste(symbol, collapse = ",")
        # ), by = .(tag)]


        #! gAD.DEL
        x.annotated.gad.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! gAD.CNV.DEL
        x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.cnv.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! DECIPHER ANNOTATION.DEL
        x.annotated.decipher.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.decipher.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! CLINVAR ANNOTATION.DEL
        x.annotated.clinvar.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.clinvar.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

        #! gAD.DUP
        x.annotated.gad.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! gAD.CNV.DUP
        x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.cnv.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! DECIPHER ANNOTATION.DUP
        x.annotated.decipher.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.decipher.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! CLINVAR ANNOTATION.DUP
        x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.clinvar.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

        x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.clinvar.multi.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.clinvar.multi.dup), use.names = TRUE, fill = TRUE)

        x.annotated.all[, length := end - start]
        ng100k.annotated.all.collapsed <- copy(x.annotated.all)
        ng100k.annotated.all.collapsed$tool <- "NG100K"
    }


    # ####! ------------------------------------------ NG10K
    # if(length(input[[6]]) == 1 && input[[6]] == ''){
    #     ng10k.annotated.all.collapsed <- NULL
    # } else {
    #     x <- input[[6]]
    #     if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]
    #     # aggregalo, non po esse che ogni gene ha una riga differente anche se viene dalla stessa SV
    #     x[, tag := paste(seqnames, start, end, sep = '.')]

    #     # x <- x[, .(
    #     #     seqnames = first(seqnames),
    #     #     start = first(start),
    #     #     end = first(end),
    #     #     boh1 = first(boh1),
    #     #     CNV.Type = first(CNV.Type),
    #     #     Copies = first(Copies),
    #     #     Freq = first(Freq),
    #     #     Freq.Error = first(Freq.Error),
    #     #     Shift = first(Shift),
    #     #     start.1 = first(start.1),
    #     #     end.1 = first(end.1),
    #     #     gene_biotype = first(gene_biotype),
    #     #     symbol = paste(symbol, collapse = ",")
    #     # ), by = .(tag)]


    #     #! gAD.DEL
    #     x.annotated.gad.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! gAD.CNV.DEL
    #     x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.cnv.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! DECIPHER ANNOTATION.DEL
    #     x.annotated.decipher.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.decipher.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! CLINVAR ANNOTATION.DEL
    #     x.annotated.clinvar.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.clinvar.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #     #! gAD.DUP
    #     x.annotated.gad.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! gAD.CNV.DUP
    #     x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.cnv.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! DECIPHER ANNOTATION.DUP
    #     x.annotated.decipher.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.decipher.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! CLINVAR ANNOTATION.DUP
    #     x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.clinvar.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #     x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.clinvar.multi.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.clinvar.multi.dup), use.names = TRUE, fill = TRUE)
    #     x.annotated.all[, length := end - start]

    #     ng10k.annotated.all.collapsed <- copy(x.annotated.all)
    #     ng10k.annotated.all.collapsed$tool <- "NG10K"

    # }


    #! FINE
    return(list(GASOLINE = gasoline.annotated.all.collapsed, INTERTRANS = x.intertrans.aggregated, CUTESV = cutesv.annotated.all.collapsed, SNIFFLES2 = sniffles.annotated.all.collapsed, NG100K = ng100k.annotated.all.collapsed))

} 

collapse.annotation.6 <- function(x, overlap.names, reciprocal.overlap.threshold = 0.8){

    def.col <- c(
        "gene.group.specific",
        "counts",
        "tot.score",
        "HAPLOINSUFFICIENCY",
        "TRIPLOSENSITIVITY",
        "ONLINE REPORT",
        "DepMap.CNS.Min",
        "DepMap.CNS.Q1",
        "DepMap.CNS.Median",
        "DepMap.CNS.Q3",
        "DepMap.CNS.Max",
        "DepMap.Pediatric.Min",
        "DepMap.Pediatric.Q1",
        "DepMap.Pediatric.Median",
        "DepMap.Pediatric.Q3",
        "DepMap.Pediatric.Max",
        "decipher.#population_cnv_id",
        "decipher.deletion_observations",
        "decipher.deletion_frequency",
        "decipher.duplication_observations",
        "decipher.duplication_frequency",
        "decipher.observations",
        "decipher.frequency",
        "decipher.type",
        "decipher.sample_size",
        "decipher.study",
        "decipher.reciprocal.overlap",
        "gAD.CNV.sv.CNV",
        "gAD.CNV.symbol.CNV",
        "gAD.CNV.SF",
        "gAD.CNV.SF_nfe",
        "gAD.CNV.reciprocal.overlap",
        "gnomAD.svtype",
        "gnomAD.AF",
        "gnomAD.FREQ_HOMREF",
        "gnomAD.FREQ_HET",
        "gnomAD.FREQ_HOMALT",
        "gnomAD.AF_nfe",
        "gnomAD.FREQ_HOMREF_nfe",
        "gnomAD.FREQ_HET_nfe",
        "gAD.reciprocal.overlap",
        "clinvar.dbVAR.SV",
        "clinvar.dbVAR.id.variant",
        "clinvar.dbVAR.clinical",
        "clinvar.dbVAR.id.clinical",
        "clinvar.reciprocal.overlap",
        "tool"
    )


    # GASOLINE
    to.aggregate <- c(
        "Chr", "start", "end", "length", "NumReadsClustered", "SVLabelMostProb", "CallLabel", "symbol", def.col
    )

    cur1 <- x[[1]]
    mat <- t(apply(cur1[, ..overlap.names],1, function(x) {x}))
    ind.db <- apply(mat, 1, function(y) {
        all(y < 0.7 | is.na(y))
    })

    # qui trovo quelle senza overlap con qualche db
    cur1.not.overlapping <- cur1[ind.db, ]
    cur1.not.overlapping <- cur1.not.overlapping[, .SD[which.max(rowSums(.SD[, ..overlap.names], na.rm = TRUE))], by = .(tag, symbol)]

    # ora devo collassarle tra di loro --> trovo quelle che overlappano tra di loro
    tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.not.overlapping$Chr), start = gsub(", .*", "", cur1.not.overlapping$start), end = gsub(", .*", "", cur1.not.overlapping$end)))
    hits <- findOverlaps(tmp, tmp)
    gr1 <- tmp[queryHits(hits)]
    gr2 <- tmp[subjectHits(hits)]
    intersection <- pintersect(gr1, gr2)
    length.intersection <- width(intersection)
    width.gr1 <- width(gr1)
    width.gr2 <- width(gr2)
    reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

    # quelle che non overlappano
    hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
    hits.toadd.df <- data.frame(hits.toadd)
    hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
    ind.rm <- which(duplicated(hits.toadd.df))
    hits.toadd <- hits.toadd[-ind.rm]
    hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

    # quelle che overlappano (trovo anche le non collassabili, ovvero quelle che overlappavano solo tra loro stesse, che sarebbero dovute finire in hits.to.add)
    hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
    hit.counts1 <- table(queryHits(hits))
    hit.counts2 <- table(subjectHits(hits))
    isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
    cur1.not.overlapped.not.collasable <- cur1.not.overlapping[isolated.indices, ]

    # e quelle che overlappano con altre
    hits <- hits[which(queryHits(hits) != subjectHits(hits))]
    hits.df <- data.frame(hits)
    hits.df <- t(apply(hits.df, 1, sort))
    ind.rm <- which(duplicated(hits.df))
    hits <- hits[-ind.rm]

    hits.dt <- data.table(data.frame(hits))
    hits.query <- hits.dt$queryHits
    hits.query.universe <- unique(hits.query)
    ind <- sapply(hits.query.universe, function(x){
        subjectHits(hits)[which(hits.query %in% x)]
    }, simplify = FALSE)
    names(ind) <- hits.query.universe
    aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
        rbind(cur1.not.overlapping[query.x], cur1.not.overlapping[ind[[as.character(query.x)]]])
    })

    tags <- c()
    cur1.not.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
        out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
        if(!(out$tag[1] %in% tags)){
            tags <<- c(tags, out$tag) 
            return(out)
        } else { 
            return()
        }
    }))

    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))

    colm <- c("tag", to.aggregate)

    gasoline.not.overlapping.db.out <- rbind(
        cur1.not.overlapping.aggregated.dt,
        cur1.not.overlapping[hits.toadd.universe, ..colm],
        cur1.not.overlapped.not.collasable[,..colm]
    )

    # mentre quelli con l'overlap > 0.7 e lo stesso tag le collasso in singola riga
    cur1.overlapping <- cur1[-which(ind.db), ]
    if(nrow(cur1.overlapping) == 0){
        cur1.overlapping <- NULL
        cur1.not.not.collasable <- NULL
        cur1.overlapping.aggregated.dt <- NULL
    } else {
    # qui collasso le regioni tra loro stesse
        tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.overlapping$Chr), start = gsub(", .*", "", cur1.overlapping$start), end = gsub(", .*", "", cur1.overlapping$end)))
        hits <- findOverlaps(tmp, tmp)
        gr1 <- tmp[queryHits(hits)]
        gr2 <- tmp[subjectHits(hits)]
        intersection <- pintersect(gr1, gr2)
        length.intersection <- width(intersection)
        width.gr1 <- width(gr1)
        width.gr2 <- width(gr2)
        reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
        reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
        reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2, 2)

        hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
        hits.toadd.df <- data.frame(hits.toadd)
        hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
        ind.rm <- which(duplicated(hits.toadd.df))
        hits.toadd <- hits.toadd[-ind.rm]
        hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

        hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
        hit.counts1 <- table(queryHits(hits))
        hit.counts2 <- table(subjectHits(hits))
        isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
        cur1.not.not.collasable <- cur1.overlapping[isolated.indices, ]

        hits <- hits[which(queryHits(hits) != subjectHits(hits))]
        hits.df <- data.frame(hits)
        hits.df <- t(apply(hits.df, 1, sort))
        ind.rm <- which(duplicated(hits.df))
        hits <- hits[-ind.rm]

        hits.dt <- data.table(data.frame(hits))
        hits.query <- hits.dt$queryHits
        hits.query.universe <- unique(hits.query)
        ind <- sapply(hits.query.universe, function(x){
            subjectHits(hits)[which(hits.query %in% x)]
        }, simplify = FALSE)
        names(ind) <- hits.query.universe
        aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
            rbind(cur1.overlapping[query.x], cur1.overlapping[ind[[as.character(query.x)]]])
        })

        tags <- c()
        cur1.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
            out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
            if(!(out$tag[1] %in% tags)){
                tags <<- c(tags, out$tag) 
                return(out)
            } else { 
                return()
            }
        }))
    }

    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    # colm <- colnames(cur1.overlapping.aggregated.dt)
    colm <- c("tag", to.aggregate)

    list.of.dts <- list(
            cur1.overlapping.aggregated.dt,
            if (!is.null(cur1.overlapping) && length(hits.toadd.universe) > 0) cur1.overlapping[hits.toadd.universe, ..colm],
            if (!is.null(cur1.not.not.collasable)) cur1.not.not.collasable[,..colm],
            gasoline.not.overlapping.db.out
        )
    gasoline.out <- rbindlist(list.of.dts, use.names = TRUE, fill = TRUE)


    # INTERTRANS --> NO


    # CUTESV
    to.aggregate <- c(
        "seqnames", "start", "end", "length", "reads", "sv", "symbol", def.col    
     )

    cur1 <- x[[3]]
    mat <- t(apply(cur1[, ..overlap.names],1, function(x) {x}))
    ind.db <- apply(mat, 1, function(y) {
        all(y < 0.7 | is.na(y))
    })

    # qui trovo quelle senza overlap con qualche db
    cur1.not.overlapping <- cur1[ind.db, ]
    cur1.not.overlapping <- cur1.not.overlapping[, .SD[which.max(rowSums(.SD[, ..overlap.names], na.rm = TRUE))], by = .(tag, symbol)]

    # ora devo collassarle tra di loro --> trovo quelle che overlappano tra di loro
    tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.not.overlapping$seqnames), start = gsub(", .*", "", cur1.not.overlapping$start), end = gsub(", .*", "", cur1.not.overlapping$end)))
    hits <- findOverlaps(tmp, tmp)
    gr1 <- tmp[queryHits(hits)]
    gr2 <- tmp[subjectHits(hits)]
    intersection <- pintersect(gr1, gr2)
    length.intersection <- width(intersection)
    width.gr1 <- width(gr1)
    width.gr2 <- width(gr2)
    reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

    # quelle che non overlappano
    hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
    hits.toadd.df <- data.frame(hits.toadd)
    hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
    ind.rm <- which(duplicated(hits.toadd.df))
    hits.toadd <- hits.toadd[-ind.rm]
    hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

    # quelle che overlappano (trovo anche le non collassabili, ovvero quelle che overlappavano solo tra loro stesse, che sarebbero dovute finire in hits.to.add)
    hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
    hit.counts1 <- table(queryHits(hits))
    hit.counts2 <- table(subjectHits(hits))
    isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
    cur1.not.overlapped.not.collasable <- cur1.not.overlapping[isolated.indices, ]

    # e quelle che overlappano con altre
    hits <- hits[which(queryHits(hits) != subjectHits(hits))]
    hits.df <- data.frame(hits)
    hits.df <- t(apply(hits.df, 1, sort))
    ind.rm <- which(duplicated(hits.df))
    hits <- hits[-ind.rm]

    hits.dt <- data.table(data.frame(hits))
    hits.query <- hits.dt$queryHits
    hits.query.universe <- unique(hits.query)
    ind <- sapply(hits.query.universe, function(x){
        subjectHits(hits)[which(hits.query %in% x)]
    }, simplify = FALSE)
    names(ind) <- hits.query.universe
    aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
        rbind(cur1.not.overlapping[query.x], cur1.not.overlapping[ind[[as.character(query.x)]]])
    })

    tags <- c()
    cur1.not.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
        out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
        if(!(out$tag[1] %in% tags)){
            tags <<- c(tags, out$tag) 
            return(out)
        } else { 
            return()
        }
    }))

    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    colm <- c("tag", to.aggregate)
    cutesv.not.overlapping.db.out <- rbind(
        cur1.not.overlapping.aggregated.dt,
        cur1.not.overlapping[hits.toadd.universe, ..colm],
        cur1.not.overlapped.not.collasable[,..colm]
    )

    # mentre quelli con l'overlap > 0.7 e lo stesso tag le collasso in singola riga
    cur1.overlapping <- cur1[-which(ind.db), ]
    if(nrow(cur1.overlapping) == 0){
        cur1.overlapping <- NULL
        cur1.not.not.collasable <- NULL
        cur1.overlapping.aggregated.dt <- NULL
    } else {

        # qui collasso le regioni tra loro stesse
        tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.overlapping$seqnames), start = gsub(", .*", "", cur1.overlapping$start), end = gsub(", .*", "", cur1.overlapping$end)))
        hits <- findOverlaps(tmp, tmp)
        gr1 <- tmp[queryHits(hits)]
        gr2 <- tmp[subjectHits(hits)]
        intersection <- pintersect(gr1, gr2)
        length.intersection <- width(intersection)
        width.gr1 <- width(gr1)
        width.gr2 <- width(gr2)
        reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
        reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
        reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

        hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
        hits.toadd.df <- data.frame(hits.toadd)
        hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
        ind.rm <- which(duplicated(hits.toadd.df))
        hits.toadd <- hits.toadd[-ind.rm]
        hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

        hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
        hit.counts1 <- table(queryHits(hits))
        hit.counts2 <- table(subjectHits(hits))
        isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
        cur1.not.not.collasable <- cur1.overlapping[isolated.indices, ]

        hits <- hits[which(queryHits(hits) != subjectHits(hits))]
        hits.df <- data.frame(hits)
        hits.df <- t(apply(hits.df, 1, sort))
        ind.rm <- which(duplicated(hits.df))
        hits <- hits[-ind.rm]

        hits.dt <- data.table(data.frame(hits))
        hits.query <- hits.dt$queryHits
        hits.query.universe <- unique(hits.query)
        ind <- sapply(hits.query.universe, function(x){
            subjectHits(hits)[which(hits.query %in% x)]
        }, simplify = FALSE)
        names(ind) <- hits.query.universe
        aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
            rbind(cur1.overlapping[query.x], cur1.overlapping[ind[[as.character(query.x)]]])
        })

        tags <- c()
        cur1.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
            out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
            if(!(out$tag[1] %in% tags)){
                tags <<- c(tags, out$tag) 
                return(out)
            } else { 
                return()
            }
        }))
    }

    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    colm <- c("tag", to.aggregate)

    list.of.dts <- list(
            cur1.overlapping.aggregated.dt,
            if (!is.null(cur1.overlapping) && length(hits.toadd.universe) > 0) cur1.overlapping[hits.toadd.universe, ..colm],
            if (!is.null(cur1.not.not.collasable)) cur1.not.not.collasable[,..colm],
            cutesv.not.overlapping.db.out
        )
    cutesv.out <- rbindlist(list.of.dts, use.names = TRUE, fill = TRUE)

    # SNIFFLES2
    to.aggregate <- c(
        "seqnames", "start", "end", "length", "svtype", "SUPPORT", "COVERAGE", "VAF", "hap", "symbol", def.col
    )

    cur1 <- x[[4]]
    mat <- t(apply(cur1[, ..overlap.names],1, function(x) {x}))
    ind.db <- apply(mat, 1, function(y) {
        all(y < 0.7 | is.na(y))
    })

    # qui trovo quelle senza overlap con qualche db
    cur1.not.overlapping <- cur1[which(ind.db), ]
    cur1.not.overlapping <- cur1.not.overlapping[, .SD[which.max(rowSums(.SD[, ..overlap.names], na.rm = TRUE))], by = .(tag, symbol)]

    # ora devo collassarle tra di loro --> trovo quelle che overlappano tra di loro
    tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.not.overlapping$seqnames), start = gsub(", .*", "", cur1.not.overlapping$start), end = gsub(", .*", "", cur1.not.overlapping$end)))
    hits <- findOverlaps(tmp, tmp)
    gr1 <- tmp[queryHits(hits)]
    gr2 <- tmp[subjectHits(hits)]
    intersection <- pintersect(gr1, gr2)
    length.intersection <- width(intersection)
    width.gr1 <- width(gr1)
    width.gr2 <- width(gr2)
    reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

    # quelle che non overlappano
    hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
    hits.toadd.df <- data.frame(hits.toadd)
    hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
    ind.rm <- which(duplicated(hits.toadd.df))
    hits.toadd <- hits.toadd[-ind.rm]
    hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

    # quelle che overlappano (trovo anche le non collassabili, ovvero quelle che overlappavano solo tra loro stesse, che sarebbero dovute finire in hits.to.add)
    hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
    hit.counts1 <- table(queryHits(hits))
    hit.counts2 <- table(subjectHits(hits))
    isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
    cur1.not.overlapped.not.collasable <- cur1.not.overlapping[isolated.indices, ]

    # e quelle che overlappano con altre
    hits <- hits[which(queryHits(hits) != subjectHits(hits))]
    hits.df <- data.frame(hits)
    hits.df <- t(apply(hits.df, 1, sort))
    ind.rm <- which(duplicated(hits.df))
    hits <- hits[-ind.rm]

    hits.dt <- data.table(data.frame(hits))
    hits.query <- hits.dt$queryHits
    hits.query.universe <- unique(hits.query)
    ind <- sapply(hits.query.universe, function(x){
        subjectHits(hits)[which(hits.query %in% x)]
    }, simplify = FALSE)
    names(ind) <- hits.query.universe
    aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
        rbind(cur1.not.overlapping[query.x], cur1.not.overlapping[ind[[as.character(query.x)]]])
    })

    tags <- c()
    cur1.not.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
        out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
        if(!(out$tag[1] %in% tags)){
            tags <<- c(tags, out$tag) 
            return(out)
        } else { 
            return()
        }
    }))

    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    colm <- c("tag", to.aggregate)
    sniffles2.not.overlapping.db.out <- rbind(
        cur1.not.overlapping.aggregated.dt,
        cur1.not.overlapping[hits.toadd.universe, ..colm],
        cur1.not.overlapped.not.collasable[,..colm]
    )

    # mentre quelli con l'overlap > 0.7 e lo stesso tag le collasso in singola riga
    cur1.overlapping <- cur1[-which(ind.db), ]

    if(nrow(cur1.overlapping) == 0){
        cur1.overlapping <- NULL
        cur1.not.not.collasable <- NULL
        cur1.overlapping.aggregated.dt <- NULL
    } else {
        # qui collasso le regioni tra loro stesse
        tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.overlapping$seqnames), start = gsub(", .*", "", cur1.overlapping$start), end = gsub(", .*", "", cur1.overlapping$end)))
        hits <- findOverlaps(tmp, tmp)
        gr1 <- tmp[queryHits(hits)]
        gr2 <- tmp[subjectHits(hits)]
        intersection <- pintersect(gr1, gr2)
        length.intersection <- width(intersection)
        width.gr1 <- width(gr1)
        width.gr2 <- width(gr2)
        reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
        reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
        reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

        hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
        hits.toadd.df <- data.frame(hits.toadd)
        hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
        ind.rm <- which(duplicated(hits.toadd.df))
        hits.toadd <- hits.toadd[-ind.rm]
        hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

        hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
        hit.counts1 <- table(queryHits(hits))
        hit.counts2 <- table(subjectHits(hits))
        isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
        cur1.not.not.collasable <- cur1.overlapping[isolated.indices, ]

        hits <- hits[which(queryHits(hits) != subjectHits(hits))]
        hits.df <- data.frame(hits)
        hits.df <- t(apply(hits.df, 1, sort))
        ind.rm <- which(duplicated(hits.df))
        hits <- hits[-ind.rm]

        hits.dt <- data.table(data.frame(hits))
        hits.query <- hits.dt$queryHits
        hits.query.universe <- unique(hits.query)
        ind <- sapply(hits.query.universe, function(x){
            subjectHits(hits)[which(hits.query %in% x)]
        }, simplify = FALSE)
        names(ind) <- hits.query.universe
        aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
            rbind(cur1.overlapping[query.x], cur1.overlapping[ind[[as.character(query.x)]]])
        })

        tags <- c()
        cur1.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
            out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
            if(!(out$tag[1] %in% tags)){
                tags <<- c(tags, out$tag) 
                return(out)
            } else { 
                return()
            }
        }))
    }
    hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    colm <- c("tag", to.aggregate)

    list.of.dts <- list(
            cur1.overlapping.aggregated.dt,
            if (!is.null(cur1.overlapping) && length(hits.toadd.universe) > 0) cur1.overlapping[hits.toadd.universe, ..colm],
            if (!is.null(cur1.not.not.collasable)) cur1.not.not.collasable[,..colm],
            sniffles2.not.overlapping.db.out
        )
    sniffles2.out <- rbindlist(list.of.dts, use.names = TRUE, fill = TRUE)


    # NG100K
    to.aggregate <- c(
        "seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "Freq", "symbol", def.col
    )

    cur1 <- x[[5]]
    if(is.null(nrow(cur1))){
        ng100k.out <- NULL
    } else {
        mat <- t(apply(cur1[, ..overlap.names],1, function(x) {x}))
        ind.db <- apply(mat, 1, function(y) {
            all(y < 0.7 | is.na(y))
        })

        # qui trovo quelle senza overlap con qualche db
        cur1.not.overlapping <- cur1[ind.db, ]
        cur1.not.overlapping <- cur1.not.overlapping[, .SD[which.max(rowSums(.SD[, ..overlap.names], na.rm = TRUE))], by = .(tag, symbol)]

        # ora devo collassarle tra di loro --> trovo quelle che overlappano tra di loro
        tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.not.overlapping$seqnames), start = gsub(", .*", "", cur1.not.overlapping$start), end = gsub(", .*", "", cur1.not.overlapping$end)))
        hits <- findOverlaps(tmp, tmp)
        gr1 <- tmp[queryHits(hits)]
        gr2 <- tmp[subjectHits(hits)]
        intersection <- pintersect(gr1, gr2)
        length.intersection <- width(intersection)
        width.gr1 <- width(gr1)
        width.gr2 <- width(gr2)
        reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
        reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
        reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

        # quelle che non overlappano al di sopra della soglia del reciprocal overlap
        hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
        hits.toadd.df <- data.frame(hits.toadd)
        hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
        ind.rm <- which(duplicated(hits.toadd.df))
        hits.toadd <- hits.toadd[-ind.rm]
        hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

        # quelle che overlappano (trovo anche le non collassabili, ovvero quelle che overlappavano solo tra loro stesse, che sarebbero dovute finire in hits.to.add)
        hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
        hit.counts1 <- table(queryHits(hits))
        hit.counts2 <- table(subjectHits(hits))
        isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
        cur1.not.overlapped.not.collasable <- cur1.not.overlapping[isolated.indices, ]

        # e quelle che overlappano con altre
        hits <- hits[which(queryHits(hits) != subjectHits(hits))]
        hits.df <- data.frame(hits)
        hits.df <- t(apply(hits.df, 1, sort))
        ind.rm <- which(duplicated(hits.df))
        hits <- hits[-ind.rm]

        hits.dt <- data.table(data.frame(hits))
        hits.query <- hits.dt$queryHits
        hits.query.universe <- unique(hits.query)
        ind <- sapply(hits.query.universe, function(x){
            subjectHits(hits)[which(hits.query %in% x)]
        }, simplify = FALSE)
        names(ind) <- hits.query.universe
        aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
            rbind(cur1.not.overlapping[query.x], cur1.not.overlapping[ind[[as.character(query.x)]]])
        })

        tags <- c()
        cur1.not.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
            out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
            if(!(out$tag[1] %in% tags)){
                tags <<- c(tags, out$tag) 
                return(out)
            } else { 
                return()
            }
        }))

        hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
        colm <- to.aggregate
        ng100k.not.overlapping.db.out <- rbind(
            cur1.not.overlapping.aggregated.dt,
            cur1.not.overlapping[hits.toadd.universe, ..colm],
            cur1.not.overlapped.not.collasable[,..colm], fill =
        )

        # mentre quelli con l'overlap > 0.7 e lo stesso tag le collasso in singola riga
        cur1.overlapping <- cur1[-which(ind.db), ]
        if(nrow(cur1.overlapping) == 0){
            cur1.overlapping <- NULL
            cur1.not.not.collasable <- NULL
            cur1.overlapping.aggregated.dt <- NULL
        } else {
            hold.columns <- c("seqnames", "start", "end", "boh1", "CNV.Type", "Copies", "Freq", "Freq.Error", "Shift", "start.1", "end.1", "gene_biotype", "symbol") # RIMOSSO 'tag'
            paste.columns <- setdiff(names(cur1.overlapping), c("tag", hold.columns))
            cols.to.process <- c(hold.columns, paste.columns)
            cur1.overlapping <- cur1.overlapping[, {
                results.list <- lapply(cols.to.process, function(col.name) {
                    if (col.name %in% hold.columns) {
                    .SD[[col.name]][1]
                    } else {
                    smart.paste(.SD[[col.name]])
                    }
                })
                setNames(results.list, cols.to.process)
            }, by = tag]

            # qui collasso le regioni tra loro stesse
            tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.overlapping$seqnames), start = gsub(", .*", "", cur1.overlapping$start), end = gsub(", .*", "", cur1.overlapping$end)))
            hits <- findOverlaps(tmp, tmp)
            gr1 <- tmp[queryHits(hits)]
            gr2 <- tmp[subjectHits(hits)]
            intersection <- pintersect(gr1, gr2)
            length.intersection <- width(intersection)
            width.gr1 <- width(gr1)
            width.gr2 <- width(gr2)
            reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
            reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
            reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

            hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
            hits.toadd.df <- data.frame(hits.toadd)
            hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
            ind.rm <- which(duplicated(hits.toadd.df))
            hits.toadd <- hits.toadd[-ind.rm]
            hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

            hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
            hit.counts1 <- table(queryHits(hits))
            hit.counts2 <- table(subjectHits(hits))
            isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
            cur1.not.not.collasable <- cur1.overlapping[isolated.indices, ]

            hits <- hits[which(queryHits(hits) != subjectHits(hits))]
            hits.df <- data.frame(hits)
            hits.df <- t(apply(hits.df, 1, sort))
            ind.rm <- which(duplicated(hits.df))
            hits <- hits[-ind.rm]

            hits.dt <- data.table(data.frame(hits))
            hits.query <- hits.dt$queryHits
            hits.query.universe <- unique(hits.query)
            ind <- sapply(hits.query.universe, function(x){
                subjectHits(hits)[which(hits.query %in% x)]
            }, simplify = FALSE)
            names(ind) <- hits.query.universe
            aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
                rbind(cur1.overlapping[query.x], cur1.overlapping[ind[[as.character(query.x)]]])
            })

            tags <- c()
            cur1.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
                out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
                if(!(out$tag[1] %in% tags)){
                    tags <<- c(tags, out$tag) 
                    return(out)
                } else { 
                    return()
                }
            }))
        }

        hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
        colm <- to.aggregate

        list.of.dts <- list(
                cur1.overlapping.aggregated.dt,
                if (!is.null(cur1.overlapping) && length(hits.toadd.universe) > 0) cur1.overlapping[hits.toadd.universe, ..colm],
                if (!is.null(cur1.not.not.collasable)) cur1.not.not.collasable[,..colm],
                ng100k.not.overlapping.db.out
            )
        ng100k.out <- rbindlist(list.of.dts, use.names = TRUE, fill = TRUE)
    }


    # # NG10K

    # to.aggregate <- c(
    #     "seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "Freq", "symbol", def.col
    # )

    # cur1 <- x[[6]]
    # if(is.null(nrow(cur1))){
    #     ng10k.out <- NULL
    # } else {

    #     mat <- t(apply(cur1[, ..overlap.names],1, function(x) {x}))
    #     ind.db <- apply(mat, 1, function(y) {
    #         all(y < 0.7 | is.na(y))
    #     })

    #     # qui trovo quelle senza overlap con qualche db
    #     cur1.not.overlapping <- cur1[ind.db, ]
    #     cur1.not.overlapping <- cur1.not.overlapping[, .SD[which.max(rowSums(.SD[, ..overlap.names], na.rm = TRUE))], by = .(tag, symbol)]

    #     # ora devo collassarle tra di loro --> trovo quelle che overlappano tra di loro
    #     tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.not.overlapping$seqnames), start = gsub(", .*", "", cur1.not.overlapping$start), end = gsub(", .*", "", cur1.not.overlapping$end)))
    #     hits <- findOverlaps(tmp, tmp)
    #     gr1 <- tmp[queryHits(hits)]
    #     gr2 <- tmp[subjectHits(hits)]
    #     intersection <- pintersect(gr1, gr2)
    #     length.intersection <- width(intersection)
    #     width.gr1 <- width(gr1)
    #     width.gr2 <- width(gr2)
    #     reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    #     reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    #     reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

    #     # quelle che non overlappano
    #     hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
    #     hits.toadd.df <- data.frame(hits.toadd)
    #     hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
    #     ind.rm <- which(duplicated(hits.toadd.df))
    #     hits.toadd <- hits.toadd[-ind.rm]
    #     hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

    #     # quelle che overlappano (trovo anche le non collassabili, ovvero quelle che overlappavano solo tra loro stesse, che sarebbero dovute finire in hits.to.add)
    #     hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
    #     hit.counts1 <- table(queryHits(hits))
    #     hit.counts2 <- table(subjectHits(hits))
    #     isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
    #     cur1.not.overlapped.not.collasable <- cur1.not.overlapping[isolated.indices, ]

    #     # e quelle che overlappano con altre
    #     hits <- hits[which(queryHits(hits) != subjectHits(hits))]
    #     hits.df <- data.frame(hits)
    #     hits.df <- t(apply(hits.df, 1, sort))
    #     ind.rm <- which(duplicated(hits.df))
    #     hits <- hits[-ind.rm]

    #     hits.dt <- data.table(data.frame(hits))
    #     hits.query <- hits.dt$queryHits
    #     hits.query.universe <- unique(hits.query)
    #     ind <- sapply(hits.query.universe, function(x){
    #         subjectHits(hits)[which(hits.query %in% x)]
    #     }, simplify = FALSE)
    #     names(ind) <- hits.query.universe
    #     aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
    #         rbind(cur1.not.overlapping[query.x], cur1.not.overlapping[ind[[as.character(query.x)]]])
    #     })

    #     tags <- c()
    #     cur1.not.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
    #         out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
    #         if(!(out$tag[1] %in% tags)){
    #             tags <<- c(tags, out$tag) 
    #             return(out)
    #         } else { 
    #             return()
    #         }
    #     }))

    #     hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    #     colm <- c("tag", to.aggregate)
    #     ng10k.not.overlapping.db.out <- rbind(
    #         cur1.not.overlapping.aggregated.dt,
    #         cur1.not.overlapping[hits.toadd.universe, ..colm],
    #         cur1.not.overlapped.not.collasable[,..colm], fill = TRUE
    #     )

    #     # mentre quelli con l'overlap > 0.7 e lo stesso tag le collasso in singola riga
    #     cur1.overlapping <- cur1[-which(ind.db), ]
    #     if(nrow(cur1.overlapping) == 0){
    #         cur1.overlapping <- NULL
    #         cur1.not.not.collasable <- NULL
    #         cur1.overlapping.aggregated.dt <- NULL
    #     } else {
    #         hold.columns <- c("seqnames", "start", "end", "boh1", "CNV.Type", "Copies", "Freq", "Freq.Error", "Shift", "start.1", "end.1", "gene_biotype", "symbol") # RIMOSSO 'tag'
    #         paste.columns <- setdiff(names(cur1.overlapping), c("tag", hold.columns))
    #         cols.to.process <- c(hold.columns, paste.columns)
    #         cur1.overlapping <- cur1.overlapping[, {
    #             results.list <- lapply(cols.to.process, function(col.name) {
    #                 if (col.name %in% hold.columns) {
    #                 .SD[[col.name]][1]
    #                 } else {
    #                 smart.paste(.SD[[col.name]])
    #                 }
    #             })
    #             setNames(results.list, cols.to.process)
    #         }, by = tag]

    #         # qui collasso le regioni tra loro stesse
    #         tmp <- makeGRangesFromDataFrame(data.frame(seqnames = gsub(", .*", "", cur1.overlapping$seqnames), start = gsub(", .*", "", cur1.overlapping$start), end = gsub(", .*", "", cur1.overlapping$end)))
    #         hits <- findOverlaps(tmp, tmp)
    #         gr1 <- tmp[queryHits(hits)]
    #         gr2 <- tmp[subjectHits(hits)]
    #         intersection <- pintersect(gr1, gr2)
    #         length.intersection <- width(intersection)
    #         width.gr1 <- width(gr1)
    #         width.gr2 <- width(gr2)
    #         reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    #         reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    #         reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2,2)

    #         hits.toadd <- hits[reciprocal.overlap.mean < reciprocal.overlap.threshold]
    #         hits.toadd.df <- data.frame(hits.toadd)
    #         hits.toadd.df <- t(apply(hits.toadd.df, 1, sort))
    #         ind.rm <- which(duplicated(hits.toadd.df))
    #         hits.toadd <- hits.toadd[-ind.rm]
    #         hits.toadd.universe <- unique(c(queryHits(hits.toadd), subjectHits(hits.toadd))) # devo togliere quelle righe che invece overlappano > 0.7 (cioe qui ho magari anche la stessa riga che overlappa con una al 0.5 ma in hits con un altra a 0.7)

    #         hits <- hits[reciprocal.overlap.mean > reciprocal.overlap.threshold]
    #         hit.counts1 <- table(queryHits(hits))
    #         hit.counts2 <- table(subjectHits(hits))
    #         isolated.indices <- unique(as.numeric(names(hit.counts1[hit.counts1 == 1])), as.numeric(names(hit.counts2[hit.counts2 == 1])))
    #         cur1.not.not.collasable <- cur1.overlapping[isolated.indices, ]

    #         hits <- hits[which(queryHits(hits) != subjectHits(hits))]
    #         hits.df <- data.frame(hits)
    #         hits.df <- t(apply(hits.df, 1, sort))
    #         ind.rm <- which(duplicated(hits.df))
    #         hits <- hits[-ind.rm]

    #         hits.dt <- data.table(data.frame(hits))
    #         hits.query <- hits.dt$queryHits
    #         hits.query.universe <- unique(hits.query)
    #         ind <- sapply(hits.query.universe, function(x){
    #             subjectHits(hits)[which(hits.query %in% x)]
    #         }, simplify = FALSE)
    #         names(ind) <- hits.query.universe
    #         aggregated.overlapped.hits <- lapply(as.numeric(names(ind)), function(query.x){
    #             rbind(cur1.overlapping[query.x], cur1.overlapping[ind[[as.character(query.x)]]])
    #         })

    #         tags <- c()
    #         cur1.overlapping.aggregated.dt <- do.call(rbind, lapply(aggregated.overlapped.hits, function(aoh){
    #             out <- aoh[, lapply(.SD, paste, collapse = ", "), by = tag, .SDcols = to.aggregate]
    #             if(!(out$tag[1] %in% tags)){
    #                 tags <<- c(tags, out$tag) 
    #                 return(out)
    #             } else { 
    #                 return()
    #             }
    #         }))
    #     }

    #     hits.toadd.universe <- setdiff(hits.toadd.universe, unique(c(queryHits(hits), subjectHits(hits))))
    #     colm <- colnames(cur1.overlapping.aggregated.dt)

    #     list.of.dts <- list(
    #             cur1.overlapping.aggregated.dt,
    #             if (!is.null(cur1.overlapping) && length(hits.toadd.universe) > 0) cur1.overlapping[hits.toadd.universe, ..colm],
    #             if (!is.null(cur1.not.not.collasable)) cur1.not.not.collasable[,..colm],
    #             ng10k.not.overlapping.db.out
    #         )
    #     ng10k.out <- rbindlist(list.of.dts, use.names = TRUE, fill = TRUE)

    # }

    if(nrow(data.frame(gasoline.out)) > 0){
        colnames(gasoline.out)[2] <- "seqnames"
        # colnames.gasoline.out <- c("tag", "seqnames", "start", "end", "length", "SVLabelMostProb", "NumReadsClustered", "CallLabel", "symbol", grep("decipher", colnames(gasoline.out), value = T), grep("gAD", colnames(gasoline.out), value = T), grep("gnomAD", colnames(gasoline.out), value = T), grep("clinvar", colnames(gasoline.out), value = T), "tool")
        colnames.gasoline.out <- c("tag", "seqnames", "start", "end", "length", "SVLabelMostProb", "NumReadsClustered", "CallLabel", "symbol", def.col)
        gasoline.out.formatted <- gasoline.out[, ..colnames.gasoline.out]
        colnames(gasoline.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "hap", "symbol", def.col)
    } else {
        gasoline.out.formatted <- NULL
    }

    if(nrow(data.frame(cutesv.out)) > 0){
        colnames.cutesv.out <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
        cutesv.out.formatted <- cutesv.out[, ..colnames.cutesv.out]
        colnames(cutesv.out.formatted) <- c("tag", "seqnames", "start", "end", "length", "sv", "reads", "symbol", def.col)
    } else {
        cutesv.out.formatted <- NULL
    }

    if(nrow(data.frame(sniffles2.out)) > 0){
        colnames.sniffles2.out <- c("tag", "seqnames", "start", "end", "length", "svtype", "SUPPORT", "COVERAGE", "VAF", "hap", "symbol", def.col)
        sniffles2.out.formatted <- sniffles2.out[, ..colnames.sniffles2.out]
        colnames(sniffles2.out.formatted) <-  c("tag", "seqnames", "start", "end", "length", "sv", "reads", "coverage", "VAF", "hap", "symbol", def.col)
    } else {
        sniffles2.out.formatted <- NULL
    }

    if(nrow(data.frame(ng100k.out)) > 0){
        colnames.ng100k.out <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
        ng100k.out.formatted <- ng100k.out[, ..colnames.ng100k.out]
        colnames(ng100k.out.formatted) <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    } else {
        ng100k.out.formatted <- NULL
    }

    # if(nrow(data.frame(ng10k.out)) > 0){
    #     colnames.ng10k.out <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    #     ng10k.out.formatted <- ng10k.out[, ..colnames.ng10k.out]
    #     colnames(ng10k.out.formatted) <- c("seqnames", "start", "end", "length", "CNV.Type", "Shift", "Copies", "symbol", def.col)
    # } else {
    #     ng10k.out.formatted <- NULL
    # }

    combined.dt <- rbindlist(list(gasoline.out.formatted , cutesv.out.formatted, sniffles2.out.formatted, ng100k.out.formatted), fill = TRUE)
    combined.dt

}

clean <- function(column) {
    sapply(column, function(cell.value) {
        # 1. NUOVO CONTROLLO: Gestisce celle che sono NULL o vettori di lunghezza zero.
        #    Questo è il controllo che risolve l'errore "argument is of length zero".
        if (is.null(cell.value) || length(cell.value) == 0) {
            return(NA_character_)
        }

        # 2. CONTROLLO MIGLIORATO PER NA: Gestisce il caso in cui una cella
        #    possa contenere uno o più NA (es. list(NA, NA)).
        if (all(is.na(cell.value))) {
            return(NA_character_)
        }
        
        # Rimuoviamo gli NA espliciti se la cella è un vettore misto (es. c("A", NA))
        cell.value <- cell.value[!is.na(cell.value)]
        
        # Se dopo aver tolto gli NA non resta nulla, trattiamolo come vuoto/NA
        if (length(cell.value) == 0) {
            return(NA_character_)
        }

        # 3. Convertiamo la cella in carattere per sicurezza
        cell.string <- as.character(cell.value)

        # 4. Dividiamo la stringa in un vettore di elementi
        #    Il paste iniziale gestisce il caso in cui cell.value sia un vettore (es. c("A", "B"))
        elements <- unlist(strsplit(paste(cell.string, collapse=", "), ",[ ]*"))

        # 5. Troviamo gli elementi unici e non vuoti
        unique.elements <- unique(elements[elements != ""])

        # 6. Riuniamo gli elementi unici in una singola stringa
        return(paste(unique.elements, collapse = ", "))

    }, USE.NAMES = FALSE)
}

rm.duplicated.info <- function(to.clean){
    to.clean$symbol <- sapply(to.clean$symbol, function(x){
        paste0(unique(unlist(strsplit(x, ', '))), collapse = ', ')
    })
    to.clean
}

collapse.ranges.tools <- function(to.collapse){
    # Parte 1: Trova tutti gli overlap e separa gli intervalli isolati
    hits <- findOverlaps(makeGRangesFromDataFrame(to.collapse), makeGRangesFromDataFrame(to.collapse))
    hits <- as.data.frame(hits)
    hits <- t(apply(hits, 1, sort))
    hits <- as.data.frame(unique(hits))
    hits <- subset(hits, V1 != V2)

    if(nrow(hits) == 0) {
        # ... (gestione del caso senza overlap, come prima)
        return(to.collapse)
    }

    all.overlapping.indices <- unique(c(hits$V1, hits$V2))
    not.collapsable <- to.collapse[-all.overlapping.indices,]
    collapsable <- to.collapse[all.overlapping.indices,]

    # Parte 2: Lavora sul subset 'collapsable'
    hits2 <- findOverlaps(makeGRangesFromDataFrame(collapsable), makeGRangesFromDataFrame(collapsable))
    hits2 <- as.data.frame(hits2)
    hits2 <- t(apply(hits2, 1, sort))
    hits2 <- as.data.frame(unique(hits2))
    hits2 <- subset(hits2, V1 != V2)

    intersection.width <- width(pintersect(makeGRangesFromDataFrame(collapsable[hits2$V1,]), makeGRangesFromDataFrame(collapsable[hits2$V2,])))
    perc.from <- intersection.width / width(makeGRangesFromDataFrame(collapsable[hits2$V1,]))
    perc.to   <- intersection.width / width(makeGRangesFromDataFrame(collapsable[hits2$V2,]))
    
    is.strong.hit <- perc.from > 0.7 & perc.to > 0.7
    valid.edges <- hits2[is.strong.hit, ]
    
    # --- Creazione dell'oggetto intermedio per ispezione ---
    # Questo è l'oggetto che volevi, contenente le *coppie* di overlap deboli.
    # Viene creato ma NON usato per l'output finale, per evitare duplicati.
    hit.not.collapsable <- hits2[!is.strong.hit, ]
    
    # In questo punto dello script, puoi ispezionare l'oggetto come preferisci:
    # print("Controllo gli hit non collassabili:")
    # print(head(hit.not.collapsable))
    # Oppure usa browser() per mettere in pausa l'esecuzione qui:
    # browser()
    # --------------------------------------------------------

    # Parte 3: Raggruppa e collassa usando la logica corretta (solo 'valid.edges')
    graph.obj <- graph_from_data_frame(valid.edges, directed = FALSE, vertices = 1:nrow(collapsable))
    collapsable$gruppo.collasso <- components(graph.obj)$membership
    collapsable.final <- collapsable[, lapply(.SD, smart.paste), by = gruppo.collasso][, -1]

    # # Parte 4: Uniforma i tipi di 'not.collapsable'
    # if (nrow(not.collapsable) > 0) {
    #   colonne.lista <- names(which(sapply(not.collapsable, is.list)))
    #   if (length(colonne.lista) > 0) {
    #     not.collapsable[, (colonne.lista) := lapply(.SD, function(col) sapply(col, smart.paste)), .SDcols = colonne.lista]
    #   }
    # }

    # return(rbind(collapsable.final, not.collapsable, fill = TRUE))

    final.dt <- rbind(collapsable.final, not.collapsable, fill = TRUE)

    #    Questo garantisce un output pulito e consistente.
    colonne.lista.finali <- names(which(sapply(final.dt, is.list)))

    if (length(colonne.lista.finali) > 0) {
      for (col in colonne.lista.finali) {
        final.dt[[col]] <- sapply(final.dt[[col]], smart.paste)
      }
    }

    return(final.dt)
}

filter.by.reads <- function(x){
    x[[1]][x[[1]]$NumReadsClustered > 9, ]
    x[[2]][x[[2]]$Reads > 9, ]
    x[[3]][x[[3]]$reads > 9, ]
    x[[4]][x[[4]]$SUPPORT > 9, ]
    x[[5]][x[[5]]$Shift > 0.4 | x[[5]]$Shift < -0.8 ]
    # x[[6]][x[[6]]$Shift > 0.4 | x[[6]]$Shift < -0.8 ]
    x
}

smart.paste <- function(x) {
    tryCatch({
        unlisted.x <- unlist(x, recursive = TRUE)
        if (is.null(unlisted.x)) return(NA_character_)
        unlisted.x <- unlisted.x[!is.na(unlisted.x)]
        if (length(unlisted.x) == 0) return(NA_character_)
        return(paste(unique(unlisted.x), collapse = ", "))
    }, error = function(e) {
        return(paste("ERRORE DI UNLIST/PASTE:", e$message))
    })
}

get.max.or.min.from.string <- function(s, max.or.min) {
  
  if (is.na(s) || s == "") {
    return(NA_real_) # NA di tipo numerico
  }
  
  char.values <- strsplit(s, ", ")[[1]]
  num.values <- suppressWarnings(as.numeric(char.values))
  num.values <- num.values[!is.na(num.values)]
  if (length(num.values) == 0) {
    return(NA_real_)
  }

  if(max.or.min == 'max'){
    return(max(num.values))
  } else if (max.or.min == 'min') {
    return(min(num.values))
  } else {
    print("set max or min minchione")
  }
}

get.value.from.vector <- function(values, operation) {
  if (is.null(values) || length(values) == 0 || all(is.na(values))) {
    return(NA_real_)
  }
  num.values <- suppressWarnings(as.numeric(values))
  num.values <- num.values[!is.na(num.values)]
  if (length(num.values) == 0) {
    return(NA_real_)
  }
  if (operation == 'max') {
    return(max(num.values))
  } else if (operation == 'min') {
    return(min(num.values))
  } else {
    stop("L'argomento 'operation' deve essere 'max' o 'min'.")
  }
}

sort.chrom <- function(x){

    dt <- x
    dt[, min.start.for.sort := sapply(start, get.value.from.vector, operation = "min")]
    chromosome.levels <- mixedsort(unique(dt$seqnames))
    dt[, seqnames.for.sort := factor(seqnames, levels = chromosome.levels)]
    x.sorted <- dt[order(seqnames.for.sort, min.start.for.sort)]

    # Rimuovi le colonne temporanee dal nuovo data.table ordinato
    x.sorted[, c("min.start.for.sort", "seqnames.for.sort") := NULL]
    return(x.sorted)
}

check.cell <- function(x){
    check <- c(0,0)
    for(vert in 1:nrow(x)){
        for(oriz in 1:ncol(x)){
            cur <- nchar(x[[oriz]][vert])
            if(cur > 1000 & !is.na(cur)){
                check <- rbind(check, c(vert, oriz))
            }
        }
    }
    return(check)
}

replace.long.cells <- function(dt, replacement_value = "Too Long") {
  # Controlla che l'input sia un data.table
  if (!is.data.table(dt)) {
    stop("L'input deve essere un data.table")
  }
  
  # Trova le celle che necessitano di sostituzione
  ind_list <- check.cell(dt)
  
  # Se ci sono più di una riga negli indici, procedi con la sostituzione
  if(nrow(ind_list) > 1) {
    for(i in 2:nrow(ind_list)) {
      # Estrai gli indici
      row_idx <- ind_list[i, 1]
      col_idx <- ind_list[i, 2]
      
      # Assegna il valore usando la sintassi data.table
      set(dt, 
          i = row_idx, 
          j = col_idx, 
          value = replacement_value)
    }
  }
  
  # Restituisce il data.table modificato
  return(dt)
}

check.row.length <- function(x) {
    row_lengths <- numeric(nrow(x))
    
    for(vert in 1:nrow(x)) {
        total_chars <- 0
        for(oriz in 1:ncol(x)) {
            cur <- nchar(as.character(x[[oriz]][vert]))
            if(!is.na(cur)) {
                total_chars <- total_chars + cur
            }
        }
        row_lengths[vert] <- total_chars
    }
    
    return(row_lengths)
}

replace.long.row <- function(dt, replacement_value = "Too Long") {
  # Controlla che l'input sia un data.table
  if (!is.data.table(dt)) {
    stop("L'input deve essere un data.table")
  }
  
  # Trova le celle che necessitano di sostituzione
  ind_list <- check.row.length(dt)
  
  # Se ci sono più di una riga negli indici, procedi con la sostituzione
  if(nrow(ind_list) > 1) {
    for(i in 2:nrow(ind_list)) {
      # Estrai gli indici
      row_idx <- ind_list[i, 1]
      col_idx <- ind_list[i, 2]
      
      # Assegna il valore usando la sintassi data.table
      set(dt, 
          i = row_idx, 
          j = col_idx, 
          value = replacement_value)
    }
  }
  
  # Restituisce il data.table modificato
  return(dt)
}

flatten.data.table <- function(input.dt) {
  dt.flat <- copy(input.dt)
  cols.to.flatten <- names(which(sapply(dt.flat, is.list)))
  
  if (length(cols.to.flatten) > 0) {
    dt.flat[, (cols.to.flatten) := lapply(.SD, function(col) {
      unlist(lapply(col, paste, collapse = ", "))
    }), .SDcols = cols.to.flatten]
  }
  
  return(dt.flat)
}

exclude.from.universe <- function(df, universe, soglia.overlap = 0.90){

    hits <- findOverlaps(makeGRangesFromDataFrame(df), makeGRangesFromDataFrame(universe))
    intersection <- pintersect(makeGRangesFromDataFrame(df)[queryHits(hits)], makeGRangesFromDataFrame(universe)[subjectHits(hits)])
    length.intersection <- width(intersection)
    width.gr1 <- width(makeGRangesFromDataFrame(df)[queryHits(hits)])
    width.gr2 <- width(makeGRangesFromDataFrame(universe)[subjectHits(hits)])
    reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
    reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
    reciprocal.overlap.mean <- (reciprocal.overlap1 + reciprocal.overlap2) / 2
    hits.subset <- hits[which(reciprocal.overlap.mean <= soglia.overlap )]
    df.out <- df[queryHits(hits.subset),]

    unique.hits <- unique(queryHits(hits))
    to.check <- 1:nrow(df)
    not.hits <- setdiff(to.check, unique.hits)
    df.out <- rbind(df.out, df[not.hits,])

    df.out

}

auto.annotate <- function(x, sv.label, sv.type, db, db.tag, seqnames.col, start.col, end.col){
    sv.label.pos <- match(sv.label, colnames(x))

    x.sub <- x[x[[sv.label.pos]] == sv.type]

    if(sum(x[[sv.label.pos]] == sv.type) == 0) { 
        print(paste0("non ci sono ", sv.type," da annotare --> INTERROMPO"))
        out <- NULL
        return(out)
    } else {
        db.colnames <- colnames(db)
        templato <- db[0,]
        x.sub <- x[x[[sv.label.pos]] == sv.type]
        x.sub[, (db.colnames) := templato]         # cosi ogni colonna e della classe che mi serve, senno mi sbrocca
        x.sub[,(paste0(db.tag, ".reciprocal.overlap")) := NA_real_] # NA per valori numerici
        col.pos <- c(match(seqnames.col, colnames(x)), match(start.col, colnames(x)), match(end.col, colnames(x)))
        x.sub[, tag := do.call(paste, c(.SD, sep = ".")), .SDcols = col.pos]

        hits <- findOverlaps(makeGRangesFromDataFrame(x.sub), makeGRangesFromDataFrame(data.frame(db)))
        if(length(hits) > 0 ){
            gr1 <- makeGRangesFromDataFrame(x.sub[queryHits(hits),])
            gr2 <- makeGRangesFromDataFrame(db[subjectHits(hits)])
            intersection <- pintersect(gr1, gr2)
            length.intersection <- width(intersection)
            width.gr1 <- width(gr1)
            width.gr2 <- width(gr2)
            reciprocal.overlap1 <- ifelse(width.gr1 > 0, length.intersection / width.gr1, 0)
            reciprocal.overlap2 <- ifelse(width.gr2 > 0, length.intersection / width.gr2, 0)
            reciprocal.overlap.mean <- round((reciprocal.overlap1 + reciprocal.overlap2) / 2, 3)

            x.sub.to.add <- x.sub[setdiff(seq_len(nrow(x.sub)), queryHits(hits)),]
            x.sub.annotated <- x.sub[queryHits(hits),]
            x.sub.annotated[ , (paste0(db.tag, ".reciprocal.overlap")) := reciprocal.overlap.mean]
            x.sub.annotated[ , (db.colnames) := db[subjectHits(hits),]]
            x.sub.to.add[ , (paste0(db.tag, ".reciprocal.overlap")) := NA]
            x.sub.to.add[ , (db.colnames) := NA]
            out <- rbind(x.sub.annotated, x.sub.to.add)
        } else {
            print("non mi trova hits in gasoline del gAD")
            out <- x.sub # se non lo crea allora mettilo come x.del.gdv (x.del che potrebbe essere come x.del se anche prima non ha trovato nulla) (mi serve per gAD)
        }

        # qui aggrego (prendo solo la prima, sticazzi delle varie annoatzioni che non ci servono) cio che ha un reciprocal overlap basso, senno ho troppe righe
        out.to.add <- out[out[[(paste0(db.tag, ".reciprocal.overlap"))]] < 0.7]
        col.overlap <- paste0(db.tag, ".reciprocal.overlap")
        out.to.add <- out.to.add[, .SD[which.max(get(col.overlap))], by = .(tag, symbol)]
        out <- rbind(out[out[[(paste0(db.tag, ".reciprocal.overlap"))]] >= 0.7], out.to.add, out[is.na(out[[(paste0(db.tag, ".reciprocal.overlap"))]])])

        return(out)
    }
}

annotate.results6.list <- function(input){

    ####! ------------------------------------------ GASOLINE
    x <- input[[1]]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Del", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Del", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Del", db = decipher, db.tag = "decipher", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Del", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "Chr", start.col = "start", end.col = "end")


    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "TDup", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "TDup", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "TDup", db = decipher, db.tag = "decipher", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "TDup", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "Chr", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Inv", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Ins", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "Chr", start.col = "start", end.col = "end")
    #! CLINVAR.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "SVLabelMostProb", sv.type = "Ins", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "Chr", start.col = "start", end.col = "end")

    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)

    gasoline.annotated.all.collapsed <- copy(x.annotated.all)
    gasoline.annotated.all.collapsed$tool <- "Gasoline"


    ####! ------------------------------------------ INTERTRANS
    x <- input[[2]]
    x[, Chr := paste0(seqnames.1, "-", seqnames.2)]
    x[, start := paste0(start.1, "-", start.2)]
    x[, gene := paste0(gene1, "-", gene2)]
    x.intertrans.aggregated <- x
    x.intertrans.aggregated$tools <- "Intertrans"


    ####! ------------------------------------------ CUTESV

    x <- input[[3]]
    if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "sv", sv.type = "DEL", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DUP
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "sv", sv.type = "DUP", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "sv", sv.type = "INV", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "sv", sv.type = "INS", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "sv", sv.type = "INS", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)


    cutesv.annotated.all.collapsed <- copy(x.annotated.all)
    cutesv.annotated.all.collapsed$tool <- "CuteSV"


    ####! ------------------------------------------ SNIFFLES2

    x <- input[[4]]
    colnames(x)[14:18] <- paste0("ens.", colnames(x)[14:18])
    if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]

    #! gAD.DEL
    x.annotated.gad.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DEL
    x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DEL
    x.annotated.decipher.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DEL
    x.annotated.clinvar.del <- auto.annotate(x, sv.label = "svtype", sv.type = "DEL", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.DUP
    x.annotated.gad.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.CNV.DUP
    x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! DECIPHER ANNOTATION.DUP
    x.annotated.decipher.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.DUP
    x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "svtype", sv.type = "DUP", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #! gAD.INV
    x.annotated.gad.inv <- auto.annotate(x, sv.label = "svtype", sv.type = "INV", db = gnomad.sv.inv, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! gAD.INS
    x.annotated.gad.ins <- auto.annotate(x, sv.label = "svtype", sv.type = "INS", db = gnomad.sv.ins, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #! CLINVAR ANNOTATION.INS
    x.annotated.clinvar.ins <- auto.annotate(x, sv.label = "svtype", sv.type = "INS", db = clinvar.ins.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
   
    x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.gad.ins, x.annotated.clinvar.ins, x.annotated.gad.inv), use.names = TRUE, fill = TRUE)

    x.annotated.all[, SUPPORT := gsub(";.*", "", gsub(".*SUPPORT=", "", x.annotated.all$info))]
    x.annotated.all[, COVERAGE := gsub(";.*", "", gsub(".*COVERAGE=", "", x.annotated.all$info))]
    x.annotated.all[, VAF := gsub(";.*", "", gsub(".*AF=", "", x.annotated.all$info))]
    x.annotated.all[, hap := gsub(":.*", "", x.annotated.all$sample)]

    sniffles.annotated.all.collapsed <- copy(x.annotated.all)
    sniffles.annotated.all.collapsed$tool <- "Sniffles2"


    ####! ------------------------------------------ NG100K
    if(length(input[[5]]) == 1 && input[[5]] == ''){
        ng100k.annotated.all.collapsed <- NULL
    } else {
        x <- input[[5]]
        if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]
        # aggregalo, non po esse che ogni gene ha una riga differente anche se viene dalla stessa SV
        x[, tag := paste(seqnames, start, end, sep = '.')]

        x <- x[, .(
            seqnames = first(seqnames),
            start = first(start),
            end = first(end),
            boh1 = first(boh1),
            CNV.Type = first(CNV.Type),
            Copies = first(Copies),
            Freq = first(Freq),
            Freq.Error = first(Freq.Error),
            Shift = first(Shift),
            start.1 = first(start.1),
            end.1 = first(end.1),
            gene_biotype = first(gene_biotype),
            symbol = paste(symbol, collapse = ",")
        ), by = .(tag)]


        #! gAD.DEL
        x.annotated.gad.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! gAD.CNV.DEL
        x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.cnv.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! DECIPHER ANNOTATION.DEL
        x.annotated.decipher.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.decipher.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! CLINVAR ANNOTATION.DEL
        x.annotated.clinvar.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.clinvar.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

        #! gAD.DUP
        x.annotated.gad.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! gAD.CNV.DUP
        x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.gad.cnv.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! DECIPHER ANNOTATION.DUP
        x.annotated.decipher.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.decipher.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        #! CLINVAR ANNOTATION.DUP
        x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
        x.annotated.clinvar.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

        x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.clinvar.multi.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.clinvar.multi.dup), use.names = TRUE, fill = TRUE)

        x.annotated.all[, length := end - start]
        ng100k.annotated.all.collapsed <- copy(x.annotated.all)
        ng100k.annotated.all.collapsed$tool <- "NG100K"
    }


    # ####! ------------------------------------------ NG10K
    # if(length(input[[6]]) == 1 && input[[6]] == ''){
    #     ng10k.annotated.all.collapsed <- NULL
    # } else {
    #     x <- input[[6]]
    #     if(!grepl("chr", x$seqnames[1])) x[, seqnames := paste0("chr", seqnames)]
    #     # aggregalo, non po esse che ogni gene ha una riga differente anche se viene dalla stessa SV
    #     x[, tag := paste(seqnames, start, end, sep = '.')]

    #     x <- x[, .(
    #         seqnames = first(seqnames),
    #         start = first(start),
    #         end = first(end),
    #         boh1 = first(boh1),
    #         CNV.Type = first(CNV.Type),
    #         Copies = first(Copies),
    #         Freq = first(Freq),
    #         Freq.Error = first(Freq.Error),
    #         Shift = first(Shift),
    #         start.1 = first(start.1),
    #         end.1 = first(end.1),
    #         gene_biotype = first(gene_biotype),
    #         symbol = paste(symbol, collapse = ",")
    #     ), by = .(tag)]


    #     #! gAD.DEL
    #     x.annotated.gad.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.sv.del, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! gAD.CNV.DEL
    #     x.annotated.gad.cnv.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.cnv.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = gnomad.cnv.del, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! DECIPHER ANNOTATION.DEL
    #     x.annotated.decipher.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.decipher.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! CLINVAR ANNOTATION.DEL
    #     x.annotated.clinvar.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.clinvar.multi.del <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Double Deletion", db = clinvar.del.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #     #! gAD.DUP
    #     x.annotated.gad.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.sv.dup, db.tag = "gAD", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! gAD.CNV.DUP
    #     x.annotated.gad.cnv.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.gad.cnv.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = gnomad.cnv.dup, db.tag = "gAD.CNV", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! DECIPHER ANNOTATION.DUP
    #     x.annotated.decipher.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.decipher.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = decipher, db.tag = "decipher", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     #! CLINVAR ANNOTATION.DUP
    #     x.annotated.clinvar.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")
    #     x.annotated.clinvar.multi.dup <- auto.annotate(x, sv.label = "CNV.Type", sv.type = "Multi-Copy Duplication", db = clinvar.dup.formatted, db.tag = "clinvar", seqnames.col = "seqnames", start.col = "start", end.col = "end")

    #     x.annotated.all <- rbindlist(list(x.annotated.decipher.del, x.annotated.gad.cnv.del, x.annotated.gad.del, x.annotated.clinvar.del, x.annotated.clinvar.multi.del, x.annotated.gad.dup, x.annotated.gad.cnv.dup, x.annotated.decipher.dup, x.annotated.clinvar.dup, x.annotated.clinvar.multi.dup), use.names = TRUE, fill = TRUE)
    #     x.annotated.all[, length := end - start]

    #     ng10k.annotated.all.collapsed <- copy(x.annotated.all)
    #     ng10k.annotated.all.collapsed$tool <- "NG10K"

    # }


    #! FINE
    return(list(GASOLINE = gasoline.annotated.all.collapsed, INTERTRANS = x.intertrans.aggregated, CUTESV = cutesv.annotated.all.collapsed, SNIFFLES2 = sniffles.annotated.all.collapsed, NG100K = ng100k.annotated.all.collapsed))

} 

annotate_encode_blacklist <- function(dt, blacklist_path){
    
    blacklist_dt <- fread(blacklist_path, col.names = c("seqnames", "start", "end", "type"))
    blacklist_gr <- makeGRangesFromDataFrame(blacklist_dt, seqnames.field = "seqnames", start.field = "start", end.field = "end")
    
    columns.match <- colnames(dt)[which(colnames(dt) %in% c("seqnames", "start", "end"))]
    if(length(columns.match) < 3) {
        print(columns.match); print("skipping one annotation...maybe intertrans?")
    } else {

        cnv_gr <- makeGRangesFromDataFrame(dt, seqnames.field = columns.match[1], start.field = columns.match[2], end.field = columns.match[3], keep.extra.columns = FALSE)
        
        hits <- findOverlaps(cnv_gr, blacklist_gr)
        
        # Calcola overlap per ogni match
        intersection <- pintersect(cnv_gr[queryHits(hits)], blacklist_gr[subjectHits(hits)])
        
        # Crea data.table con overlap per ogni hit
        overlap_dt <- data.table(cnv_idx = queryHits(hits),overlap_width = width(intersection))
        
        # Somma overlap per CNV (gestisce overlap multipli)
        overlap_summary <- overlap_dt[, .(total_overlap = sum(overlap_width)), by = cnv_idx]
        
        # Calcola percentuale
        overlap_summary[, cnv_width := width(cnv_gr[cnv_idx])]
        overlap_summary[, overlap_pct := (total_overlap / cnv_width) * 100]
        
        # Categorizza
        overlap_summary[, encode_category := fcase(
            overlap_pct < 25, 1L,
            overlap_pct < 50, 2L,
            overlap_pct < 75, 3L,
            overlap_pct >= 75, 4L
        )]
        
        # Crea risultato per TUTTE le CNV
        dt_result <- copy(dt)
        dt_result[, ENCODE.blacklist.overlap.pct := 0]
        dt_result[, ENCODE.blacklist.category := 1L]
        
        # Assegna valori per CNV con overlap
        dt_result[overlap_summary$cnv_idx, ENCODE.blacklist.overlap.pct := round(overlap_summary$overlap_pct, 2)]
        dt_result[overlap_summary$cnv_idx, ENCODE.blacklist.category := overlap_summary$encode_category]
        
        return(dt_result)
    }
}


