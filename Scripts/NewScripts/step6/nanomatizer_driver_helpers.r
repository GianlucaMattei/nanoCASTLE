# ! ------------------FIX INS (COME HO FATTO PRIMA VEDO LE INS CHE POSSONO AVERE UN EFFETTO SULLO SPLICEING, MA COSI AD OCCHIO... INVECE COSI VEDO SOLO LE INS CHE RICADONO SU ESONI SENZA INVENTRARMI NULLA)

exonsEnsDb <- GenomicFeatures::exonsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by = 'gene')

remove.ins.in.introns <- function(input, sv.tag, sv = "INS", chr = "seqnames"){
    input.ins <- input[get(sv.tag) == sv]
    if(nrow(input.ins)<1) return(input)
    new.chr <-   gsub("chr", "", input.ins[[chr]])
    # first.col <- names(input.ins)[1]
    input.ins[, (chr) := new.chr]
    input.ins[, end := start + 1]
    tmp.ins <- makeGRangesFromDataFrame(input.ins, seqnames.field = chr)
    hits3 <- findOverlaps(tmp.ins, unlist(exonsEnsDb))
    input.ins <- input.ins[queryHits(hits3),]
    input.ins[[chr]] <-   paste0("chr", input.ins[[chr]])
    input.ins <- unique(input.ins, by = c("start", "end", "symbol"))
    input.other <- input[get(sv.tag) != sv]

    return(rbind(input.other, input.ins))
}

build_sv_list <- function(name, env = parent.frame()) {
  # Carica data.table
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' non è installato.")
  }
  library(data.table)

  # Helper per recuperare gli oggetti dall'ambiente
  get_obj <- function(obj_name) {
    if (!exists(obj_name, envir = env, inherits = FALSE)) {
      stop("Oggetto non trovato: ", obj_name,
           ". Verifica che sia stato creato.")
    }
    get(obj_name, envir = env, inherits = FALSE)
  }

  ## -----------------------------
  ## gasoline.<name>.exons
  ## -----------------------------
  gasoline_exons_name <- paste0("gasoline.", name, ".exons")
  gasoline_exons <- get_obj(gasoline_exons_name)

  # colnames(gasoline.Pec.exons)[grep("SYMBOL", colnames(gasoline.Pec.exons))] <- 'symbol'
  colnames(gasoline_exons)[grep("SYMBOL", colnames(gasoline_exons))] <- "symbol"

  # gasoline.Pec.exons.dt <- data.table(gasoline.Pec.exons, key = "symbol")
  gasoline_exons_dt <- data.table(gasoline_exons, key = "symbol")

  # gasoline.Pec.exons.dt$tag <- apply(gasoline.Pec.exons.dt, 1, function(x){
  #    gsub(" +", "", paste0(x["seqnames"], ":", x["start"], "-", x["end"]))
  # }) 
  gasoline_exons_dt$tag <- apply(gasoline_exons_dt, 1, function(x) {
    gsub(" +", "", paste0(x["seqnames"], ":", x["start"], "-", x["end"]))
  })

  # gasoline.Pec.exons.dt <- gasoline.Pec.exons.dt[,-1]
  gasoline_exons_dt <- gasoline_exons_dt[, -1]

  # gasoline.Pec.exons.dt <- gasoline.Pec.exons.dt[!duplicated(gasoline.Pec.exons.dt),]
  gasoline_exons_dt <- gasoline_exons_dt[!duplicated(gasoline_exons_dt), ]

  ## -----------------------------
  ## gasoline.<name>.intertrans
  ## -----------------------------
  gasoline_intertrans_name <- paste0("gasoline.", name, ".intertrans")
  gasoline_intertrans <- get_obj(gasoline_intertrans_name)

  # gasoline.Pec.intertrans.dt <- data.table(gasoline.Pec.intertrans, key = c("gene1", "gene2"))
  gasoline_intertrans_dt <- data.table(gasoline_intertrans, key = c("gene1", "gene2"))

  # gasoline.Pec.intertrans.dt$tag <- apply(gasoline.Pec.intertrans.dt, 1, function(x){
  #    gsub(" +", "", paste0(x[1], "-", x[2], ":", x[3], "-", x[4]))
  # }) 
  gasoline_intertrans_dt$tag <- apply(gasoline_intertrans_dt, 1, function(x) {
    gsub(" +", "", paste0(x[1], "-", x[2], ":", x[3], "-", x[4]))
  })

  ## -----------------------------
  ## cutesv.<name>.exons
  ## -----------------------------
  cutesv_exons_name <- paste0("cutesv.", name, ".exons")
  cutesv_exons <- get_obj(cutesv_exons_name)

  # colnames(cutesv.Pec.exons)[grep("SYMBOL", colnames(cutesv.Pec.exons))] <- 'symbol'
  colnames(cutesv_exons)[grep("SYMBOL", colnames(cutesv_exons))] <- "symbol"

  # cutesv.Pec.exons.dt <- data.table(cutesv.Pec.exons, key = "symbol")
  cutesv_exons_dt <- data.table(cutesv_exons, key = "symbol")

  # cutesv.Pec.exons.dt$tag <- apply(cutesv.Pec.exons.dt, 1, function(x){
  #    gsub(" +", "", paste0("chr", x["seqnames"], ":", x["start"], "-", x["end"]))
  # }) 
  cutesv_exons_dt$tag <- apply(cutesv_exons_dt, 1, function(x) {
    gsub(" +", "", paste0("chr", x["seqnames"], ":", x["start"], "-", x["end"]))
  })

  ## -----------------------------
  ## sniffles2.<name>.exons
  ## -----------------------------
  sniffles_exons_name <- paste0("sniffles2.", name, ".exons")
  sniffles_exons <- get_obj(sniffles_exons_name)

  # colnames(sniffles2.Pec.exons)[grep("SYMBOL", colnames(sniffles2.Pec.exons))] <- 'symbol'
  colnames(sniffles_exons)[grep("SYMBOL", colnames(sniffles_exons))] <- "symbol"

  # sniffles2.Pec.exons.dt <- data.table(sniffles2.Pec.exons, key = "symbol")
  sniffles_exons_dt <- data.table(sniffles_exons, key = "symbol")

  # sniffles2.Pec.exons.dt$tag <- apply(sniffles2.Pec.exons.dt, 1, function(x){
  #    gsub(" +", "", paste0("chr", x["seqnames"], ":", x["start"], "-", x["end"]))
  # }) 
  sniffles_exons_dt$tag <- apply(sniffles_exons_dt, 1, function(x) {
    gsub(" +", "", paste0("chr", x["seqnames"], ":", x["start"], "-", x["end"]))
  })

  ## -----------------------------
  ## ng100k.<name>.exons
  ## -----------------------------
  ng100k_exons_name <- paste0("ng100k.", name, ".exons")
  ng100k_exons <- get_obj(ng100k_exons_name)

  # ng100k.Pec.exons.dt <- data.table(ng100k.Pec.exons, key = "symbol")
  ng100k_exons_dt <- data.table(ng100k_exons, key = "symbol")

  # ng100k.Pec.exons.dt$tag <- apply(ng100k.Pec.exons.dt, 1, function(x){
  #    gsub(" +", "", paste0("chr", x[1], ":", x[2], "-", x[3]))
  # }) 
  ng100k_exons_dt$tag <- apply(ng100k_exons_dt, 1, function(x) {
    gsub(" +", "", paste0("chr", x[1], ":", x[2], "-", x[3]))
  })

  ## -----------------------------
  ## Lista finale (come list.Pec)
  ## -----------------------------
  result_list <- list(
    gasoline   = remove.ins.in.introns(gasoline_exons_dt, "SV.Type", "Ins", chr = "seqnames"),
    intertrans = gasoline_intertrans_dt,
    cutesv     = remove.ins.in.introns(cutesv_exons_dt, "sv"),
    sniffles2  = remove.ins.in.introns(sniffles_exons_dt, "svtype"),
    NG100K     = ng100k_exons_dt
  )

  return(result_list)
}

create.helping.col <- function(input){
    cur <- input[,grep("DepMap\\.CNS", colnames(input)), with = FALSE]
    input[, CNS.Sum.Score := apply(cur, 1, function(x){sum(as.numeric(x), na.rm = TRUE)})]
    CNS.Sum.Score.match <- NA
    CNS.Sum.Score.match[input$CNS.Sum.Score > 0 & (grepl("DUP", input$sv) | grepl("TDup", input$sv) | grepl("uplicat", input$CNV.Type) )] <- 1
    CNS.Sum.Score.match[input$CNS.Sum.Score < 0 & (grepl("DUP", input$sv) | grepl("TDup", input$sv) | grepl("uplicat", input$CNV.Type) )] <- -1
    input$CNS.Sum.Score.match <- CNS.Sum.Score.match

    cur <- input[,grep("DepMap\\.Pediatric", colnames(input)), with = FALSE]
    input[, Pediatric.Sum.Score := apply(cur, 1, function(x){sum(as.numeric(x), na.rm = TRUE)})]
    Pediatric.Sum.Score.match <- NA
    Pediatric.Sum.Score.match[input$Pediatric.Sum.Score > 0 & (grepl("DUP", input$sv) | grepl("TDup", input$sv) | grepl("uplicat", input$CNV.Type) )] <- 1
    Pediatric.Sum.Score.match[input$Pediatric.Sum.Score < 0 & (grepl("DUP", input$sv) | grepl("TDup", input$sv) | grepl("uplicat", input$CNV.Type) )] <- -1
    input$Pediatric.Sum.Score.match <- Pediatric.Sum.Score.match

    input
}
