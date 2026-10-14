tf <- data.table::fread(paste0(dir, "/trrust_rawdata.human.tsv"))
onkoDB <- data.table::fread(paste0(dir, "/OncoKB.tsv"))
macrocateg <- read.table(paste0(dir, "/h.all.v2024.1.Hs.symbols.CON.MACROCATEGORIE.liberzon.et.al.gmt"), fill = T, header = T, sep = '\t')

gdv <- data.table::fread(paste0(dir, "/DGV.GS.hg38.gff3"), sep = '\t')
gdv.revamped <- gdv[, {

  # Estrazione di variant_sub_type (invariato)
  variant.sub.type <- sub(".*variant_sub_type=([^;]+).*", "\\1", V9)
  
  # Estrazione di PopulationSummary (invariato)
  pop.summary.string <- sub(".*PopulationSummary=([^;]+).*", "\\1", V9)
  
  # Ordine delle popolazioni (invariato)
  pop.order <- c("European", "African", "Asian", "Mexican", "Admixed", "MiddleEast", "Unknown")
  
  # --- LOGICA DI ESTRAZIONE CORRETTA ---
  pop.values <- lapply(pop.order, function(pop) {
    # 1. Controlla se la popolazione ESISTE nella stringa originale
    if (grepl(paste0("\\b", pop, " "), pop.summary.string)) {
      # 2. Se esiste, estrai il numero
      val <- sub(paste0(".*\\b", pop, " (\\d+)\\b.*"), "\\1", pop.summary.string)
      as.integer(val)
    } else {
      # 3. Se non esiste, il valore è 0
      0L # 0L è un modo efficiente per specificare un intero 0
    }
  })
  
  # Assegna i nomi alle colonne (invariato)
  names(pop.values) <- pop.order
  
  ID = gsub(";.*", "", gsub("ID=", "", V9))


  # Combina i risultati (invariato)
  c(
    .(chr = V1, 
      start = V4, 
      end = V5,
      type = fifelse(variant.sub.type == "Gain", "DUP", 
                     fifelse(variant.sub.type == "Loss", "DEL", NA_character_))
    ), 
    pop.values, ID = ID
  )
}, by = 1:nrow(gdv)]
gdv.revamped[, nrow := NULL]
colnames(gdv.revamped) <- paste0("GDV.", colnames(gdv.revamped))
gdv.revamped.del <- subset(gdv.revamped, gdv.revamped$GDV.type == "DEL")
gdv.revamped.dup <- subset(gdv.revamped, gdv.revamped$GDV.type == "DUP")

clinvar.del <- data.table::fread(paste0(dir, "/GRCh38.nr_deletions.pathogenic.tsv"))
clinvar.dup <- data.table::fread(paste0(dir, "/GRCh38.nr_duplications.somatic.tsv"))
clinvar.ins <- data.table::fread(paste0(dir, "/GRCh38.nr_insertions.somatic.tsv"))

if ("#chr" %in% names(clinvar.del)) {
    setnames(clinvar.del, "#chr", "chr")
}

clinvar.del.formatted <- clinvar.del[, .(
    seqnames = paste0("chr", chr),
    start = outermost_start,
    end = outermost_stop, # <-- Uso 'outermost_stop'
    dbVAR.SV = "DEL",
    dbVAR.id.variant = variant,
    dbVAR.clinical = clinical_assertion,
    dbVAR.id.clinical = clinvar_accession
)]

if ("#chr" %in% names(clinvar.dup)) {
    setnames(clinvar.dup, "#chr", "chr")
}

clinvar.dup.formatted <- clinvar.dup[, .(
    seqnames = paste0("chr", chr),
    start = outermost_start,
    end = outermost_stop, # <-- Uso 'outermost_stop' 
    dbVAR.SV = "dup",
    dbVAR.id.variant = variant,
    dbVAR.clinical = clinical_assertion,
    dbVAR.id.clinical = clinvar_accession
)]

if ("#chr" %in% names(clinvar.ins)) {
    setnames(clinvar.ins, "#chr", "chr")
}

clinvar.ins.formatted <- clinvar.ins[, .(
    seqnames = paste0("chr", chr),
    start = outermost_start,
    end = outermost_stop, # <-- Uso 'outermost_stop'
    dbVAR.SV = "ins",
    dbVAR.id.variant = variant,
    dbVAR.clinical = clinical_assertion,
    dbVAR.id.clinical = clinvar_accession
)]


decipher <- data.table::fread(paste0(dir, "/population_cnv_grch38.txt"))
decipher[, chr := paste0("chr", chr)]
ind <- c(6,7,9,10,12,13)
decipher[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]

colnames(clinvar.del.formatted) <- paste0("clinvar.", colnames(clinvar.del.formatted))
colnames(clinvar.dup.formatted) <- paste0("clinvar.", colnames(clinvar.dup.formatted))
colnames(clinvar.ins.formatted) <- paste0("clinvar.", colnames(clinvar.ins.formatted))
colnames(decipher) <- paste0("decipher.", colnames(decipher))


gnomad.cnv.source <- data.table::fread(cmd = paste("grep -v '^#'", paste0(dir, "/gnomad.v4.1.cnv.all.vcf")), sep = '\t')

column8 = lapply(gnomad.cnv.source$V8, function(x) {unlist(strsplit(x, ";"))})

gnomad.cnv <- data.table(
    seqnames = gnomad.cnv.source$V1, 
    start = gnomad.cnv.source$V2, 
    end = lapply(column8, function(x) {gsub("END=", "", x[1])}),
    sv.gnomAD.CNV = lapply(column8, function(x) {gsub("SVTYPE=", "", x[3])}),
    symbol.gnomAD.CNV = lapply(column8, function(x) {gsub("Genes=", "", x[8])}),
    SF = lapply(column8, function(x) {gsub("SF=", "", x[31])}),
    SF_afr = lapply(column8, function(x) {gsub("SF_afr=", "", x[32])}),
    SF_amr = lapply(column8, function(x) {gsub("SF_amr=", "", x[33])}),
    SF_asj = lapply(column8, function(x) {gsub("SF_asj=", "", x[34])}),
    SF_eas = lapply(column8, function(x) {gsub("SF_eas=", "", x[35])}),
    SF_fin = lapply(column8, function(x) {gsub("SF_fin=", "", x[36])}),
    SF_mid = lapply(column8, function(x) {gsub("SF_mid=", "", x[37])}),
    SF_nfe = lapply(column8, function(x) {gsub("SF_nfe=", "", x[38])}),
    SF_remaining = lapply(column8, function(x) {gsub("SF_remaining=", "", x[39])}),
    SF_sas = lapply(column8, function(x) {gsub("SF_sas=", "", x[40])})
)

gnomad.cnv.del <- subset(gnomad.cnv, gnomad.cnv$sv == "DEL")
colnames(gnomad.cnv.del) <- gsub("gnomAD.", "", colnames(gnomad.cnv.del))
colnames(gnomad.cnv.del) <- paste0("gAD.CNV.", colnames(gnomad.cnv.del))
gnomad.cnv.dup <- subset(gnomad.cnv, gnomad.cnv$sv == "DUP")
colnames(gnomad.cnv.dup) <- gsub("gnomAD.", "", colnames(gnomad.cnv.dup))
colnames(gnomad.cnv.dup) <- paste0("gAD.CNV.", colnames(gnomad.cnv.dup))


gnomad.sv <- fread(paste0(dir, "/gnomad.v4.1.sv.sites.selected.columns.bed"))
gnomad.sv <- subset(gnomad.sv, gnomad.sv$svtype != "BND")
colnames(gnomad.sv) <- paste0("gnomAD.", colnames(gnomad.sv))
gnomad.sv.del <- subset(gnomad.sv, gnomad.sv$gnomAD.svtype == "DEL")
gnomad.sv.dup <- subset(gnomad.sv, gnomad.sv$gnomAD.svtype == "DUP")
gnomad.sv.ins <- subset(gnomad.sv, gnomad.sv$gnomAD.svtype == "INS")
gnomad.sv.inv <- subset(gnomad.sv, gnomad.sv$gnomAD.svtype == "INV")

ind <- c(8,9,10,12,13,14)
gnomad.sv.del[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]
gnomad.sv.dup[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]
gnomad.sv.ins[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]
gnomad.sv.inv[, (ind) := lapply(.SD, function(z) round(as.numeric(z), 3)), .SDcols = ind]


# setDT(gdv.revamped.del)
# setDT(gdv.revamped.dup)
setDT(clinvar.del.formatted)
setDT(clinvar.dup.formatted)
setDT(clinvar.ins.formatted)
# gnomad.cnv.del
# gnomad.cnv.dup
# gnomad.bed.del
# gnomad.bed.dup
# gnomad.bed.cnv
# gnomad.bed.ins
# gnomad.bed.inv
# decipher
# gnomad.bed
string <- fread(paste0(dir, "/9606.FIXED.protein.links.v12.0.txt"))
voc <- data.table(biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, keys = string$protein1, keytype = "PROTEINID", columns = c("SYMBOL", "PROTEINID")))
voc <- na.omit(voc, cols = "PROTEINID")
setnames(string, 'protein1', "PROTEINID")
string <- merge(string, voc, all.x = TRUE, by = 'PROTEINID', allow.cartesian = TRUE)
setnames(string, "SYMBOL", "symbol")
string.abs <- copy(string)
string.abs.high <- subset(string.abs, combined_score >= 700)
string.abs.high[, counts := .N, by = symbol]
string.abs.high[, tot.score := sum(combined_score, na.rm = TRUE), by = symbol]
string.abs.medium <- subset(string.abs, combined_score >= 400)
string.abs.medium[, counts := .N, by = symbol]
string.abs.medium[, tot.score := sum(combined_score, na.rm = TRUE), by = symbol]

clingen <- fread(paste0(dir, "/Clingen-Dosage-Sensitivity-2025-08-30.csv"), skip = 4)
clingen <- clingen[-1]
setnames(clingen, "GENE SYMBOL", "symbol")
clingen <- clingen[,c(-2, -6)]

## --- SEZIONE DEPMAP: caricata SOLO per depmap="default".
##     `depmap.mode` viene assegnato da load_annotation_dbs (l'unico punto di
##     ingresso di questo blocco) nell'ambiente di valutazione, prima dell'eval.
##     Per none/custom i due CSV non sono necessari e non devono esistere
##     (es. immagine Docker senza CRISPRGeneEffect.csv/Model.csv). ---
if (depmap.mode == "default") {
depmap <- fread(paste0(dir, "/CRISPRGeneEffect.csv"))
colnames(depmap) <- gsub(" .*", "", colnames(depmap))
depmap.models <- fread(paste0(dir, "/Model.csv"))
depmap.models.cns <- subset(depmap.models, OncotreeLineage  == 'CNS/Brain')
depmap.models.cns.pediatric <- subset(depmap.models.cns, depmap.models.cns$AgeCategory == 'Pediatric')
depmap.cns <- subset(depmap, depmap$V1 %in% depmap.models.cns$ModelID)
depmap.cns.pediatric <- subset(depmap.cns, depmap.cns$V1 %in% depmap.models.cns.pediatric$ModelID)


depmap.models.mbl <- subset(depmap.models, depmap.models$DepmapModelType  == 'MBL')
depmap.mbl <- subset(depmap, depmap$V1 %in% depmap.models.mbl$ModelID)
colnames(depmap.mbl) <- gsub(" .*", "", colnames(depmap.mbl))

colnames(depmap.cns) <- gsub(" .*", "", colnames(depmap.cns))
colnames(depmap.cns.pediatric) <- gsub(" .*", "", colnames(depmap.cns.pediatric))
colnames(depmap.models.mbl) <- gsub(" .*", "", colnames(depmap.models.mbl))
}

