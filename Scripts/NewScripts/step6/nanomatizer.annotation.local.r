library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(IRanges)
library(EnsDb.Hsapiens.v86)
library(writexl)

# this file is needed -- do I need to upgrade it?
# /opt/nanomatizer/Nanomatizer/Databsese/EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3


# Ottieni CDS e info geni (#! tolgo i cazzo di UTR)
cdsEnsDb <- cdsBy(EnsDb.Hsapiens.v86, by = 'gene')
geneInfo <- genes(EnsDb.Hsapiens.v86)

# Unlist
cdsUnlisted <- unlist(cdsEnsDb, use.names = FALSE)

# Crea vettore gene_id ripetuto per ogni range
gene_ids <- rep(names(cdsEnsDb), lengths(cdsEnsDb))

# Match vettorizzato
m <- match(gene_ids, geneInfo$gene_id)

# Aggiungi metadata
mcols(cdsUnlisted)$gene_id <- geneInfo$gene_id[m]
mcols(cdsUnlisted)$symbol <- geneInfo$symbol[m]
mcols(cdsUnlisted)$gene_name <- geneInfo$gene_name[m]
mcols(cdsUnlisted)$gene_biotype <- geneInfo$gene_biotype[m]

# Relist se vuoi mantenere la struttura per gene
cdsEnsDb <- relist(cdsUnlisted, cdsEnsDb)

#! creating annotations
genesEnsDb <- GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
# exonsEnsDb <- GenomicFeatures::exonsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by = 'gene')
exonsEnsDb <- cdsEnsDb
regulatory <- # path derivato da P6 (definito dal driver che lo source); fallback container
if (exists("P6")) {
  # P6 = step6/  →  4 livelli su = T (…/Nanomatizer/Scripts/NewScripts/step6)
  T_ROOT <- normalizePath(file.path(P6, "..", "..", "..", ".."), winslash = "/")
} else {
  T_ROOT <- "/opt/nanomatizer"
}
regulatory <- read.table(file.path(T_ROOT, "Nanomatizer/Databsese/EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3"), header = F, sep = "\t", fill = T)
regulatory.coord <- GenomicRanges::makeGRangesFromDataFrame(data.frame(seqnames = paste0("chr", regulatory[,1]), start = regulatory[,4], end = regulatory[,5]))


#! CLAIRE
annotateClair <- function(x, sample.name.output, output.dir){
    x <- subset(x, x$V7 == "PASS")
    x <- subset(x, x$V6 >= 28)

    coords.df <- data.table(
        seqnames = gsub("chr", "", x[[1]]),
        start = pmin(x[[2]], x[[2]] + (nchar(gsub("\\.", "", x[[4]]))) - nchar(gsub("\\.", "", x[[3]]))),
        end = pmax(x[[2]], x[[2]] + (nchar(gsub("\\.", "", x[[4]]))) - nchar(gsub("\\.", "", x[[3]])))
    )

    fields <- strsplit(x$V9[1], ":")[[1]]  # Prendi i nomi dal primo row (assumendo siano uguali)
    values <- tstrsplit(x$V10, ":", fixed = TRUE)
    names(values) <- fields
    x[, (fields) := tstrsplit(V10, ":", fixed = TRUE)]
    x <- x[,c(-8,-9,-10)]

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(coords.df), unlist(exonsEnsDb))
    toadd <- data.frame(data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    colnames(toadd) <- paste0("exon.", colnames(toadd))
    sample.genes.exons <- cbind(x[S4Vectors::queryHits(hits3),], toadd)
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$exon.EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "exon.EnsID", by.y = "GENEID", all.x = TRUE)
    colnames(sample.genes.exons)[1:8] <- c("EnsgID", "Chr", "Pos", "Strand", "Ref", "Alt", "Qual", "Filter")
    sample.genes.exons <- cbind(sample.genes.exons[,c(-1, -4)], sample.genes.exons[,1])


    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/Claire/", sample.name.output, "/", sample.name.output, "germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/Claire/", sample.name.output, "/", sample.name.output, "germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    return(sample.genes.exons)
}


annotateClair2 <- function(x, sample.name.output, output.dir){

# -gene
# -NM_gene e proteina
# -Func.refGene: localizzazione (esone, introne splicing..)
# -exonicFunc.refGene: nonsynonymous SNV,   synonymous SNV,nonframeshift deletion, ...unknown
# -database: GnomaAD exome e genome...con all e popolazioni, se possibile COSMIC, ONCOKB, clinvar, altri database....OMIM_ID, OMIM_phenotype
# -gnomAD_nhomalt_total_homozygotes_exomes e genome
# -tool bioinformatici: CADD e reveal e per splicing spliceAI

    x <- subset(x, x$V7 == "PASS")
    x <- subset(x, x$V6 >= 28)

    coords.df <- data.table(
        seqnames = gsub("chr", "", x[[1]]),
        start = pmin(x[[2]], x[[2]] + (nchar(gsub("\\.", "", x[[4]]))) - nchar(gsub("\\.", "", x[[3]]))),
        end = pmax(x[[2]], x[[2]] + (nchar(gsub("\\.", "", x[[4]]))) - nchar(gsub("\\.", "", x[[3]])))
    )

    fields <- strsplit(x$V9[1], ":")[[1]]  # Prendi i nomi dal primo row (assumendo siano uguali)
    values <- tstrsplit(x$V10, ":", fixed = TRUE)
    names(values) <- fields
    x[, (fields) := tstrsplit(V10, ":", fixed = TRUE)]
    x <- x[,c(-8,-9,-10)]

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(coords.df), genesEnsDb)
    x.sub <- x[S4Vectors::queryHits(hits3),]

    coords.df <- data.table(
        seqnames = gsub("chr", "", x.sub[[1]]),
        start = pmin(x[[2]], x.sub[[2]] + (nchar(gsub("\\.", "", x.sub[[4]]))) - nchar(gsub("\\.", "", x.sub[[3]]))),
        end = pmax(x[[2]], x.sub[[2]] + (nchar(gsub("\\.", "", x.sub[[4]]))) - nchar(gsub("\\.", "", x.sub[[3]])))
    )
    


    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(coords.df), unlist(exonsEnsDb))
    toadd <- data.frame(data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    colnames(toadd) <- paste0("exon.", colnames(toadd))
    sample.genes.exons <- cbind(x[S4Vectors::queryHits(hits3),], toadd)
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$exon.EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "exon.EnsID", by.y = "GENEID", all.x = TRUE)
    colnames(sample.genes.exons)[1:8] <- c("EnsgID", "Chr", "Pos", "Strand", "Ref", "Alt", "Qual", "Filter")
    sample.genes.exons <- cbind(sample.genes.exons[,c(-1, -4)], sample.genes.exons[,1])


    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/Claire/", sample.name.output, "/", sample.name.output, "germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/Claire/", sample.name.output, "/", sample.name.output, "germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    return(sample.genes.exons)
}

#! CUTESV

annotateCuteSV <- function(sample.file, sample.name, output.dir, sample.name.output = NA){
    sample <- data.table::fread(cmd = paste0("grep -v '^##' ", sample.file))
    sample <- subset(sample, sample$FILTER == 'PASS')
    type <- gsub("\\..*", "", gsub(".*cuteSV\\." , "", sample$ID))
    sample.split <- split(sample, type)
    
    current.sv <- sample.split$DEL
    nchar.ref <- nchar(current.sv$REF)
    nchar.alt <- nchar(current.sv$ALT)
    sv.length <- abs(nchar.ref - nchar.alt)
    sv.start <- current.sv$POS
    sv.end <- sv.start + sv.length
    reads <-  gsub(":.*$", "", gsub(".\\/.:.:", "", data.frame(current.sv)[,10]))
    del.df <- data.frame(seqnames = current.sv[,1], start = sv.start, end = sv.end, sv = "DEL", reads = reads, length = sv.length)

    current.sv <- sample.split$INS 
    nchar.ref <- nchar(current.sv$REF)
    nchar.alt <- nchar(current.sv$ALT)
    sv.length <- abs(nchar.ref - nchar.alt)
    sv.start <- current.sv$POS
    sv.end <- sv.start + sv.length
    reads <-  gsub(":.*$", "", gsub(".\\/.:.:", "", data.frame(current.sv)[,10]))
    ins.df <- data.frame(seqnames = current.sv[,1], start = sv.start, end = sv.end, sv = "INS", reads = reads, length = sv.length)

    current.sv <- sample.split$DUP
    sv.end <- as.numeric(gsub(";.*", "", gsub(".*END=", "", current.sv$INFO)))
    sv.length <- as.numeric(gsub(";.*", "", gsub(".*SVLEN=", "", current.sv$INFO)))
    sv.start <- sv.end - sv.length
    reads <-  gsub(":.*$", "", gsub(".\\/.:.:", "", data.frame(current.sv)[,10]))
    dup.df <- data.frame(seqnames = current.sv[,1], start = sv.start, end = sv.end, sv = "DUP", reads = reads, length = sv.length)

    current.sv <- sample.split$INV
    sv.end <- as.numeric(gsub(";.*", "", gsub(".*END=", "", current.sv$INFO)))
    sv.length <- as.numeric(gsub(";.*", "", gsub(".*SVLEN=", "", current.sv$INFO)))
    sv.start <- sv.end - sv.length
    reads <-  gsub(":.*$", "", gsub(".\\/.:.:", "", data.frame(current.sv)[,10]))
    inv.df <- data.frame(seqnames = current.sv[,1], start = sv.start, end = sv.end, sv = "INV", reads = reads, length = sv.length)

    coords.df <- rbind(del.df, dup.df, ins.df, inv.df)
    colnames(coords.df)[1] <- "seqnames"

    coords.df <- subset(coords.df, coords.df$reads >= 7)

    #! Genes and Exons
    coords.df[,1] <- gsub("chr", '', coords.df[,1])
    hits2 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(coords.df), genesEnsDb)
    sample.genes <- data.frame(coords.df[S4Vectors::queryHits(hits2),], symbol = genesEnsDb[S4Vectors::subjectHits(hits2)]$symbol, ensID = genesEnsDb[S4Vectors::subjectHits(hits2)]$gene_id)
    sample.genes <- sample.genes[order(sample.genes$reads),]

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(coords.df), unlist(exonsEnsDb))
    sample.genes.exons <- cbind(coords.df[S4Vectors::queryHits(hits3),], data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "EnsID", by.y = "GENEID", all.x = TRUE)
    sample.genes.exons <- subset(sample.genes.exons, !duplicated(sample.genes.exons[c(1:5, 14)]))

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/CuteSV.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/CuteSV.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    writexl::write_xlsx(sample.genes, paste0(output.dir, "/CuteSV.", sample.name.output, ".germline.SVs.genes.xlsx"))
    write.csv(sample.genes, file = paste0(output.dir, "/CuteSV.", sample.name.output, ".germline.SVs.genes.csv"), row.names = FALSE, quote = FALSE)

}

#! GASOLINE

annotateGasoline <- function(sample.file, sample.name, output.dir, sample.name.output = NA, reads.thr = 7, sv.length.thr = 100000){

    sample.source <- data.frame(data.table::fread(sample.file))

    colnames(sample.source)[c(1, 4, 5)] <- c("seqnames", "start", "end")
    sample.source[,1] <- gsub("chr", "", sample.source[,1])
    colnames(sample.source) <- c("seqnames", "Chr", "NumberEventClustering", "start", "end", "length", "SDStart", "SDEnd", "Coesion", "NumReadsClustered", "AlnIns", "AlnDel", "SplitIns", "SplitDup", "SplitDel", "SplitInv", "ReadInv", "ReadIns", "ProbDel", "ProbIns", "ProbTDup", "ProbInv", "ProbCPIns", "ProbIDup", "SVLabelMostProb", "n", "k0", "LNorm1", "LNorm2", "LNorm3", "CallLabel")

    #! GERMINAL
    sample.subset <- subset(sample.source, sample.source$NumReadsClustered>reads.thr)
    sample.subset.tmp <- subset(sample.subset, sample.subset$length < sv.length.thr)
    sample.subset.long <- subset(sample.subset, sample.subset$length >= sv.length.thr)
    sample.subset <- sample.subset.tmp

    sample.subset.splitted <- split(sample.subset, sample.subset$SVLabelMostProb)
    sample.subset.splitted$Ins$end  <- sample.subset.splitted$Ins$start + 1

    sample.subset.for.genes <- do.call("rbind", sample.subset.splitted)

    cur.to.fix <- sample.subset.splitted$Inv
    cur.to.fix$start <- cur.to.fix$end + 1
    sample.subset.splitted$Inv$end <- sample.subset.splitted$Inv$start + 1
    sample.subset.splitted$Inv <- rbind(sample.subset.splitted$Inv, cur.to.fix)
    sample.subset <- do.call("rbind", sample.subset.splitted)
    sample.subset <- subset(sample.subset, sample.subset$NumReadsClustered > reads.thr)
    sample.subset.for.genes <- subset(sample.subset.for.genes, sample.subset.for.genes$NumReadsClustered > reads.thr)

    #! Genes and Exons
    hits2 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(sample.subset.for.genes[,c(1, 4, 5)]), genesEnsDb)
    sample.genes <- data.frame(sample.subset.for.genes[S4Vectors::queryHits(hits2),], symbol = genesEnsDb[S4Vectors::subjectHits(hits2)]$symbol, ensID = genesEnsDb[S4Vectors::subjectHits(hits2)]$gene_id)
    sample.genes <- sample.genes[order(sample.genes$NumReadsClustered),]

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(sample.subset[c(1,4,5)]), unlist(exonsEnsDb))
    sample.genes.exons <- cbind(sample.subset[S4Vectors::queryHits(hits3),], data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "EnsID", by.y = "GENEID", all.x = TRUE)
    sample.genes.exons <- subset(sample.genes.exons, !duplicated(sample.genes.exons[,c(-37:-33)])) # tolgo le colonne che non mi fanno riconoscere le righe duplicate (seqnames start ect fino a exon_id)

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    writexl::write_xlsx(sample.genes, paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.genes.xlsx"))
    write.csv(sample.genes, file = paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.genes.csv"), row.names = FALSE, quote = FALSE)

}



annotateGasoline2 <- function(sample.file, sample.name, output.dir, sample.name.output = NA, reads.thr = 7, sv.length.thr = 100000){

    sample.source <- data.frame(data.table::fread(sample.file))

    colnames(sample.source)[c(1, 2, 3)] <- c("seqnames", "start", "end")
    sample.source[,1] <- gsub("chr", "", sample.source[,1])
    # colnames(sample.source) <- c("seqnames", "NumberEventClustering", "start", "end", "length", "SDStart", "SDEnd", "Coesion", "NumReadsClustered", "AlnIns", "AlnDel", "SplitIns", "SplitDup", "SplitDel", "SplitInv", "ReadInv", "ReadIns", "ProbDel", "ProbIns", "ProbTDup", "ProbInv", "ProbCPIns", "ProbIDup", "SVLabelMostProb", "n", "k0", "LNorm1", "LNorm2", "LNorm3", "CallLabel")

    #! GERMINAL
    sample.subset <- subset(sample.source, sample.source$Clustered.Signatures>reads.thr)
    sample.subset.tmp <- subset(sample.subset, sample.subset$Size < sv.length.thr)
    sample.subset.long <- subset(sample.subset, sample.subset$Size >= sv.length.thr)
    sample.subset <- sample.subset.tmp

    sample.subset.splitted <- split(sample.subset, sample.subset$SV.Type)
    sample.subset.splitted$Ins$end  <- sample.subset.splitted$Ins$start + 1

    sample.subset.for.genes <- do.call("rbind", sample.subset.splitted)

    cur.to.fix <- sample.subset.splitted$Inv
    cur.to.fix$start <- cur.to.fix$end + 1
    sample.subset.splitted$Inv$end <- sample.subset.splitted$Inv$start + 1
    sample.subset.splitted$Inv <- rbind(sample.subset.splitted$Inv, cur.to.fix)
    sample.subset <- do.call("rbind", sample.subset.splitted)
    sample.subset <- subset(sample.subset, sample.subset$Clustered.Signatures > reads.thr)
    sample.subset.for.genes <- subset(sample.subset.for.genes, sample.subset.for.genes$Clustered.Signatures > reads.thr)

    #! Genes and Exons
    hits2 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(sample.subset.for.genes[,c(1, 2, 3)]), genesEnsDb)
    sample.genes <- data.frame(sample.subset.for.genes[S4Vectors::queryHits(hits2),], symbol = genesEnsDb[S4Vectors::subjectHits(hits2)]$symbol, ensID = genesEnsDb[S4Vectors::subjectHits(hits2)]$gene_id)
    sample.genes <- sample.genes[order(sample.genes$Clustered.Signatures),]

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(sample.subset[c(1,2,3)]), unlist(exonsEnsDb))
    sample.genes.exons <- cbind(sample.subset[S4Vectors::queryHits(hits3),], data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "EnsID", by.y = "GENEID", all.x = TRUE)
    sample.genes.exons <- sample.genes.exons[,c(-32:-28)]
    sample.genes.exons <- subset(sample.genes.exons, !duplicated(sample.genes.exons)) # tolgo le colonne che non mi fanno riconoscere le righe duplicate (seqnames start ect fino a exon_id)

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    writexl::write_xlsx(sample.genes, paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.genes.xlsx"))
    write.csv(sample.genes, file = paste0(output.dir, "/Gasoline.", sample.name.output, ".germline.SVs.genes.csv"), row.names = FALSE, quote = FALSE)

}


annotateInterTrans <- function(sample.file, sample.name, output.dir, sample.name.output = NA, reads.thr = 7){

    tabl <- read.table(
        sample.file,
        sep = "",        # uno o più spazi considerati come separatore
        fill = TRUE,     # riempi le colonne mancanti con NA
        na.strings = "", # eventuali stringhe vuote come NA
        col.names = paste0("V", 1:17),  # forza 17 colonne
        comment.char = "",
        stringsAsFactors = FALSE
    )
    colnames(tabl) <- c("seqnames.1", "seqnames.2", "V3", "start.1", "V5", "start.2", "V7", "V8", "V9", "V10", "V11", "Reads", "V13", "Prob.Homo.WT", "Prob.Ete", "Prom.Homo", "Prob.Gen")
    tabl.subset <- subset(tabl, tabl$Prob.Gen == '0/1' | tabl$Prob.Gen == '1/1')
    tabl.subset <- subset(tabl.subset, tabl.subset$Reads > reads.thr)

    tabl.subset$gene1 <- ''
    tabl.subset$regulatory1 <- ''
    tabl.subset$regulatory.type1 <- ''

    tabl.subset$gene2 <- ''
    tabl.subset$regulatory2 <- ''
    tabl.subset$regulatory.type2 <- ''

    #! gene 1
    tabl.subset.coord1 <- tabl.subset[,c(1,4,4)]
    tabl.subset.coord1[,3] <- as.numeric(tabl.subset.coord1[,3]) + 1
    tabl.subset.coord1[,1] <- gsub("chr", "", tabl.subset.coord1[,1])
    colnames(tabl.subset.coord1) <- c('seqnames', "start", "end")
    hits1.genes <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), genesEnsDb)
    if(length(hits1.genes) > 0){
        for(i in 1:length(hits1.genes)){
            tabl.subset[S4Vectors::queryHits(hits1.genes[i]), "gene1"] <- genesEnsDb[S4Vectors::subjectHits(hits1.genes[i])]$symbol
        }
    }

    hits.dist1 <- distanceToNearest(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), genesEnsDb, ignore.strand = TRUE)
    if(length(hits1.genes) > 0){
        for(i in 1:length(hits.dist1)){
            tabl.subset[S4Vectors::queryHits(hits.dist1[i]), "nearest.gene1"] <- genesEnsDb[S4Vectors::subjectHits(hits.dist1[i])]$symbol
            tabl.subset[S4Vectors::queryHits(hits.dist1[i]), "nearest.gene1.distance"] <- ifelse(is.null(data.frame(hits.dist1[i])[,'distance']), 0, data.frame(hits.dist1[i])[,'distance'])
        }
    }
    
    #! gene 2
    tabl.subset.coord2 <- tabl.subset[,c(2,6,6)]
    tabl.subset.coord2[,3] <- as.numeric(tabl.subset.coord2[,3]) + 1
    tabl.subset.coord2[,1] <- gsub("chr", "", tabl.subset.coord2[,1])
    colnames(tabl.subset.coord2) <- c('seqnames', "start", "end")
    hits2.genes <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), genesEnsDb)
    if(length(hits2.genes) > 0){
        for(i in 1:length(hits2.genes)){
            tabl.subset[S4Vectors::queryHits(hits2.genes[i]), "gene2"] <- genesEnsDb[S4Vectors::subjectHits(hits2.genes[i])]$symbol
        }
    }

    hits.dist2 <- distanceToNearest(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), genesEnsDb, ignore.strand = TRUE)
    if(length(hits2.genes) > 0){
        for(i in 1:length(hits2.genes)){
            tabl.subset[S4Vectors::queryHits(hits.dist2[i]), "nearest.gene2"] <- genesEnsDb[S4Vectors::subjectHits(hits.dist2[i])]$symbol
            tabl.subset[S4Vectors::queryHits(hits.dist2[i]), "nearest.gene2.distance"] <- ifelse(is.null(data.frame(hits.dist2[i])[,'distance']), 0, data.frame(hits.dist2[i])[,'distance'])

        }
    }



    #! reg 1
    tabl.subset.coord1 <- tabl.subset[,c(1,4,4)]
    tabl.subset.coord1[,3] <- as.numeric(tabl.subset.coord1[,3]) + 1
    tabl.subset.coord1[,1] <- gsub("chr", "", tabl.subset.coord1[,1])
    colnames(tabl.subset.coord1) <- c('seqnames', "start", "end")
    hits1.proms <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), GenomicRanges::makeGRangesFromDataFrame(regulatory.coord))
    if(length(hits1.proms)){
        for(i in 1:length(hits1.proms)){
            tabl.subset[S4Vectors::queryHits(hits1.proms[i]), "regulatory1"] <- regulatory[S4Vectors::subjectHits(hits1.proms[i]), "V9"]
            tabl.subset[S4Vectors::queryHits(hits1.proms[i]), "regulatory.type1"] <- regulatory[S4Vectors::subjectHits(hits1.proms[i]), "V3"]
        }
    }
    tabl.subset$regulatory1 <- gsub(";.*$", "", gsub("^.*gene_name=", "", tabl.subset$regulatory1))

    #! reg 2
    tabl.subset.coord2 <- tabl.subset[,c(2,6,6)]
    tabl.subset.coord2[,3] <- as.numeric(tabl.subset.coord2[,3]) + 1
    tabl.subset.coord2[,1] <- gsub("chr", "", tabl.subset.coord2[,1])
    colnames(tabl.subset.coord2) <- c('seqnames', "start", "end")
    hits2.proms <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), GenomicRanges::makeGRangesFromDataFrame(regulatory.coord))
    if(length(hits2.proms)){
        for(i in 1:length(hits2.proms)){
            tabl.subset[S4Vectors::queryHits(hits2.proms[i]), "regulatory2"] <- regulatory[S4Vectors::subjectHits(hits2.proms[i]), "V9"]
            tabl.subset[S4Vectors::queryHits(hits2.proms[i]), "regulatory.type2"] <- regulatory[S4Vectors::subjectHits(hits2.proms[i]), "V3"]
        }
    }
    tabl.subset$regulatory2 <- gsub(";.*$", "", gsub("^.*gene_name=", "", tabl.subset$regulatory2))

    tabl.subset <- tabl.subset[,c(-3,-5,-7,-8, -9, -10, -11 )]

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(tabl.subset, paste0(output.dir, "/GasolineInterTrans.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(tabl.subset, file = paste0(output.dir, "/GasolineInterTrans.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

}



annotateInterTrans2 <- function(sample.file, sample.name, output.dir, sample.name.output = NA, reads.thr = 7){

    tabl <- read.table(sample.file, header = T, sep = '\t', stringsAsFactors = FALSE, comment.char = "", na.strings = "")

    colnames(tabl) <- c("seqnames.1", "seqnames.2", "start.1", "start.2", "Break1.SD", "Break2.SD", "Coesion", "Clustered.Signatures", "Reads", "Alt.Reads", "Prob.Homo.WT", "Prob.Ete", "Prom.Homo", "Prob.Gen")
    tabl.subset <- subset(tabl, tabl$Prob.Gen == '0/1' | tabl$Prob.Gen == '1/1')
    tabl.subset <- subset(tabl.subset, tabl.subset$Reads > reads.thr)

    tabl.subset$gene1 <- ''
    tabl.subset$regulatory1 <- ''
    tabl.subset$regulatory.type1 <- ''

    tabl.subset$gene2 <- ''
    tabl.subset$regulatory2 <- ''
    tabl.subset$regulatory.type2 <- ''

    #! gene 1
    tabl.subset.coord1 <- tabl.subset[,c(1,3,3)]
    tabl.subset.coord1[,3] <- as.numeric(tabl.subset.coord1[,3]) + 1
    tabl.subset.coord1[,1] <- gsub("chr", "", tabl.subset.coord1[,1])
    colnames(tabl.subset.coord1) <- c('seqnames', "start", "end")
    hits1.genes <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), genesEnsDb)
    if(length(hits1.genes) > 0){
        for(i in 1:length(hits1.genes)){
            tabl.subset[S4Vectors::queryHits(hits1.genes[i]), "gene1"] <- genesEnsDb[S4Vectors::subjectHits(hits1.genes[i])]$symbol
        }
    }

    hits.dist1 <- distanceToNearest(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), genesEnsDb, ignore.strand = TRUE)
    if(length(hits1.genes) > 0){
        for(i in 1:length(hits.dist1)){
            tabl.subset[S4Vectors::queryHits(hits.dist1[i]), "nearest.gene1"] <- genesEnsDb[S4Vectors::subjectHits(hits.dist1[i])]$symbol
            tabl.subset[S4Vectors::queryHits(hits.dist1[i]), "nearest.gene1.distance"] <- ifelse(is.null(data.frame(hits.dist1[i])[,'distance']), 0, data.frame(hits.dist1[i])[,'distance'])
        }
    }
    
    #! gene 2
    tabl.subset.coord2 <- tabl.subset[,c(2,4,4)]
    tabl.subset.coord2[,3] <- as.numeric(tabl.subset.coord2[,3]) + 1
    tabl.subset.coord2[,1] <- gsub("chr", "", tabl.subset.coord2[,1])
    colnames(tabl.subset.coord2) <- c('seqnames', "start", "end")
    hits2.genes <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), genesEnsDb)
    if(length(hits2.genes) > 0){
        for(i in 1:length(hits2.genes)){
            tabl.subset[S4Vectors::queryHits(hits2.genes[i]), "gene2"] <- genesEnsDb[S4Vectors::subjectHits(hits2.genes[i])]$symbol
        }
    }

    hits.dist2 <- distanceToNearest(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), genesEnsDb, ignore.strand = TRUE)
    if(length(hits2.genes) > 0){
        for(i in 1:length(hits2.genes)){
            tabl.subset[S4Vectors::queryHits(hits.dist2[i]), "nearest.gene2"] <- genesEnsDb[S4Vectors::subjectHits(hits.dist2[i])]$symbol
            tabl.subset[S4Vectors::queryHits(hits.dist2[i]), "nearest.gene2.distance"] <- ifelse(is.null(data.frame(hits.dist2[i])[,'distance']), 0, data.frame(hits.dist2[i])[,'distance'])

        }
    }



    #! reg 1
    tabl.subset.coord1 <- tabl.subset[,c(1,3,3)]
    tabl.subset.coord1[,3] <- as.numeric(tabl.subset.coord1[,3]) + 1
    tabl.subset.coord1[,1] <- gsub("chr", "", tabl.subset.coord1[,1])
    colnames(tabl.subset.coord1) <- c('seqnames', "start", "end")
    hits1.proms <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord1), GenomicRanges::makeGRangesFromDataFrame(regulatory.coord))
    if(length(hits1.proms)){
        for(i in 1:length(hits1.proms)){
            tabl.subset[S4Vectors::queryHits(hits1.proms[i]), "regulatory1"] <- regulatory[S4Vectors::subjectHits(hits1.proms[i]), "V9"]
            tabl.subset[S4Vectors::queryHits(hits1.proms[i]), "regulatory.type1"] <- regulatory[S4Vectors::subjectHits(hits1.proms[i]), "V3"]
        }
    }
    tabl.subset$regulatory1 <- gsub(";.*$", "", gsub("^.*gene_name=", "", tabl.subset$regulatory1))

    #! reg 2
    tabl.subset.coord2 <- tabl.subset[,c(2,4,4)]
    tabl.subset.coord2[,3] <- as.numeric(tabl.subset.coord2[,3]) + 1
    tabl.subset.coord2[,1] <- gsub("chr", "", tabl.subset.coord2[,1])
    colnames(tabl.subset.coord2) <- c('seqnames', "start", "end")
    hits2.proms <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.coord2), GenomicRanges::makeGRangesFromDataFrame(regulatory.coord))
    if(length(hits2.proms)){
        for(i in 1:length(hits2.proms)){
            tabl.subset[S4Vectors::queryHits(hits2.proms[i]), "regulatory2"] <- regulatory[S4Vectors::subjectHits(hits2.proms[i]), "V9"]
            tabl.subset[S4Vectors::queryHits(hits2.proms[i]), "regulatory.type2"] <- regulatory[S4Vectors::subjectHits(hits2.proms[i]), "V3"]
        }
    }
    tabl.subset$regulatory2 <- gsub(";.*$", "", gsub("^.*gene_name=", "", tabl.subset$regulatory2))

    # tabl.subset <- tabl.subset[,c(-3,-5,-7,-8, -9, -10, -11 )]
    tabl.subset <- tabl.subset[, (c(5,6,7,8) * -1)]

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(tabl.subset, paste0(output.dir, "/GasolineInterTrans.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(tabl.subset, file = paste0(output.dir, "/GasolineInterTrans.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

}


#! SNIFFLES2

annotateSniffles <- function(sample.file, sample.name, output.dir, sample.name.output = NA, qual.thr = 20, sv.length.thr = 100000){
    tabl <- data.table::fread(cmd = paste0("grep -v '^##' ", sample.file))
    tabl <- data.frame(subset(tabl, tabl$QUAL > qual.thr))
    tabl <- data.frame(seqnames = tabl[,1], start = as.numeric(tabl[,2]), end = as.numeric(tabl[,2] + as.numeric(gsub(";.*", "", gsub(".*SVLEN=-?", "", tabl[,"INFO"])))), length = as.numeric(gsub(";.*", "", gsub(".*SVLEN=-?", "", tabl[,"INFO"]))), svtype = gsub("\\..*$", "", gsub("Sniffles2\\.", "", tabl[,"ID"])), filter = tabl[,"FILTER"], info = tabl[,"INFO"], format = tabl[,"FORMAT"], sample = tabl[,"SAMPLE"])
    tabl.subset <- subset(tabl, tabl$filter == "PASS")
    tabl.subset$end[is.na(tabl.subset$end)] <- tabl.subset$start[is.na(tabl.subset$end)] + 1
    tabl.subset[,1] <- gsub("chr", "", tabl.subset[,1])
    tabl.subset$info <- gsub(",", ";", tabl.subset$info)

    support <- as.numeric(gsub(";.*", "", gsub(".*SUPPORT=", "", tabl.subset[,7])))
    tabl.subset <- tabl.subset[support >=7, ]
    
    ind <- which(width(GenomicRanges::makeGRangesFromDataFrame(tabl.subset[,c(1, 2, 3)])) < sv.length.thr)
    tabl.subset <- tabl.subset[ind, ]
    tabl.subset.splitted <- split(tabl.subset, tabl.subset$svtype)

    tabl.subset.splitted$INS$end  <- tabl.subset.splitted$INS$start + 1
    tabl.subset.for.genes <- do.call("rbind", tabl.subset.splitted)

    cur.to.fix <- tabl.subset.splitted$INV
    cur.to.fix$start <- cur.to.fix$end + 1
    tabl.subset.splitted$INV$end <- tabl.subset.splitted$INV$start + 1
    tabl.subset.splitted$INV <- rbind(tabl.subset.splitted$INV, cur.to.fix)
    tabl.subset <- do.call("rbind", tabl.subset.splitted)

    hits2 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset.for.genes[,c(1, 2, 3)]), genesEnsDb)
    sample.genes <- data.frame(tabl.subset.for.genes[S4Vectors::queryHits(hits2),], symbol = genesEnsDb[S4Vectors::subjectHits(hits2)]$symbol, ensID = genesEnsDb[S4Vectors::subjectHits(hits2)]$gene_id)
    sample.genes$SUPPORT <- gsub(";.*", "", gsub(".*SUPPORT=", "", sample.genes$info))

    hits3 <- IRanges::findOverlaps(GenomicRanges::makeGRangesFromDataFrame(tabl.subset), unlist(exonsEnsDb))
    sample.genes.exons <- cbind(tabl.subset[S4Vectors::queryHits(hits3),], data.frame(unlist(exonsEnsDb)[S4Vectors::subjectHits(hits3)]), EnsID = names(unlist(exonsEnsDb))[S4Vectors::subjectHits(hits3)])
    annot <- biomaRt::select(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, key = "GENEID", keys = sample.genes.exons$EnsID, columns = c("SYMBOL", "GENEID"))
    sample.genes.exons <- merge(sample.genes.exons, annot, by.x = "EnsID", by.y = "GENEID", all.x = TRUE)
    sample.genes.exons <- subset(sample.genes.exons, !duplicated(sample.genes.exons[,c(11:15) * -1]))

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(sample.genes.exons, paste0(output.dir, "/Sniffles.", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(sample.genes.exons, file = paste0(output.dir, "/Sniffles.", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

    writexl::write_xlsx(sample.genes, paste0(output.dir, "/Sniffles.", sample.name.output, ".germline.SVs.genes.xlsx"))
    write.csv(sample.genes, file = paste0(output.dir, "/Sniffles.", sample.name.output, ".germline.SVs.genes.csv"), row.names = FALSE, quote = FALSE)

}


#! nanogladiator
annotateGladiator <- function(sample.file, sample.name, w, output.dir, sample.name.output = NA, shift.thr = 0.4, freq.error.thr = 0.4){
    nanogladiator.results <- fread(sample.file)
    colnames(nanogladiator.results) <- c("seqnames", "start", "end", "boh1", "CNV.Type", "Copies", "Freq", "Freq.Error", "Shift")
    nanogladiator.results[, seqnames :=  gsub("chr", "", nanogladiator.results$seqnames)]
    nanogladiator.results.filt <- subset(nanogladiator.results, nanogladiator.results$Freq - nanogladiator.results$Freq.Error >= freq.error.thr)
    nanogladiator.results.filt <- subset(nanogladiator.results.filt, abs(nanogladiator.results.filt$Shift) >= shift.thr)
    nanogladiator.results.gr <- GenomicRanges::makeGRangesFromDataFrame(nanogladiator.results.filt)
    hits <- IRanges::findOverlaps(nanogladiator.results.gr, genesEnsDb)
    nanogladiator.annotated.genes <- cbind(nanogladiator.results.filt[S4Vectors::queryHits(hits),], start = data.frame(genesEnsDb[S4Vectors::subjectHits(hits)])$start, end = data.frame(genesEnsDb[S4Vectors::subjectHits(hits)])$end, symbol = data.frame(genesEnsDb[S4Vectors::subjectHits(hits)])$symbol,  gene_biotype = data.frame(genesEnsDb[S4Vectors::subjectHits(hits)])$gene_biotype )
    nanogladiator.annotated.genes.protein.coding <- subset(nanogladiator.annotated.genes, nanogladiator.annotated.genes$gene_biotype == "protein_coding")

    if(is.na(sample.name.output)) sample.name.output <- sample.name

    writexl::write_xlsx(nanogladiator.annotated.genes.protein.coding, paste0(output.dir, "/NanoGladiator.", w, ".", sample.name.output, ".germline.SVs.exons.xlsx"))
    write.csv(nanogladiator.annotated.genes.protein.coding, file = paste0(output.dir, "/NanoGladiator.", w, ".", sample.name.output, ".germline.SVs.exons.csv"), row.names = FALSE, quote = FALSE)

}
