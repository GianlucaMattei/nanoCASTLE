# =============================================================================
# nanocastle_genome.r — reference-dependent gene annotation for step 6.
# NANOCASTLE_ASSEMBLY = hg38 (default) | chm13, set by run_pipeline_sequential.sh
# from the genome fingerprint of the reference.
#   hg38  : EnsDb.Hsapiens.v86 (as in 1.0.0)
#   chm13 : the same EnsDb.Hsapiens.v86 genes/CDS/exons lifted to T2T-CHM13v2.0
#           (EnsDb.Hsapiens.v86_chm13.rds built by build_chm13_dbs.r, in
#           NANOCASTLE_DBDIR). Same gene ids, symbols and seqnames style ("1", not "chr1").
# Ensembl regulatory build: same rule (GRCh38 gff3 or its liftover in NANOCASTLE_DBDIR).
# =============================================================================
NANOCASTLE_ASSEMBLY <- Sys.getenv("NANOCASTLE_ASSEMBLY", "hg38")
if (!NANOCASTLE_ASSEMBLY %in% c("hg38", "chm13")) stop("NANOCASTLE_ASSEMBLY must be hg38 or chm13, not ", NANOCASTLE_ASSEMBLY)
.nc_chm13 <- NULL
.nc_chm13_obj <- function() {
    if (is.null(.nc_chm13)) {
        f <- file.path(Sys.getenv("NANOCASTLE_DBDIR", "/opt/nanomatizer/Nanomatizer/Databsese_chm13"),
                       "EnsDb.Hsapiens.v86_chm13.rds")
        if (!file.exists(f)) stop("chm13 gene annotation missing: ", f, " (run build_chm13_dbs.r)")
        .nc_chm13 <<- readRDS(f)
    }
    .nc_chm13
}
ens_genes <- function() {
    if (NANOCASTLE_ASSEMBLY == "chm13") .nc_chm13_obj()$genes
    else GenomicFeatures::genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
}
ens_cds_by_gene <- function() {
    if (NANOCASTLE_ASSEMBLY == "chm13") .nc_chm13_obj()$cds
    else GenomicFeatures::cdsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by = "gene")
}
ens_exons_by_gene <- function() {
    if (NANOCASTLE_ASSEMBLY == "chm13") .nc_chm13_obj()$exons
    else GenomicFeatures::exonsBy(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86, by = "gene")
}
regulatory_gff3 <- function(t_root) {
    if (NANOCASTLE_ASSEMBLY == "chm13")
        file.path(Sys.getenv("NANOCASTLE_DBDIR", "/opt/nanomatizer/Nanomatizer/Databsese_chm13"),
                  "EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3")
    else file.path(t_root, "Nanomatizer/Databsese/EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3")
}
