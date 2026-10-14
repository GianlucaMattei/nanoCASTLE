# Licensing, third-party software and data

## Licensing of the nanoCASTLE image

| Component | Licence |
|---|---|
| nanoCASTLE pipeline and annotation code (`/opt/nanomatizer/Nanomatizer`) | MIT — `Nanomatizer/LICENSE` |
| NanoGLADIATOR 1.0 (`/opt/nanomatizer/NanoGLADIATOR_1.0`) | PolyForm Noncommercial 1.0.0 — `NanoGLADIATOR_1.0/LICENSE` |
| GASOLINE 2.0 (`/opt/nanomatizer/GASOLINE_2.0`) | PolyForm Noncommercial 1.0.0 — `GASOLINE_2.0/LICENSE` |
| everything else | its own licence (tables below) |

NanoGLADIATOR and GASOLINE may be used for noncommercial purposes only;
commercial use requires a separate licence from the copyright holder
(Azienda Ospedaliero-Universitaria Meyer IRCCS). nanoCASTLE itself does not need them: for
commercial use run the pipeline with the other callers only, e.g.
`-t sniffles,cutesv`.

## Third-party software and data

The nanoCASTLE image bundles third-party tools and public annotation
resources. Each of them remains under its own licence / terms of use; please
cite the original authors when you use the corresponding results.

### Software

| Component | Where in the image | Source |
|---|---|---|
| NanoGLADIATOR 1.0 (trimmed: code + window targets only) | `/opt/nanomatizer/NanoGLADIATOR_1.0` | Magi et al., *Nano-GLADIATOR*, Bioinformatics 2019 — PolyForm Noncommercial 1.0.0 |
| GASOLINE 2.0 | `/opt/nanomatizer/GASOLINE_2.0` | Magi et al. — PolyForm Noncommercial 1.0.0 |
| Sniffles2 2.8.0 | conda env `sv-sniffles` | https://github.com/fritzsedlazeck/Sniffles (MIT) |
| cuteSV 2.1.4 | conda env `sv-cutesv` | https://github.com/tjiangHIT/cuteSV (MIT) |
| samtools 1.24, bedtools 2.31.1, minimap2, R, Bioconductor/CRAN packages, Python | conda envs in `/opt/nanomatizer/mamba-root/envs` | conda-forge / bioconda, each under its own licence (see `conda-meta/` in each env) |
| micromamba | `/opt/micromamba` | https://github.com/mamba-org/mamba (BSD-3-Clause) |
| Ubuntu 24.04 base system | `/` | https://ubuntu.com/legal/open-source-licences |

The exact package list of every conda environment is in
`/opt/nanomatizer/Nanomatizer/docker/envspecs/`.

### Annotation data shipped in the image (`/opt/nanomatizer/Nanomatizer/Databsese`)

| File | Resource | Terms |
|---|---|---|
| `gnomad.v4.1.cnv.all.vcf`, `gnomad.v4.1.sv.sites.selected.columns.bed` | gnomAD v4.1 SV / CNV | CC0 — https://gnomad.broadinstitute.org/policies |
| `GRCh38.nr_deletions.pathogenic.tsv`, `GRCh38.nr_duplications.somatic.tsv`, `GRCh38.nr_insertions.somatic.tsv` | NCBI dbVar non-redundant SV sets | NCBI public data — https://www.ncbi.nlm.nih.gov/home/about/policies/ |
| `Clingen-Dosage-Sensitivity-2025-08-30.csv` | ClinGen Dosage Sensitivity | CC0 — https://clinicalgenome.org/docs/terms-of-use/ |
| `9606.FIXED.protein.links.v12.0.txt` | STRING v12.0 (reformatted) | CC BY 4.0 — https://string-db.org/cgi/access |
| `h.all.v2024.1.Hs.symbols.CON.MACROCATEGORIE.liberzon.et.al.gmt` | MSigDB Hallmark gene sets v2024.1 + process categories from Liberzon et al. 2015 | CC BY 4.0 — https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp |
| `ENCODE.blacklist.CNV.bed` | ENCODE blacklist regions | https://www.encodeproject.org |
| `EnsemblRegulatory/Homo_sapiens.GRCh38.regulatory_features.v112.gff3` | Ensembl Regulatory Build, release 112 | https://www.ensembl.org/info/about/legal/disclaimer.html |
| `GASOLINE_2.0/DB/*.rds` | GENCODE, HGSVC, IGSR (1000 Genomes), gnomAD — reformatted for GASOLINE | open data, see the respective projects |

### Data NOT shipped (downloaded by the user)

These resources are fetched from the original provider with
`Nanomatizer/download_databases.sh` and mounted in `/opt/nanomatizer/userdb`.
By downloading them you accept the provider's terms of use.

| File | Resource | Terms |
|---|---|---|
| `OncoKB.tsv` | OncoKB Cancer Gene List (MSKCC) | https://www.oncokb.org/terms — redistribution not permitted without a licence |
| `population_cnv_grch38.txt` | DECIPHER population CNV frequencies | https://www.deciphergenomics.org/about/downloads/data — redistribution requires an agreement with DECIPHER |
| `trrust_rawdata.human.tsv` | TRRUST v2 | https://www.grnpedia.org/trrust/ — see the provider's terms of use |
| `DGV.GS.hg38.gff3` | Database of Genomic Variants, Gold Standard (hg38) | http://dgv.tcag.ca — see the provider's terms of use |
| `CRISPRGeneEffect.csv`, `Model.csv` (optional) | DepMap Public (Broad Institute) | CC BY 4.0 — https://depmap.org/portal/data_page/ |
