#!/usr/bin/env Rscript
# Shared gene axis for a benchmark with no common targeted panel (e.g. two whole-transcriptome
# technologies): every gene that passed QC (detected in >= --min-cells-per-gene cells, qc_seurat*.R)
# in EVERY listed sample. Written in the same one-column gene_symbol format as
# reference/shared_panel_gene_list.csv, so gene_list_for_set("shared_panel_genes", ...) reads it
# unchanged once the benchmark's registry points SHARED_PANEL_GENE_LIST_PATH at it.
#
# Usage: Rscript shared_gene_list.R --qc-rds <a>_qc_seurat.rds,<b>_qc_seurat.rds --out-csv <path>

suppressPackageStartupMessages({
  library(Seurat)
  library(optparse)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--qc-rds", type = "character", help = "comma-separated QCed Seurat .rds paths"),
  make_option("--out-csv", type = "character")
)))

paths <- strsplit(opt$`qc-rds`, ",")[[1]]
stopifnot("--qc-rds needs at least 2 samples" = length(paths) >= 2)

shared <- NULL
for (p in paths) {
  genes <- rownames(readRDS(p))
  cat(basename(p), ":", length(genes), "genes after QC\n")
  shared <- if (is.null(shared)) genes else intersect(shared, genes)
}
shared <- sort(shared)
cat("Shared across all", length(paths), "samples:", length(shared), "genes\n")
stopifnot("No genes shared across samples - check gene symbols match (e.g. mouse 'Gfap' vs human 'GFAP')" = length(shared) > 0)

dir.create(dirname(opt$`out-csv`), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(gene_symbol = shared), opt$`out-csv`, row.names = FALSE)
cat("Saved:", opt$`out-csv`, "\n")
