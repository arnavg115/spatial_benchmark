# Mouse brain benchmark (StrataMap vs. Trekker, same animal) - sourced at the end of
# code/R/plot_common.R when BENCHMARK=mouse_brain, overriding the breast benchmark's sample list,
# labels, colors, path helpers and gene-list paths. Everything else in plot_common.R (themes, scale
# bars, save_post(), ...) is shared unchanged, as is every analysis script that sources it.
#
# The path helpers below must stay in sync with benchmarks/mouse_brain/config.sh, which builds the
# same paths on the shell side.

MB_PROCESSED <- file.path(REPRO_DATA, "processed/mouse_brain")

TECH_DISPLAY <- c(mouse_stratamap = "StrataMap", mouse_trekker = "Trekker")

label_for_tech <- function(tech, samp = NA_character_) {
  if (tech %in% names(TECH_DISPLAY)) return(TECH_DISPLAY[[tech]])
  tech
}

CELL_LEVEL_SAMPLES <- data.frame(
  technology = c("mouse_stratamap", "mouse_trekker"),
  is_stratamap = c(TRUE, FALSE),
  stringsAsFactors = FALSE
)
CELL_LEVEL_SAMPLES$label <- vapply(CELL_LEVEL_SAMPLES$technology, label_for_tech, character(1))

TECH_ORDER <- CELL_LEVEL_SAMPLES$technology

# StrataMap keeps the breast benchmark's StrataMap color; Trekker gets a new one.
SAMPLE_COLORS <- c("StrataMap" = "#A9906F", "Trekker" = "#3B7DB8")
tech_family <- function(tech) label_for_tech(tech)
TECH_FAMILY_COLORS <- SAMPLE_COLORS

SPATIAL_FLIP_XY <- character(0)

qc_rds_path <- function(technology) {
  file.path(MB_PROCESSED, "std_seurat_analysis", technology, paste0(technology, "_qc_seurat.rds"))
}

rctd_dir_path <- function(technology) {
  file.path(MB_PROCESSED, "cell_type_annotation", technology, "abc_wmb_class/rctd_class_aware/Level2/single_cell")
}

reductions_path <- function(technology, gene_set) {
  file.path(MB_PROCESSED, "std_seurat_analysis", technology, "embeddings", paste0("reductions_", gene_set, ".rds"))
}

# Allen Brain Cell Atlas class colors, written by code/python/build_reference_abc_mouse.py next to
# the reference it builds (RCTD's Level2 = ABC "class", so first_type values match these names).
CELLTYPE_PALETTE_PATH <- file.path(REPRO_DATA, "reference/mouse_brain/abc_wmb_class_palette.csv")

# "Shared panel genes" here = every gene passing QC in BOTH samples (code/R/shared_gene_list.R,
# 02_preprocess.sh STAGE 2b) - there's no targeted panel, both are whole-transcriptome.
SHARED_PANEL_GENE_LIST_PATH <- file.path(REPRO_RESULTS, "shared_genes/shared_panel_gene_list.csv")
PANEL_GENE_LIST_PATH <- NA_character_
shared_gene_set_for <- function(tech) "shared_panel_genes"

SEPARATION_METRIC_CSV <- setNames(paste0("separation_metric_", TECH_ORDER, "_cell.csv"), TECH_ORDER)

# Figure 4's single-row layout is sized for 8 panels in the breast benchmark.
POST_FIGURE_SIZES$threshold_sensitivity_row_tall$width <- 3.0 * length(TECH_ORDER)
