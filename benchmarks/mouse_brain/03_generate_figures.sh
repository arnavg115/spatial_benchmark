#!/bin/bash
# Mouse brain benchmark figures - the same 8 plotting scripts as the breast benchmark's
# 03_generate_figures.sh, run against this benchmark's results (BENCHMARK=mouse_brain, exported by
# config.sh). Run after benchmarks/mouse_brain/02_preprocess.sh.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

OUT="$REPRO_RESULTS/figures"

dispatch mb_fig1_spatial_grid 32G 00:20:00 2 \
  "$R_CMD $CODE_R/spatial_grid.R --out-dir $OUT/spatial --format png --dpi 300 \
    --point-size 0.15 --alpha 0.7 --panel-width 5 --panel-height 5 --out-name spatial_grid_all"

dispatch mb_fig2_umap_by_geneset 16G 00:20:00 2 \
  "$R_CMD $CODE_R/umap_by_geneset.R --out-dir $OUT/umap"
dispatch mb_fig2_umap_self_grid 16G 00:20:00 2 \
  "$R_CMD $CODE_R/umap_self_grid.R --out-dir $OUT/umap --out-name umap_self_grid_all"

dispatch mb_fig3_composition_grid 16G 00:15:00 2 \
  "$R_CMD $CODE_R/composition_grid.R --summary-csv $REPRO_RESULTS/composition/composition_summary.csv --out-dir $OUT/composition/sorted_by_abundance"

dispatch mb_fig4_threshold_sensitivity 16G 00:20:00 2 \
  "$R_CMD $CODE_R/celltype_threshold_sensitivity.R --extract-dir $REPRO_RESULTS/count_distributions --out-dir $OUT/threshold_sensitivity --layout row_tall"

dispatch mb_fig5_summary_grid 16G 00:20:00 2 \
  "$R_CMD $CODE_R/summary_grid_shared_panel.R \
    --ngenes-csv $REPRO_RESULTS/ngenes_ncells_barplot/ngenes_barplot.csv \
    --ncells-csv $REPRO_RESULTS/ngenes_ncells_barplot/ncells_barplot.csv \
    --overall-csv $REPRO_RESULTS/sensitivity/sensitivity_overall_all_technologies.csv \
    --ari-csv $REPRO_RESULTS/separation_ari/separation_ari_comparison.csv \
    --out-dir $OUT/summary_grid"

dispatch mb_fig6_w2_boxplot 16G 00:15:00 2 \
  "$R_CMD $CODE_R/w2_boxplot.R --cache-root $REPRO_RESULTS/spillover_cache --out-dir $OUT/w2"

dispatch mb_fig7_spillover_metrics 16G 00:20:00 2 \
  "$R_CMD $CODE_R/spillover_metrics.R --cache-root $REPRO_RESULTS/spillover_cache --out-dir $OUT/spillover"

dispatch mb_fig8_gene_sensitivity_scatter 16G 00:15:00 2 \
  "$R_CMD $CODE_R/gene_sensitivity_scatter.R --stats-csv $REPRO_RESULTS/gene_sensitivity/csv/gene_sensitivity_stats.csv --out-dir $OUT/gene_sensitivity_scatter"

wait_all

echo "Done. Figures under $OUT/ (same filenames as the breast benchmark's 03_generate_figures.sh)."
