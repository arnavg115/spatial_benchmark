#!/bin/bash
# Mouse brain benchmark (StrataMap vs. Trekker, same animal): raw data -> every processed artifact
# 03_generate_figures.sh needs. Same stages and the same analysis scripts as the breast benchmark's
# 02_preprocess.sh - only the samples, the RCTD reference (Allen Brain Cell Atlas, not Janesick
# Chromium Flex) and the shared gene axis (genes passing QC in both samples, not a targeted panel)
# differ. code/R/plot_common.R picks up this benchmark's registry.R via BENCHMARK=mouse_brain
# (exported by config.sh).
#
# Run as: bash benchmarks/mouse_brain/02_preprocess.sh, after setting the MB_* inputs in config.sh.

set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/config.sh"

for v in MB_TREKKER_COORDS MB_TREKKER_MICRONS_PER_UNIT; do
  if [ -z "${!v}" ]; then echo "$v is not set - see benchmarks/mouse_brain/RUNNING.md" >&2; exit 1; fi
done
for f in matrix.mtx.gz features.tsv.gz barcodes.tsv.gz; do
  [ -f "$MB_STRATAMAP_DIR/Raw_Matrix_Files/$f" ] || { echo "StrataMap file not found: $MB_STRATAMAP_DIR/Raw_Matrix_Files/$f" >&2; exit 1; }
done
[ -f "$MB_STRATAMAP_SEG_CSV" ] ||{ echo "StrataMap segmentation CSV not found: $MB_STRATAMAP_SEG_CSV" >&2; exit 1; }
[ -f "$MB_TREKKER_H5AD" ] || { echo "Trekker h5ad not found: $MB_TREKKER_H5AD" >&2; exit 1; }

# ============================================================ STAGE 0: ABC atlas reference
echo "=== STAGE 0: Allen Brain Cell Atlas reference ==="
if [ ! -f "$MB_REFERENCE_H5AD" ]; then
  dispatch mb_build_reference 64G 08:00:00 4 \
    "$PYTHON_CMD $CODE_PY/build_reference_abc_mouse.py --out $MB_REFERENCE_H5AD --palette-out $MB_PALETTE_CSV \
      --download-dir $REPRO_DATA/raw/abc_atlas --regions $MB_REFERENCE_REGIONS"
  wait_all
fi

# ============================================================ STAGE 1: ingestion
echo "=== STAGE 1: ingestion ==="
dispatch mb_aggregate_stratamap "$MB_MEM_LARGE" 04:00:00 4 \
  "$PYTHON_CMD $CODE_PY/aggregate_stratamap_barcodes_to_cells.py \
    --raw-dir $MB_STRATAMAP_DIR/Raw_Matrix_Files --segmentation-csv $MB_STRATAMAP_SEG_CSV \
    --out-dir $(ingest_dir_for mouse_stratamap)"
dispatch mb_convert_trekker "$MB_MEM_SMALL" 01:00:00 2 \
  "$PYTHON_CMD $CODE_PY/convert_trekker_h5ad.py --h5ad $MB_TREKKER_H5AD --out-dir $(ingest_dir_for mouse_trekker) \
    --counts $MB_TREKKER_COUNTS --gene-symbols $MB_TREKKER_GENE_SYMBOLS \
    --coords $MB_TREKKER_COORDS --microns-per-unit $MB_TREKKER_MICRONS_PER_UNIT"
wait_all

# ============================================================ STAGE 2: QC (same nCount >= 10 /
# gene in >= 10 cells thresholds as the breast benchmark)
echo "=== STAGE 2: QC ==="
dispatch mb_qc_stratamap "$MB_MEM_LARGE" 01:30:00 4 \
  "$R_CMD $CODE_R/qc_seurat_stratamap.R --sample-name mouse_stratamap \
    --aggregated-dir $(ingest_dir_for mouse_stratamap) --out-dir $(qc_dir_for mouse_stratamap)"
dispatch mb_qc_trekker "$MB_MEM_SMALL" 01:00:00 4 \
  "$R_CMD $CODE_R/qc_seurat.R --technology mouse_trekker --sample mouse_trekker --level cell \
    --input-type mtxdir --counts-path $(ingest_dir_for mouse_trekker) \
    --metadata-path $(ingest_dir_for mouse_trekker)/cell_metadata.parquet --out-dir $(qc_dir_for mouse_trekker)"
wait_all

# STAGE 2b: shared gene axis = genes passing QC in both samples. Feeds every shared_panel_genes step.
dispatch mb_shared_genes "$MB_MEM_SMALL" 00:30:00 2 \
  "$R_CMD $CODE_R/shared_gene_list.R --qc-rds $(qc_rds_for mouse_stratamap),$(qc_rds_for mouse_trekker) \
    --out-csv $MB_SHARED_GENES_CSV"
wait_all

# ============================================================ STAGE 3: RCTD
echo "=== STAGE 3: RCTD ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_rctd_$TECH" "$MB_MEM_LARGE" 04:00:00 8 \
    "$PYTHON_CMD $CODE_PY/run_rctd_doublet.py --h5 $(qc_h5_for $TECH) --reference-h5ad $MB_REFERENCE_H5AD --out-dir $(rctd_dir_for $TECH)"
done
wait_all

# ============================================================ STAGE 4: embeddings - self ("all")
# and the shared gene axis (filtered the same way as the breast benchmark's STAGE 4c).
echo "=== STAGE 4: gene-set embeddings ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_embed_$TECH" "$MB_MEM_LARGE" 02:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample $TECH \
      --qc-rds $(qc_rds_for $TECH) --gene-sets all \
      --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH) --umap-init $MB_UMAP_INIT"
  dispatch "mb_embed_shared_$TECH" "$MB_MEM_LARGE" 02:00:00 8 \
    "$R_CMD $CODE_R/gene_set_embeddings.R --technology $TECH --sample $TECH \
      --qc-rds $(qc_rds_for $TECH) --gene-sets shared_panel_genes --min-count-in-subset 10 --npcs 15 --label-suffix _filtered \
      --rctd-dir $(rctd_dir_for $TECH) --out-dir $(embed_dir_for $TECH) --umap-init $MB_UMAP_INIT"
done
wait_all

# ============================================================ STAGE 5: spatial metadata
echo "=== STAGE 5: spatial metadata ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_spatial_$TECH" "$MB_MEM_SMALL" 00:45:00 2 \
    "$R_CMD $CODE_R/extract_spatial_metadata.R --technology $TECH --level cell \
      --qc-rds $(qc_rds_for $TECH) --rctd-dir $(rctd_dir_for $TECH) --out-dir $REPRO_RESULTS/spatial_metadata"
done
wait_all

# ============================================================ STAGE 6: spillover cache (figs 6/7).
# For Trekker (dissociated nuclei) a w2-vs-neighborhood signal reflects ambient RNA/doublets rather
# than segmentation spillover - see README.md.
echo "=== STAGE 6: spillover ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_spillover_$TECH" "$MB_MEM_LARGE" 04:00:00 8 \
    "$R_CMD $CODE_R/spillover.R --technology $TECH \
      --qc-rds $(qc_rds_for $TECH) --rctd-dir $(rctd_dir_for $TECH) \
      --cache-dir $REPRO_RESULTS/spillover_cache/$TECH"
done
wait_all

# ============================================================ STAGE 7: composition summary (fig 3)
echo "=== STAGE 7: composition summary ==="
dispatch mb_composition 16G 00:15:00 2 \
  "$R_CMD $CODE_R/celltype_composition.R --out-dir $REPRO_RESULTS/composition"
wait_all

# ============================================================ STAGE 8: separation metric + ARI (fig 5)
echo "=== STAGE 8: separation metric (ARI) ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_sepmetric_$TECH" "$MB_MEM_SMALL" 00:30:00 4 \
    "$R_CMD $CODE_R/separation_metric.R --technology $TECH \
      --embeddings-dir $(embed_dir_for $TECH) --out-csv $REPRO_RESULTS/separation_metric/separation_metric_${TECH}_cell.csv"
done
wait_all
dispatch mb_separation_ari 16G 00:15:00 2 \
  "$R_CMD $CODE_R/separation_ari_comparison.R --results-dir $REPRO_RESULTS/separation_metric --out-dir $REPRO_RESULTS/separation_ari"
wait_all

# ============================================================ STAGE 9: count distributions (figs 4/5)
echo "=== STAGE 9: count distribution extraction ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_extract_counts_$TECH" "$MB_MEM_LARGE" 01:00:00 4 \
    "$R_CMD $CODE_R/extract_count_distributions.R --technology $TECH --out-dir $REPRO_RESULTS/count_distributions"
done
wait_all

# ============================================================ STAGE 10: N genes / N cells (fig 5)
echo "=== STAGE 10: ngenes/ncells barplot data ==="
dispatch mb_ngenes_ncells 16G 00:15:00 2 \
  "$R_CMD $CODE_R/ngenes_ncells_barplot.R \
    --extract-dir $REPRO_RESULTS/count_distributions \
    --ari-csv $REPRO_RESULTS/separation_ari/separation_ari_comparison.csv \
    --out-dir $REPRO_RESULTS/ngenes_ncells_barplot"
wait_all

# ============================================================ STAGE 11: per-gene-set sensitivity (fig 5)
echo "=== STAGE 11: sensitivity by gene set ==="
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_sensitivity_$TECH" "$MB_MEM_LARGE" 02:00:00 4 \
    "$R_CMD $CODE_R/sensitivity_by_geneset_celltype.R --technology $TECH --sample $TECH \
      --qc-rds $(qc_rds_for $TECH) --gene-sets all,shared_panel_genes --min-count-in-subset 10 \
      --out-dir $REPRO_RESULTS/sensitivity/$TECH"
done
wait_all
dispatch mb_sensitivity_synthesis 16G 00:15:00 2 \
  "$R_CMD $CODE_R/sensitivity_synthesis.R --sensitivity-dir $REPRO_RESULTS/sensitivity --out-dir $REPRO_RESULTS/sensitivity"
wait_all

# ============================================================ STAGE 12: per-gene sensitivity (fig 8)
echo "=== STAGE 12: per-gene sensitivity ==="
GS_OUT="$REPRO_RESULTS/gene_sensitivity"
mkdir -p "$GS_OUT"
for TECH in "${SAMPLES[@]}"; do
  dispatch "mb_gene_sensitivity_$TECH" "$MB_MEM_LARGE" 01:00:00 4 \
    "$R_CMD $CODE_R/gene_sensitivity.R --technology $TECH --out-dir $GS_OUT"
done
wait_all
COMBINED="$GS_OUT/gene_sensitivity_combined.csv"
FIRST=1
for TECH in "${SAMPLES[@]}"; do
  f="$GS_OUT/${TECH}_gene_sensitivity.csv"
  if [ "$FIRST" -eq 1 ]; then cat "$f" > "$COMBINED"; FIRST=0; else tail -n +2 "$f" >> "$COMBINED"; fi
done
dispatch mb_gene_sensitivity_synthesis 16G 00:15:00 2 \
  "$R_CMD $CODE_R/gene_sensitivity_synthesis.R --sensitivity-dir $GS_OUT --out-dir $GS_OUT"
wait_all

echo "=== Preprocessing done. Run benchmarks/mouse_brain/03_generate_figures.sh next. ==="
