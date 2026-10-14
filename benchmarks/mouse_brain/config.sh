# Mouse brain benchmark config, sourced by this directory's 02_preprocess.sh and
# 03_generate_figures.sh. Set the MB_* input variables below (export them, or edit the defaults)
# before running - see RUNNING.md in this directory.

# ------------------------------------------------------------ shared setup
export BENCHMARK=mouse_brain
_MB_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Always its own results dir, never an inherited REPRO_RESULTS - mixing with the breast benchmark's
# results would pull breast samples into the cell-type palette and every cross-sample summary.
export REPRO_RESULTS="${MB_RESULTS:-$(cd "$_MB_DIR/../.." && pwd)/results/mouse_brain}"
source "$_MB_DIR/../../code/lib.sh"

# ------------------------------------------------------------ inputs
# Default locations (after lib.sh has set REPRO_DATA) - put the data here, or point these variables
# elsewhere. See RUNNING.md.
MB_RAW="$REPRO_DATA/raw/mouse_brain"

# StrataMap BaseSpace download for the mouse brain sample - same layout as the breast samples:
# $MB_STRATAMAP_DIR/Raw_Matrix_Files/{matrix.mtx.gz,features.tsv.gz,barcodes.tsv.gz} +
# $MB_STRATAMAP_DIR/Segmentations/registered_Expanded_5um_cell_contour_coords_local.csv
MB_STRATAMAP_DIR="${MB_STRATAMAP_DIR:-$MB_RAW/stratamap}"
MB_STRATAMAP_SEG_CSV="${MB_STRATAMAP_SEG_CSV:-$MB_STRATAMAP_DIR/Segmentations/registered_Expanded_5um_cell_contour_coords_local.csv}"

# Trekker .h5ad - run `python code/python/convert_trekker_h5ad.py --h5ad <file> --inspect` first and
# set the 4 flags below from what it prints.
MB_TREKKER_H5AD="${MB_TREKKER_H5AD:-$MB_RAW/trekker/trekker.h5ad}"
MB_TREKKER_COUNTS="${MB_TREKKER_COUNTS:-X}"                    # X | raw | layer:<name> (raw integer UMIs)
MB_TREKKER_GENE_SYMBOLS="${MB_TREKKER_GENE_SYMBOLS:-var_names}" # var_names | var:<column>
MB_TREKKER_COORDS="${MB_TREKKER_COORDS:-}"                      # obsm:<key> | obs:<x>,<y>
MB_TREKKER_MICRONS_PER_UNIT="${MB_TREKKER_MICRONS_PER_UNIT:-}"  # 1 if coordinates are already microns

# ABC atlas regions sampled for the RCTD reference (code/python/build_reference_abc_mouse.py):
# 'all' until the tissue's region is known, then e.g. "Isocortex-1,Isocortex-2,HPF".
MB_REFERENCE_REGIONS="${MB_REFERENCE_REGIONS:-all}"

# 'random' if UMAP's spectral init segfaults on very large samples (see gene_set_embeddings.R).
MB_UMAP_INIT="${MB_UMAP_INIT:-spectral}"

# Per-step memory. Unknown sample sizes - the defaults match the breast benchmark's mid-sized
# StrataMap sample (Grade2); raise them for large sections.
MB_MEM_LARGE="${MB_MEM_LARGE:-96G}"
MB_MEM_SMALL="${MB_MEM_SMALL:-32G}"

SAMPLES=(mouse_stratamap mouse_trekker)

# Same paths as registry.R's path helpers.
MB_PROCESSED="$REPRO_DATA/processed/mouse_brain"
MB_REFERENCE_H5AD="$REPRO_DATA/reference/mouse_brain/abc_wmb_class_reference.h5ad"
MB_PALETTE_CSV="$REPRO_DATA/reference/mouse_brain/abc_wmb_class_palette.csv"
MB_SHARED_GENES_CSV="$REPRO_RESULTS/shared_genes/shared_panel_gene_list.csv"
qc_dir_for()    { echo "$MB_PROCESSED/std_seurat_analysis/$1"; }
qc_rds_for()    { echo "$(qc_dir_for "$1")/$1_qc_seurat.rds"; }
qc_h5_for()     { echo "$(qc_dir_for "$1")/$1_qc_counts.h5"; }
rctd_dir_for()  { echo "$MB_PROCESSED/cell_type_annotation/$1/abc_wmb_class/rctd_class_aware/Level2/single_cell"; }
embed_dir_for() { echo "$(qc_dir_for "$1")/embeddings"; }
ingest_dir_for(){ echo "$MB_PROCESSED/segmentation/$1"; }
