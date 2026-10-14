# Running the mouse brain benchmark (StrataMap vs. Trekker)

Step-by-step instructions. For what the benchmark compares and how to interpret the results, see
[README.md](README.md).

All commands are run from the **repository root** (`atera_vs_stratamap/`).

---

## 1. One-time setup

```bash
conda env create -f environment.yml
conda activate spatial-tech-comparison
R -e 'remotes::install_github("dmcable/spacexr")'
R -e 'remotes::install_github("bdsc-tds/SPLIT")'
```

Skip this if you already set up the environment for the breast benchmark.

---

## 2. Choose a data directory

Everything large (your raw data, the reference download, intermediate files) lives under one
directory, `$REPRO_DATA`. Pick a disk with **at least ~150 GB free**: the full-brain reference
download is ~95 GB, plus your two samples and intermediates.

```bash
export REPRO_DATA=/path/with/space
```

If you don't set it, it defaults to `./data` inside the repo.

---

## 3. Put your data here

```
$REPRO_DATA/raw/mouse_brain/
├── stratamap/
│   ├── Raw_Matrix_Files/
│   │   ├── matrix.mtx.gz
│   │   ├── features.tsv.gz
│   │   └── barcodes.tsv.gz
│   └── Segmentations/
│       └── registered_Expanded_5um_cell_contour_coords_local.csv
└── trekker/
    └── trekker.h5ad
```

- **StrataMap**: copy the sample folder from your BaseSpace download (the same layout as the breast
  StrataMap samples) and rename it to `stratamap/`. Other files in it (`Images/`, nuclei
  segmentation, ...) can stay; they're ignored.
- **Trekker**: rename your file to `trekker.h5ad`.

Already have the data elsewhere and don't want to move it? Point the pipeline at it instead:

```bash
export MB_STRATAMAP_DIR=/existing/path/to/stratamap_sample_folder
export MB_TREKKER_H5AD=/existing/path/to/your_file.h5ad
```

**Quick check.** Column 2 of the StrataMap features file should hold mouse gene symbols (e.g. `Gfap`,
`Snap25`):

```bash
zcat < $REPRO_DATA/raw/mouse_brain/stratamap/Raw_Matrix_Files/features.tsv.gz | head -3
```

---

## 4. Inspect the Trekker file and set 4 variables

Trekker h5ad files differ in where they store counts, gene names and coordinates, so the pipeline
asks you to say where. Run:

```bash
python code/python/convert_trekker_h5ad.py --h5ad $REPRO_DATA/raw/mouse_brain/trekker/trekker.h5ad --inspect
```

Then set these from what it prints:

| Variable | What to look for in the output | Example |
|---|---|---|
| `MB_TREKKER_COUNTS` | The matrix line that says `integer-valued=True`: `X`, `raw`, or `layer:<name>`. **Required to be raw counts**; the script refuses normalized data. | `X` |
| `MB_TREKKER_GENE_SYMBOLS` | If `var_names` look like `Gfap`, use `var_names`. If they look like `ENSMUSG...`, find the `var:` column that holds symbols. | `var_names` |
| `MB_TREKKER_COORDS` | The `obsm:` entry (or two `obs:` columns) holding x/y coordinates. | `obsm:spatial` or `obs:x,y` |
| `MB_TREKKER_MICRONS_PER_UNIT` | What to multiply the coordinates by to get microns. Compare the printed coordinate range with your section's physical size: a range of ~0–6000 on a ~6 mm section means they're already microns, so use `1`. | `1` |

```bash
export MB_TREKKER_COUNTS=X
export MB_TREKKER_GENE_SYMBOLS=var_names
export MB_TREKKER_COORDS=obsm:spatial
export MB_TREKKER_MICRONS_PER_UNIT=1
```

Get the units right: the spillover analysis uses a 15 µm neighbor radius, and the per-mm² metrics
assume microns. Wrong units won't crash anything; they'll just give wrong numbers.

---

## 5. Optional settings

| Variable | Default | When to change it |
|---|---|---|
| `MB_REFERENCE_REGIONS` | `all` | Once you know the brain region(s) in your sections. Restricting cuts the download and stops RCTD assigning cell classes from regions that aren't in your tissue. Choices: `CB, CTXsp, HPF, HY, Isocortex-1, Isocortex-2, MB, MY, OLF, P, PAL, STR, TH` (comma-separated, e.g. `Isocortex-1,Isocortex-2,HPF`). |
| `MB_MEM_LARGE` / `MB_MEM_SMALL` | `96G` / `32G` | If a step runs out of memory (large sections). |
| `MB_UMAP_INIT` | `spectral` | Set to `random` if the UMAP step segfaults (happens at ~500k+ cells). |
| `USE_SLURM` + `SLURM_ACCOUNT` | off | To run each heavy step as a SLURM job: `export USE_SLURM=true SLURM_ACCOUNT=<your account>`. Without SLURM, everything runs sequentially on the current machine. |
| `MB_RESULTS` | `./results/mouse_brain` | To put results somewhere else. |

---

## 6. Run

```bash
bash benchmarks/mouse_brain/02_preprocess.sh       # hours; downloads reference on first run
bash benchmarks/mouse_brain/03_generate_figures.sh # minutes
```

`02_preprocess.sh` checks that your input files exist and that the Trekker variables are set before
doing anything. If it stops right away, the error message names what's missing.

Re-running is safe: every step overwrites its own outputs. The reference is only built once. To
rebuild it (e.g. after changing `MB_REFERENCE_REGIONS`), delete `$REPRO_DATA/reference/mouse_brain/`.

---

## 7. What `02_preprocess.sh` runs

For reference, or for rerunning a single step by hand. Paths assume the defaults above. Every
command runs with `BENCHMARK=mouse_brain` and `REPRO_ROOT`/`REPRO_DATA`/`REPRO_RESULTS` set, which
`config.sh` does for you. To run one by hand, first do
`source benchmarks/mouse_brain/config.sh` in a bash shell.

`<sample>` is `mouse_stratamap` or `mouse_trekker`.

| Stage | Script | Arguments | Output |
|---|---|---|---|
| 0. Reference | `code/python/build_reference_abc_mouse.py` | `--out $REPRO_DATA/reference/mouse_brain/abc_wmb_class_reference.h5ad --palette-out $REPRO_DATA/reference/mouse_brain/abc_wmb_class_palette.csv --download-dir $REPRO_DATA/raw/abc_atlas --regions $MB_REFERENCE_REGIONS` | Reference h5ad + cell-type colors |
| 1. StrataMap → cells | `code/python/aggregate_stratamap_barcodes_to_cells.py` | `--raw-dir $MB_STRATAMAP_DIR/Raw_Matrix_Files --segmentation-csv $MB_STRATAMAP_SEG_CSV --out-dir $REPRO_DATA/processed/mouse_brain/segmentation/mouse_stratamap` | Cell × gene matrix |
| 1. Trekker → cells | `code/python/convert_trekker_h5ad.py` | `--h5ad $MB_TREKKER_H5AD --out-dir $REPRO_DATA/processed/mouse_brain/segmentation/mouse_trekker --counts $MB_TREKKER_COUNTS --gene-symbols $MB_TREKKER_GENE_SYMBOLS --coords $MB_TREKKER_COORDS --microns-per-unit $MB_TREKKER_MICRONS_PER_UNIT` | Cell × gene matrix |
| 2. QC (StrataMap) | `code/R/qc_seurat_stratamap.R` | `--sample-name mouse_stratamap --aggregated-dir .../segmentation/mouse_stratamap --out-dir .../std_seurat_analysis/mouse_stratamap` | `mouse_stratamap_qc_seurat.rds` + `_qc_counts.h5` |
| 2. QC (Trekker) | `code/R/qc_seurat.R` | `--technology mouse_trekker --sample mouse_trekker --level cell --input-type mtxdir --counts-path .../segmentation/mouse_trekker --metadata-path .../segmentation/mouse_trekker/cell_metadata.parquet --out-dir .../std_seurat_analysis/mouse_trekker` | same |
| 2b. Shared genes | `code/R/shared_gene_list.R` | `--qc-rds <stratamap qc rds>,<trekker qc rds> --out-csv results/mouse_brain/shared_genes/shared_panel_gene_list.csv` | Genes passing QC in both |
| 3. Cell types | `code/python/run_rctd_doublet.py` | `--h5 <sample qc h5> --reference-h5ad <reference h5ad> --out-dir .../cell_type_annotation/<sample>/abc_wmb_class/rctd_class_aware/Level2/single_cell` | RCTD results |
| 4. Embeddings | `code/R/gene_set_embeddings.R` | `--technology <sample> --sample <sample> --qc-rds ... --gene-sets all --rctd-dir ... --out-dir .../std_seurat_analysis/<sample>/embeddings`, then again with `--gene-sets shared_panel_genes --min-count-in-subset 10 --npcs 15 --label-suffix _filtered` | PCA/UMAP |
| 5. Spatial table | `code/R/extract_spatial_metadata.R` | `--technology <sample> --level cell --qc-rds ... --rctd-dir ... --out-dir results/mouse_brain/spatial_metadata` | x/y + cell type CSV |
| 6. Spillover | `code/R/spillover.R` | `--technology <sample> --qc-rds ... --rctd-dir ... --cache-dir results/mouse_brain/spillover_cache/<sample>` | Per-cell spillover table |
| 7. Composition | `code/R/celltype_composition.R` | `--out-dir results/mouse_brain/composition` | Fig 3 data |
| 8. ARI | `code/R/separation_metric.R`, then `separation_ari_comparison.R` | `--technology <sample> --embeddings-dir ... --out-csv results/mouse_brain/separation_metric/separation_metric_<sample>_cell.csv`; then `--results-dir results/mouse_brain/separation_metric --out-dir results/mouse_brain/separation_ari` | Fig 5 ARI |
| 9. Count distributions | `code/R/extract_count_distributions.R` | `--technology <sample> --out-dir results/mouse_brain/count_distributions` | Figs 4/5 data |
| 10. N genes / N cells | `code/R/ngenes_ncells_barplot.R` | `--extract-dir results/mouse_brain/count_distributions --ari-csv results/mouse_brain/separation_ari/separation_ari_comparison.csv --out-dir results/mouse_brain/ngenes_ncells_barplot` | Fig 5 data |
| 11. Sensitivity | `code/R/sensitivity_by_geneset_celltype.R`, then `sensitivity_synthesis.R` | `--technology <sample> --sample <sample> --qc-rds ... --gene-sets all,shared_panel_genes --min-count-in-subset 10 --out-dir results/mouse_brain/sensitivity/<sample>` | Fig 5 data |
| 12. Per-gene sensitivity | `code/R/gene_sensitivity.R`, then `gene_sensitivity_synthesis.R` | `--technology <sample> --out-dir results/mouse_brain/gene_sensitivity` | Fig 8 data |

(`...` = `$REPRO_DATA/processed/mouse_brain`; the exact paths are in `02_preprocess.sh`.)

---

## 8. Outputs

Figures are in `results/mouse_brain/figures/`:

| Figure | File |
|---|---|
| 1. Spatial maps by cell type | `spatial/spatial_grid_all.png` |
| 2. UMAPs | `umap/umap_self_grid_all.png` (+ one per sample) |
| 3. Cell-type composition | `composition/sorted_by_abundance/composition_grid_all_pct_linear_fixedy_tall.png` |
| 4. Cells retained vs. nCount floor | `threshold_sensitivity/celltype_threshold_sensitivity_tall.png` |
| 5. N genes / N cells / nCount / ARI | `summary_grid/summary_grid_shared_panel.png` |
| 6. Secondary-type weight (w2) | `w2/w2_boxplot_by_score.png` |
| 7. Spillover cosine similarity | `spillover/spillover_cosine.png` |
| 8. Per-gene sensitivity | `gene_sensitivity_scatter/scatter_self_vs_shared_panel_genes_transcripts_per_cell_median_log10.png` |

Each figure has a matching CSV of the plotted values next to it or in the stage's results folder.
Logs from every step are in `results/mouse_brain/logs/` (SLURM) or printed to the terminal.

---

## 9. Check these before trusting the figures

1. **Ingestion log (stage 1).** Check what % of StrataMap barcodes fell inside a segmented cell, and
   how many Trekker nuclei were dropped for having no coordinates.
2. **Fig 1.** Both sections should look like brain tissue with a sensible scale bar. A section
   that's 1000× too big or too small means `MB_TREKKER_MICRONS_PER_UNIT` is wrong.
3. **Fig 3.** Look for cell classes that can't be in your section (e.g. `29 CB Glut` cerebellar
   granule cells in a cortex section). If you see them, set `MB_REFERENCE_REGIONS`, delete
   `$REPRO_DATA/reference/mouse_brain/`, and rerun.
4. **Shared genes (stage 2b log).** The number of genes passing QC in both samples. If it's near
   zero, the gene symbols don't match between the two files (e.g. Ensembl IDs in one), so revisit
   `MB_TREKKER_GENE_SYMBOLS`.

---

## 10. Status

The pipeline hasn't been run end to end on real data yet. Expect the first run to hit problems,
most likely in the reference download/build, RCTD on the mouse reference, or memory limits. Save
the failing step's error output when that happens.
