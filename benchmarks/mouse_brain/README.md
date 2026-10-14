# Mouse brain benchmark: StrataMap vs. Trekker

A second, separate benchmark built on the same analysis code as the breast comparison: one
StrataMap and one Trekker (Takara; commercial Slide-tags) mouse brain sample from the same animal,
both annotated against the same Allen Brain Cell Atlas reference, producing the same 8 figures.

Only three things differ from the breast benchmark:

| | Breast benchmark | Mouse brain benchmark |
|---|---|---|
| Samples | 8 (6 technologies) | 2 (`mouse_stratamap`, `mouse_trekker`) |
| RCTD reference | Janesick 2023 Chromium Flex (human breast) | ABC atlas whole mouse brain 10Xv3 (Yao et al. 2023, *Nature* 624:317), Level2 = ABC class |
| "Shared panel genes" axis (figs 5, 8) | literal 353-gene list | every gene passing QC in **both** samples (computed in STAGE 2b) |

`code/R/plot_common.R` loads `registry.R` from this directory when `BENCHMARK=mouse_brain`
(`config.sh` exports it), which swaps in the sample list, labels, colors, paths and gene lists. Every
other script is shared unchanged, so the QC thresholds, RCTD settings and metric definitions are
identical to the breast benchmark's.

## Running it

See **[RUNNING.md](RUNNING.md)** for where to put the data, the Trekker inspection step, every
setting, and the exact script and arguments behind each stage. In short: put the StrataMap folder
and Trekker h5ad under `$REPRO_DATA/raw/mouse_brain/`, set the four `MB_TREKKER_*` variables, then run
`benchmarks/mouse_brain/02_preprocess.sh` and `03_generate_figures.sh`. Results go to
`results/mouse_brain/` and never touch the breast benchmark's outputs.

## Interpreting the results

- **n = 1 per platform.** Differences describe these two samples, not the platforms in general.
  Section-to-section differences (region, depth, tissue quality) are confounded with platform,
  even within one animal. If the two sections come from different brain regions, the composition
  (fig 3) and ARI (fig 5) comparisons mostly reflect anatomy.
- **Nuclei vs. segmented cells.** Trekker profiles dissociated nuclei, while StrataMap aggregates
  transcripts within 5 µm-expanded cell segmentations. Lower counts per cell and a nuclear-transcript
  bias are expected for Trekker (figs 4, 5, 8).
- **Positioned nuclei only.** Trekker assigns coordinates to only some of the nuclei it profiles.
  Cell types that position poorly will be under-represented (fig 3).
- **Spillover (figs 6, 7) means different things on the two platforms.** For StrataMap, secondary
  RCTD weight that tracks the neighborhood suggests segmentation or diffusion spillover. Trekker
  nuclei are dissociated before sequencing, so a neighborhood correlation there more likely
  reflects ambient RNA, doublets, or genuine spatial co-localization of cell types.
- **Reference.** Both samples are annotated against ABC 10Xv3 single-cell profiles at the class level
  (34 classes). Until the region is known and `MB_REFERENCE_REGIONS` is restricted, RCTD can
  assign classes from regions that aren't in the section. Check fig 3 for implausible classes.
- **Shared gene axis.** Both platforms are whole-transcriptome, so "shared panel genes" means
  genes detected (≥ 10 cells after QC) in both. That is typically most expressed genes, so this
  axis differs much less from "self" than it does in the breast benchmark.

## Status

The Trekker converter was tested on a synthetic h5ad: counts round-trip exactly through Seurat's
`Read10X`, and normalized matrices and missing flags are rejected. The ABC file layout (obs =
`cell_label`, var = Ensembl IDs, raw integer counts) and the class → Level1 mapping were checked
against the live atlas files. The full pipeline has **not** been run end to end on real data.
