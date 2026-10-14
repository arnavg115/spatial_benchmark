"""
Builds a mouse brain RCTD reference from the Allen Brain Cell Atlas (ABC atlas) whole mouse brain
10x scRNA-seq taxonomy - Yao et al. 2023, Nature 624:317-332, "A high-resolution transcriptomic and
spatial atlas of cell types in the whole mouse brain" (https://doi.org/10.1038/s41586-023-06812-z).
Public, no account needed; file layout per the atlas's own manifest
(https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/releases/20241130/manifest.json):

  - metadata/WMB-10X/20241115/views/cell_metadata_with_cluster_annotation.csv (~1.4 GB) - one row
    per cell: cell_label, feature_matrix_label, library_method, class, class_color, ...
  - metadata/WMB-10X/20241115/gene.csv - gene_identifier (Ensembl) -> gene_symbol
  - expression_matrices/WMB-10Xv3/20230630/<feature_matrix_label>-raw.h5ad - raw UMI counts, one
    file per dissection region (13 files, ~95 GB total for all of them)

Taxonomy levels used (RCTD needs a Level2 annotation + a Level1 class grouping, see
run_rctd_doublet.py):
  Level2 = ABC "class" (34 classes, e.g. "01 IT-ET Glut", "31 OPC-Oligo"). Subclass (~338 types) is
           too fine for RCTD doublet mode on a sample of unknown region, and most subclasses are
           region-restricted.
  Level1 = neurotransmitter suffix of the class name for neuronal classes (Glut/GABA/Dopa/Sero),
           Glia/Vascular/Immune for the non-neuronal ones (classes 30-34).

Cells are sampled per class (--max-cells-per-type, default 2000) across every requested region's
matrix, so the reference covers the whole brain without loading ~4M cells. Restrict --regions once
the tissue's region is known: fewer downloads, and no off-region types for RCTD to pick from.

Also writes the class -> ABC class_color palette (CSV: cell_type,color) that the mouse brain
benchmark's figures use (benchmarks/mouse_brain/registry.R's CELLTYPE_PALETTE_PATH).

Usage: python build_reference_abc_mouse.py --out <ref.h5ad> --palette-out <palette.csv> \
           --download-dir <dir> [--regions all|HPF,Isocortex-1,...] [--max-cells-per-type 2000]
"""
import argparse
import urllib.request
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

BASE_URL = "https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com"
CELL_METADATA = "metadata/WMB-10X/20241115/views/cell_metadata_with_cluster_annotation.csv"
GENE_METADATA = "metadata/WMB-10X/20241115/gene.csv"
MATRIX_TEMPLATE = "expression_matrices/WMB-10Xv3/20230630/{label}-raw.h5ad"
REGIONS = ["CB", "CTXsp", "HPF", "HY", "Isocortex-1", "Isocortex-2", "MB", "MY", "OLF", "P", "PAL", "STR", "TH"]

NEURONAL_SUFFIXES = {"Glut", "GABA", "Dopa", "Sero"}
NON_NEURONAL_LEVEL1 = {
    "30 Astro-Epen": "Glia", "31 OPC-Oligo": "Glia", "32 OEC": "Glia",
    "33 Vascular": "Vascular", "34 Immune": "Immune",
}


def download(rel_path, download_dir):
    dest = Path(download_dir) / rel_path
    if dest.exists():
        print(f"Already downloaded: {dest}")
        return dest
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Download to a .part file first, so an interrupted multi-GB download is never mistaken for a
    # complete one on the next run.
    part = dest.with_suffix(dest.suffix + ".part")
    print(f"Downloading {BASE_URL}/{rel_path} -> {dest}")
    urllib.request.urlretrieve(f"{BASE_URL}/{rel_path}", part)
    part.rename(dest)
    return dest


def level1_for_class(cls):
    if cls in NON_NEURONAL_LEVEL1:
        return NON_NEURONAL_LEVEL1[cls]
    suffix = cls.split()[-1]
    if suffix in NEURONAL_SUFFIXES:
        return suffix
    raise SystemExit(f"Class {cls!r} has no Level1 mapping - the ABC taxonomy may have changed; update "
                     "NON_NEURONAL_LEVEL1/NEURONAL_SUFFIXES in this script")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True, help="output reference .h5ad")
    ap.add_argument("--palette-out", required=True, help="output class -> color CSV")
    ap.add_argument("--download-dir", required=True, help="where ABC atlas files are cached (~95 GB for all regions)")
    ap.add_argument("--regions", default="all", help=f"'all' or a comma-separated subset of {','.join(REGIONS)}")
    ap.add_argument("--max-cells-per-type", type=int, default=2000)
    ap.add_argument("--min-cells-per-type", type=int, default=25, help="drop classes with fewer sampled cells (RCTD needs enough per type)")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    regions = REGIONS if args.regions == "all" else args.regions.split(",")
    unknown = set(regions) - set(REGIONS)
    if unknown:
        raise SystemExit(f"Unknown --regions {sorted(unknown)} - choose from {REGIONS}")
    matrix_labels = [f"WMB-10Xv3-{r}" for r in regions]

    print("Loading ABC cell metadata")
    meta = pd.read_csv(download(CELL_METADATA, args.download_dir),
                       usecols=["cell_label", "feature_matrix_label", "library_method", "class", "class_color"],
                       dtype=str)
    meta = meta[(meta["library_method"] == "10Xv3") & meta["feature_matrix_label"].isin(matrix_labels)
                & meta["class"].notna()]
    print(f"{len(meta)} annotated 10Xv3 cells in regions {regions}")

    sampled = (meta.groupby("class", group_keys=False)
               .apply(lambda g: g.sample(n=min(len(g), args.max_cells_per_type), random_state=args.seed)))
    print(f"Sampled {len(sampled)} cells across {sampled['class'].nunique()} classes "
          f"(max {args.max_cells_per_type} per class)")

    pieces = []
    for label, cells in sampled.groupby("feature_matrix_label"):
        path = download(MATRIX_TEMPLATE.format(label=label), args.download_dir)
        backed = ad.read_h5ad(path, backed="r")
        idx = np.sort(backed.obs_names.get_indexer(cells["cell_label"]))
        if (idx < 0).any():
            raise SystemExit(f"{(idx < 0).sum()} sampled cell_labels not found in {path}'s obs_names")
        piece = backed[idx].to_memory()
        backed.file.close()
        print(f"  {label}: {piece.n_obs} cells")
        pieces.append(piece)
    adata = ad.concat(pieces, join="inner")

    genes = pd.read_csv(download(GENE_METADATA, args.download_dir)).set_index("gene_identifier")
    missing = ~adata.var_names.isin(genes.index)
    if missing.any():
        raise SystemExit(f"{missing.sum()} var_names not in gene.csv - expected Ensembl IDs as var_names")
    adata.var["gene_identifier"] = adata.var_names
    adata.var_names = genes.loc[adata.var_names, "gene_symbol"].astype(str).to_numpy()
    adata.var_names_make_unique()

    obs = sampled.set_index("cell_label").loc[adata.obs_names]
    adata.obs = pd.DataFrame({"cell_type": obs["class"].to_numpy(), "Level2": obs["class"].to_numpy(),
                              "Level1": [level1_for_class(c) for c in obs["class"]],
                              "feature_matrix_label": obs["feature_matrix_label"].to_numpy()},
                             index=adata.obs_names)

    counts = adata.obs["Level2"].value_counts()
    too_few = counts[counts < args.min_cells_per_type].index
    if len(too_few):
        print(f"Dropping {len(too_few)} classes with < {args.min_cells_per_type} cells: {list(too_few)}")
        adata = adata[~adata.obs["Level2"].isin(too_few)].copy()

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.out)
    print(f"Saved: {args.out} ({adata.n_obs} cells x {adata.n_vars} genes)")
    print(adata.obs.groupby(["Level1", "Level2"], observed=True).size().to_string())

    palette = (sampled[["class", "class_color"]].drop_duplicates()
               .rename(columns={"class": "cell_type", "class_color": "color"})
               .sort_values("cell_type"))
    palette = palette[palette["cell_type"].isin(adata.obs["Level2"].unique())]
    palette.to_csv(args.palette_out, index=False)
    print(f"Saved: {args.palette_out} ({len(palette)} classes)")


if __name__ == "__main__":
    main()
