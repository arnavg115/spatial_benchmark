"""
Convert a Trekker (Takara; commercial Slide-tags) cell-level .h5ad into the 10x-style mtx directory
+ cell_metadata.parquet shape qc_seurat.R's --input-type mtxdir already reads (same shape as Xenium
5k's/VisiumHD's converted output).

Trekker h5ad layouts vary by pipeline version and by whatever processing happened afterwards, so
nothing about where the counts, gene symbols or coordinates live is assumed - each is an explicit
flag. Run with --inspect first to see what the file actually contains:

    python convert_trekker_h5ad.py --h5ad trekker.h5ad --inspect

Checks that fail loudly instead of silently producing wrong-scale or wrong-kind data:
  - counts must be non-negative integers (normalized/log values break the nCount QC floor, RCTD and
    every sensitivity metric) - a non-integer matrix is rejected, not rounded;
  - --microns-per-unit is required: spillover.R's 15um neighbor radius and gene_sensitivity.R's
    per-mm^2 tissue area both assume microns (the VisiumHD converter's pixel-vs-micron bug is the
    precedent - see convert_visiumhd_segmented_cells.py's header);
  - nuclei without spatial coordinates (Trekker can't position every nucleus) are dropped and
    counted in the log, not given fake coordinates.

Usage:
    python convert_trekker_h5ad.py --h5ad <file.h5ad> --out-dir <dir> \
        --counts X|layer:<name>|raw --gene-symbols var_names|var:<column> \
        --coords obsm:<key>|obs:<x_column>,<y_column> --microns-per-unit <float>
"""
import argparse
import gzip
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse


def summarize_matrix(name, m, n_rows=2000):
    sub = m[: min(n_rows, m.shape[0])]
    vals = sub.data if scipy.sparse.issparse(sub) else np.asarray(sub).ravel()
    vals = vals[vals != 0]
    if len(vals) == 0:
        print(f"  {name}: {m.shape}, dtype={m.dtype}, first {n_rows} rows all zero")
        return
    is_int = np.all(np.equal(np.mod(vals, 1), 0))
    print(f"  {name}: {m.shape}, dtype={m.dtype}, sparse={scipy.sparse.issparse(m)}, "
          f"nonzero min={vals.min():.4g} max={vals.max():.4g}, integer-valued={is_int}, "
          f"example values={np.unique(vals)[:8].tolist()}")


def inspect(adata):
    print(f"Shape: {adata.n_obs} cells x {adata.n_vars} genes")
    print("Count matrix candidates (first 2000 cells sampled):")
    summarize_matrix("X", adata.X)
    for k in adata.layers.keys():
        summarize_matrix(f"layer:{k}", adata.layers[k])
    if adata.raw is not None:
        summarize_matrix("raw", adata.raw.X)
    print(f"var_names (first 5): {adata.var_names[:5].tolist()}")
    print(f"var columns: {adata.var.columns.tolist()}")
    for c in adata.var.columns:
        print(f"  var:{c} (first 5): {adata.var[c].astype(str).iloc[:5].tolist()}")
    print(f"obs columns: {adata.obs.columns.tolist()}")
    print("obsm entries:")
    for k in adata.obsm.keys():
        v = np.asarray(adata.obsm[k])
        if v.ndim == 2 and v.shape[1] >= 2 and np.issubdtype(v.dtype, np.number):
            print(f"  obsm:{k}: shape={v.shape}, col0 range=[{np.nanmin(v[:, 0]):.4g}, {np.nanmax(v[:, 0]):.4g}], "
                  f"col1 range=[{np.nanmin(v[:, 1]):.4g}, {np.nanmax(v[:, 1]):.4g}], NaN rows={np.isnan(v[:, :2]).any(axis=1).sum()}")
        else:
            print(f"  obsm:{k}: shape={v.shape}, dtype={v.dtype}")
    print("Numeric obs columns (possible coordinates):")
    for c in adata.obs.columns:
        if pd.api.types.is_numeric_dtype(adata.obs[c]):
            v = adata.obs[c].to_numpy(dtype=float)
            print(f"  obs:{c}: range=[{np.nanmin(v):.4g}, {np.nanmax(v):.4g}], NaN={np.isnan(v).sum()}")


def get_counts(adata, spec):
    if spec == "X":
        m, var_names = adata.X, adata.var_names
    elif spec == "raw":
        if adata.raw is None:
            raise SystemExit("--counts raw: this h5ad has no .raw")
        m, var_names = adata.raw.X, adata.raw.var_names
    elif spec.startswith("layer:"):
        m, var_names = adata.layers[spec.split(":", 1)[1]], adata.var_names
    else:
        raise SystemExit(f"--counts must be X, raw or layer:<name>, got {spec!r}")
    m = scipy.sparse.csr_matrix(m)
    if m.nnz and m.data.min() < 0:
        raise SystemExit(f"--counts {spec}: negative values present - this is not a raw count matrix")
    if m.nnz and not np.all(np.equal(np.mod(m.data, 1), 0)):
        raise SystemExit(f"--counts {spec}: non-integer values present (normalized/log-transformed?) - "
                         "point --counts at the raw UMI counts (see --inspect)")
    m.data = m.data.astype(np.int64)
    return m, var_names


def get_symbols(adata, spec, counts_spec):
    var = adata.raw.var if counts_spec == "raw" else adata.var
    if spec == "var_names":
        return pd.Index(var.index).astype(str)
    if spec.startswith("var:"):
        return pd.Index(var[spec.split(":", 1)[1]]).astype(str)
    raise SystemExit(f"--gene-symbols must be var_names or var:<column>, got {spec!r}")


def get_coords(adata, spec):
    if spec.startswith("obsm:"):
        xy = np.asarray(adata.obsm[spec.split(":", 1)[1]])[:, :2].astype(float)
    elif spec.startswith("obs:"):
        cols = spec.split(":", 1)[1].split(",")
        if len(cols) != 2:
            raise SystemExit("--coords obs: needs exactly two columns, e.g. obs:x,y")
        xy = adata.obs[cols].to_numpy(dtype=float)
    else:
        raise SystemExit(f"--coords must be obsm:<key> or obs:<x>,<y>, got {spec!r}")
    return xy


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--h5ad", required=True)
    ap.add_argument("--inspect", action="store_true", help="print what the file contains and exit")
    ap.add_argument("--out-dir")
    ap.add_argument("--counts", default="X", help="X | raw | layer:<name> - must hold raw integer UMI counts")
    ap.add_argument("--gene-symbols", default="var_names", help="var_names | var:<column>")
    ap.add_argument("--coords", help="obsm:<key> (first 2 columns = x, y) | obs:<x_column>,<y_column>")
    ap.add_argument("--microns-per-unit", type=float, help="multiply coordinates by this to get microns (1 if already microns)")
    args = ap.parse_args()

    print(f"Loading {args.h5ad}")
    adata = ad.read_h5ad(args.h5ad)

    if args.inspect:
        inspect(adata)
        return
    for flag in ("out_dir", "coords", "microns_per_unit"):
        if getattr(args, flag) is None:
            raise SystemExit(f"--{flag.replace('_', '-')} is required (run with --inspect first to choose it)")

    counts, var_names = get_counts(adata, args.counts)
    symbols = get_symbols(adata, args.gene_symbols, args.counts)
    xy = get_coords(adata, args.coords) * args.microns_per_unit

    positioned = ~np.isnan(xy).any(axis=1)
    print(f"{positioned.sum()} / {adata.n_obs} cells have spatial coordinates - dropping {(~positioned).sum()} without")
    counts = counts[positioned]
    xy = xy[positioned]
    # dash, not underscore: Seurat::CreateSeuratObject silently replaces "_" with "-" in cell names.
    cell_ids = pd.Index(adata.obs_names[positioned]).astype(str).str.replace("_", "-", regex=False)
    if cell_ids.has_duplicates:
        raise SystemExit("Duplicate cell IDs after '_' -> '-' sanitization")
    print(f"Coordinates (um): x=[{xy[:, 0].min():.1f}, {xy[:, 0].max():.1f}], y=[{xy[:, 1].min():.1f}, {xy[:, 1].max():.1f}]")
    print(f"nCount per cell: median={np.median(np.asarray(counts.sum(axis=1)).ravel()):.0f}")

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    with gzip.open(out_dir / "barcodes.tsv.gz", "wt") as f:
        f.write("\n".join(cell_ids) + "\n")
    features = pd.DataFrame({"id": pd.Index(var_names).astype(str), "symbol": symbols, "type": "Gene Expression"})
    features.to_csv(out_dir / "features.tsv.gz", sep="\t", header=False, index=False, compression="gzip")
    mtx_path = out_dir / "matrix.mtx"
    scipy.io.mmwrite(str(mtx_path), counts.T.tocoo(), field="integer")
    with open(mtx_path, "rb") as f_in, gzip.open(f"{mtx_path}.gz", "wb") as f_out:
        f_out.writelines(f_in)
    mtx_path.unlink()

    pd.DataFrame({"cell_id": cell_ids, "x_centroid": xy[:, 0], "y_centroid": xy[:, 1]}).to_parquet(
        out_dir / "cell_metadata.parquet", index=False)
    print(f"Done. Wrote {len(cell_ids)} cells x {len(symbols)} genes to {out_dir}")


if __name__ == "__main__":
    main()
