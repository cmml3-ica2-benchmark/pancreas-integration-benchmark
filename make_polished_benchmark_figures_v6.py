#!/usr/bin/env python3
"""Create clean benchmark figures for the pancreas integration report.

Design choices in the plot-only update:
- Multi-panel UMAPs are drawn manually with matplotlib, not OmicVerse, so no per-panel
  legends, axis arrows, or duplicated labels are added automatically.
- Figure 2 keeps only the three main report metrics in panel a and uses fixed point sizes
  in the trade-off plot.
- Supplementary Figure S3 is redrawn from an existing unscaled scIB CSV only; no
  Benchmarker rerun is performed.
- Supplementary Figure S4 is redrawn from an existing flat summary CSV only and shows
  mean +/- 95% CI across the stored 20 stratified subsamples.
- Supplementary Table S1 is rendered reproducibly from the existing benchmark metrics CSV.
- If `results/benchmark_metrics.csv` is available, it is the source of truth for the main
  metric values; embedded fallback numbers are only for dry runs with no metrics CSV.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import anndata as ad
import matplotlib as mpl
import matplotlib.cbook as mpl_cbook
import numpy as np
import pandas as pd
from matplotlib.font_manager import fontManager
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

_cbook_impl_path = Path(mpl.__file__).with_name("cbook.py")
if _cbook_impl_path.exists():
    import importlib.util

    _spec = importlib.util.spec_from_file_location("matplotlib._cbook_impl", _cbook_impl_path)
    _cbook_impl = importlib.util.module_from_spec(_spec)
    assert _spec is not None and _spec.loader is not None
    _spec.loader.exec_module(_cbook_impl)
    for _name in ("_ExceptionInfo", "_is_pandas_dataframe"):
        if not hasattr(mpl_cbook, _name) and hasattr(_cbook_impl, _name):
            setattr(mpl_cbook, _name, getattr(_cbook_impl, _name))
if not hasattr(mpl_cbook, "_is_pandas_dataframe"):
    mpl_cbook._is_pandas_dataframe = lambda obj: isinstance(obj, pd.DataFrame)
import matplotlib.pyplot as plt

try:
    from adjustText import adjust_text
except Exception:  # optional
    adjust_text = None

METHOD_ORDER = ["Uncorrected", "Harmony", "Scanorama", "Seurat RPCA", "CONCORD"]
MAIN_METHODS = ["Uncorrected", "Harmony", "Scanorama", "Seurat RPCA"]

METHOD_COLORS = {
    "Uncorrected": "#7f7f7f",
    "Harmony": "#4c78a8",
    "Scanorama": "#f58518",
    "Seurat RPCA": "#54a24b",
    "CONCORD": "#b279a2",
}

UMAP_KEYS = {
    "Uncorrected": "X_umap_uncorrected",
    "Harmony": "X_umap_harmony",
    "Scanorama": "X_umap_scanorama",
    "Seurat RPCA": "X_umap_seurat_rpca",
    "CONCORD": "X_umap_concord",
}

FALLBACK_METRICS = pd.DataFrame(
    [
        ["Uncorrected", 0.325, 0.675, 0.076, 0.328, 0.662, 20, "euclidean"],
        ["Harmony", -0.047, 0.953, 0.207, 0.443, 0.662, 20, "euclidean"],
        ["Scanorama", -0.008, 0.992, 0.089, 0.333, 0.678, 20, "euclidean"],
        ["Seurat RPCA", 0.069, 0.931, 0.093, 0.405, 0.681, 20, "euclidean"],
        ["CONCORD", -0.063, 0.937, 0.169, 0.307, 0.662, 100, "cosine"],
    ],
    columns=[
        "Method",
        "Batch ASW",
        "Batch mixing score",
        "Cell-type ASW",
        "Leiden ARI",
        "Leiden NMI",
        "Dimensions",
        "Distance",
    ],
)

CORE_METRICS_MAIN = ["Batch mixing score", "Cell-type ASW", "Leiden ARI"]
CORE_METRICS_ALL = ["Batch mixing score", "Cell-type ASW", "Leiden ARI", "Leiden NMI"]

BATCH_PALETTE_FIXED = {
    "0": "#1f77b4",
    "1": "#d62728",
    "2": "#e377c2",
    "3": "#17becf",
    "Batch 0": "#1f77b4",
    "Batch 1": "#d62728",
    "Batch 2": "#e377c2",
    "Batch 3": "#17becf",
}

CELLTYPE_PALETTE_FIXED = {
    "alpha": "#1f77b4",
    "beta": "#aec7e8",
    "ductal": "#ff7f0e",
    "acinar": "#2ca02c",
    "delta": "#98df8a",
    "gamma": "#d62728",
    "endothelial": "#9467bd",
    "activated_stellate": "#c5b0d5",
    "activated stellate": "#c5b0d5",
    "quiescent_stellate": "#8c564b",
    "quiescent stellate": "#8c564b",
    "mesenchymal": "#e377c2",
    "macrophage": "#f7b6d2",
    "PSC": "#7f7f7f",
    "unclassified endocrine": "#bcbd22",
    "mast": "#dbdb8d",
    "epsilon": "#17becf",
    "mesenchyme": "#9edae5",
}


def set_plot_style() -> str:
    fonts = {f.name for f in fontManager.ttflist}
    font = "Arial" if "Arial" in fonts else "DejaVu Sans"
    mpl.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": [font, "Arial", "Helvetica", "DejaVu Sans", "Liberation Sans"],
            "font.size": 8.4,
            "axes.titlesize": 8.8,
            "axes.labelsize": 8.5,
            "xtick.labelsize": 7.8,
            "ytick.labelsize": 7.8,
            "legend.fontsize": 7.4,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "axes.edgecolor": "#cfcfcf",
            "axes.linewidth": 0.55,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )
    return font


def normalize_metrics_table(metrics: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "method": "Method",
        "batch_ASW_raw": "Batch ASW",
        "batch_mixing_score": "Batch mixing score",
        "celltype_ASW": "Cell-type ASW",
        "leiden_ARI": "Leiden ARI",
        "leiden_NMI": "Leiden NMI",
        "n_dimensions": "Dimensions",
        "metric": "Distance",
    }
    metrics = metrics.rename(columns={k: v for k, v in aliases.items() if k in metrics.columns}).copy()
    required = ["Method", "Batch ASW", "Batch mixing score", "Cell-type ASW", "Leiden ARI", "Leiden NMI", "Dimensions", "Distance"]
    missing = [c for c in required if c not in metrics.columns]
    if missing:
        raise ValueError(f"Metrics table is missing required columns: {missing}")
    metrics = metrics[required].copy()
    metrics["Method"] = pd.Categorical(metrics["Method"].astype(str), categories=METHOD_ORDER, ordered=True)
    metrics = metrics.sort_values("Method").reset_index(drop=True)
    for col in ["Batch ASW", "Batch mixing score", "Cell-type ASW", "Leiden ARI", "Leiden NMI"]:
        metrics[col] = pd.to_numeric(metrics[col], errors="coerce")
    metrics["Dimensions"] = pd.to_numeric(metrics["Dimensions"], errors="coerce").astype("Int64")
    metrics["Distance"] = metrics["Distance"].astype(str)
    return metrics


def load_metrics(metrics_csv: Optional[str]) -> pd.DataFrame:
    if metrics_csv:
        metrics_path = Path(metrics_csv)
        if not metrics_path.exists():
            raise FileNotFoundError(f"Metrics CSV not found: {metrics_path}")
        return normalize_metrics_table(pd.read_csv(metrics_path))
    return normalize_metrics_table(FALLBACK_METRICS.copy())


def resolve_obs_key(adata: ad.AnnData, requested: str, aliases: Sequence[str]) -> str:
    for key in [requested, *aliases]:
        if key in adata.obs.columns:
            return key
    raise KeyError(f"None of these obs columns were found: {[requested, *aliases]}")


def load_data(h5ad: Optional[str], metrics_csv: Optional[str], batch_key: str, celltype_key: str):
    metrics = load_metrics(metrics_csv)
    if not h5ad:
        return None, metrics, batch_key, celltype_key
    adata = ad.read_h5ad(h5ad)
    batch_key = resolve_obs_key(adata, batch_key, aliases=("batch", "Batch", "sample", "dataset"))
    celltype_key = resolve_obs_key(adata, celltype_key, aliases=("cell_type", "celltype", "CellType", "label", "labels"))
    return adata, metrics, batch_key, celltype_key


def export_main_metric_tables(metrics: pd.DataFrame, metrics_csv: Optional[str]) -> None:
    """Keep the Figure 2 summary CSVs synchronized with the main benchmark metrics table."""
    outdir = Path(metrics_csv).resolve().parent if metrics_csv else Path("results").resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    ordered = metrics.copy()
    ordered["Method"] = pd.Categorical(ordered["Method"].astype(str), categories=METHOD_ORDER, ordered=True)
    ordered = ordered.sort_values("Method").reset_index(drop=True)
    ordered["Method"] = ordered["Method"].astype(str)

    wide = ordered[["Method", "Batch mixing score", "Cell-type ASW", "Leiden ARI"]].copy()
    wide.to_csv(outdir / "Figure2_benchmark_main_wide.csv", index=False)

    long = wide.melt(id_vars="Method", var_name="Metric", value_name="Score")
    long["Method"] = pd.Categorical(long["Method"], categories=METHOD_ORDER, ordered=True)
    long = long.sort_values(["Metric", "Method"]).reset_index(drop=True)
    long["MethodColor"] = long["Method"].map(METHOD_COLORS)
    long["Method"] = long["Method"].astype(str)
    long.to_csv(outdir / "Figure2_benchmark_main_long.csv", index=False)


def validate_obsm_keys(adata: ad.AnnData, methods: Iterable[str]) -> None:
    missing = [m for m in methods if UMAP_KEYS[m] not in adata.obsm]
    if missing:
        available = list(adata.obsm.keys())
        raise KeyError(f"Missing UMAP keys: {{m: UMAP_KEYS[m] for m in missing}}. Available obsm keys: {available}")
    for m in methods:
        arr = np.asarray(adata.obsm[UMAP_KEYS[m]])
        if arr.ndim != 2 or arr.shape[1] != 2 or arr.shape[0] != adata.n_obs:
            raise ValueError(f"{UMAP_KEYS[m]} must have shape (n_obs, 2); got {arr.shape}")


def get_category_palette(categories: Sequence[str], mode: str) -> Dict[str, str]:
    cats = [str(c) for c in categories]
    if mode == "batch":
        return {cat: BATCH_PALETTE_FIXED.get(cat, "#1f77b4") for cat in cats}
    else:
        return {cat: CELLTYPE_PALETTE_FIXED.get(cat, "#999999") for cat in cats}


def _limits(coords: np.ndarray, pad_frac: float = 0.055) -> Tuple[Tuple[float, float], Tuple[float, float]]:
    x0, x1 = float(np.nanmin(coords[:, 0])), float(np.nanmax(coords[:, 0]))
    y0, y1 = float(np.nanmin(coords[:, 1])), float(np.nanmax(coords[:, 1]))
    xp = (x1 - x0) * pad_frac if x1 > x0 else 1
    yp = (y1 - y0) * pad_frac if y1 > y0 else 1
    return (x0 - xp, x1 + xp), (y0 - yp, y1 + yp)


def plot_embedding_panel(
    ax: plt.Axes,
    coords: np.ndarray,
    labels: pd.Series,
    palette: Mapping[str, str],
    point_size: float = 3.0,
    alpha: float = 0.88,
    xlim: Optional[Tuple[float, float]] = None,
    ylim: Optional[Tuple[float, float]] = None,
    show_spines: bool = False,
) -> None:
    """Manual UMAP plotting. No OmicVerse call, no automatic legend, no arrows."""
    labels = labels.astype(str)
    counts = labels.value_counts()
    # Draw large groups first, small groups later so rare groups stay visible.
    categories = list(counts.sort_values(ascending=False).index)
    for cat in categories:
        mask = labels.values == cat
        ax.scatter(
            coords[mask, 0],
            coords[mask, 1],
            s=point_size,
            c=palette.get(cat, "#999999"),
            alpha=alpha,
            linewidths=0,
            rasterized=True,
        )
    if xlim is None or ylim is None:
        xlim, ylim = _limits(coords)
        aspect_mode = "datalim"
    else:
        aspect_mode = "box"
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.set_aspect("equal", adjustable=aspect_mode)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel("")
    ax.set_ylabel("")
    for spine in ax.spines.values():
        spine.set_visible(show_spines)
        if show_spines:
            spine.set_color("#d6d6d6")
            spine.set_linewidth(0.55)


def save_figure(fig: plt.Figure, outdir: Path, stem: str, formats: Sequence[str] = ("png", "pdf")) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(outdir / f"{stem}.{ext}", bbox_inches="tight", dpi=300 if ext == "png" else None)


def legend_handles(palette: Mapping[str, str], marker_size: float = 4.2, replace_underscores: bool = True) -> List[Line2D]:
    out = []
    for label, color in palette.items():
        lab = str(label).replace("_", " ") if replace_underscores else str(label)
        out.append(Line2D([0], [0], marker="o", linestyle="", markerfacecolor=color, markeredgecolor="none", markersize=marker_size, label=lab))
    return out


def plot_reference_umap_panel(
    ax: plt.Axes,
    coords: np.ndarray,
    labels: Sequence[str] | pd.Series,
    categories: Sequence[str],
    palette: Mapping[str, str],
    title: Optional[str] = None,
    panel_label: Optional[str] = None,
    point_size: float = 2.0,
    alpha: float = 0.80,
    title_fontsize: float = 11.0,
    panel_fontsize: float = 11.5,
) -> None:
    values = pd.Series(labels).astype(str).values
    for category in categories:
        mask = values == str(category)
        if np.any(mask):
            ax.scatter(
                coords[mask, 0],
                coords[mask, 1],
                c=[palette[str(category)]],
                s=point_size,
                alpha=alpha,
                linewidths=0,
                rasterized=True,
            )
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)
    if title:
        ax.set_title(title, fontsize=title_fontsize, pad=4)
    if panel_label:
        ax.text(
            0.01,
            0.99,
            panel_label,
            transform=ax.transAxes,
            fontsize=panel_fontsize,
            fontweight="bold",
            va="top",
            ha="left",
        )


def make_figure_1(adata: ad.AnnData, outdir: Path, batch_key: str, celltype_key: str) -> None:
    validate_obsm_keys(adata, MAIN_METHODS)
    batch_cats = sorted(pd.unique(adata.obs[batch_key].astype(str)))
    celltype_cats = list(adata.obs[celltype_key].astype(str).value_counts().index)
    batch_palette = get_category_palette(batch_cats, "batch")
    celltype_palette = get_category_palette(celltype_cats, "celltype")
    coords_lookup = {method: np.asarray(adata.obsm[UMAP_KEYS[method]]) for method in MAIN_METHODS}
    fig = plt.figure(figsize=(12.8, 7.6), constrained_layout=True)
    gs = fig.add_gridspec(
        4,
        5,
        width_ratios=[0.38, 1, 1, 1, 1],
        height_ratios=[1.0, 0.20, 1.0, 0.36],
        wspace=0.05,
        hspace=0.03,
    )
    axes = np.empty((2, 4), dtype=object)
    ax_row_top = fig.add_subplot(gs[0, 0])
    ax_row_bottom = fig.add_subplot(gs[2, 0])
    ax_row_top.axis("off")
    ax_row_bottom.axis("off")

    for col, method in enumerate(MAIN_METHODS):
        coords = coords_lookup[method]
        axes[0, col] = fig.add_subplot(gs[0, col + 1])
        axes[1, col] = fig.add_subplot(gs[2, col + 1])
        plot_reference_umap_panel(
            axes[0, col],
            coords,
            adata.obs[batch_key],
            batch_cats,
            batch_palette,
            title=method,
            panel_label=chr(ord("A") + col),
            point_size=2.0,
            alpha=0.80,
            title_fontsize=14.0,
            panel_fontsize=12.8,
        )
        plot_reference_umap_panel(
            axes[1, col],
            coords,
            adata.obs[celltype_key],
            celltype_cats,
            celltype_palette,
            panel_label=chr(ord("E") + col),
            point_size=2.0,
            alpha=0.80,
            panel_fontsize=12.8,
        )

    ax_row_top.text(0.98, 0.50, "Colored by\nbatch", transform=ax_row_top.transAxes, fontsize=13.0, fontweight="bold", va="center", ha="right")
    ax_row_bottom.text(0.98, 0.50, "Colored by\ncell type", transform=ax_row_bottom.transAxes, fontsize=13.0, fontweight="bold", va="center", ha="right")

    ax_leg_top = fig.add_subplot(gs[1, 1:])
    ax_leg_top.axis("off")
    ax_leg_top.legend(
        handles=legend_handles(batch_palette, marker_size=6.8, replace_underscores=False),
        title="Batch",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.52),
        borderaxespad=0,
        handletextpad=0.56,
        labelspacing=0.62,
        columnspacing=1.34,
        ncol=4,
        title_fontsize=12.6,
        fontsize=11.6,
    )

    ax_leg_bottom = fig.add_subplot(gs[3, 1:])
    ax_leg_bottom.axis("off")
    ax_leg_bottom.legend(
        handles=legend_handles(celltype_palette, marker_size=6.2),
        title="Cell type",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.60),
        borderaxespad=0,
        handletextpad=0.52,
        labelspacing=0.56,
        columnspacing=1.12,
        ncol=4,
        title_fontsize=12.6,
        fontsize=10.8,
    )
    save_figure(fig, outdir, "Figure_1_visual_atlas_v6")
    plt.close(fig)


def _plot_lollipop_metric(
    ax: plt.Axes,
    metrics: pd.DataFrame,
    metric: str,
    show_ylabels: bool,
    x_floor: Optional[float] = None,
    method_order: Optional[Sequence[str]] = None,
) -> None:
    ordered = metrics.copy()
    if method_order is not None:
        ordered["Method"] = pd.Categorical(ordered["Method"].astype(str), categories=list(method_order), ordered=True)
        ordered = ordered.sort_values("Method").reset_index(drop=True)
        ordered["Method"] = ordered["Method"].astype(str)
    else:
        ordered = ordered.sort_values(metric, ascending=False).reset_index(drop=True)
    y = np.arange(len(ordered))
    values = ordered[metric].to_numpy(dtype=float)
    if x_floor is None:
        x_floor = min(0.0, float(np.nanmin(values)) - 0.02)
    x_max = float(np.nanmax(values)) + (0.045 if metric != "Batch mixing score" else 0.035)
    ax.set_xlim(x_floor, x_max)

    for i, row in ordered.iterrows():
        method = str(row["Method"])
        value = float(row[metric])
        ax.hlines(i, x_floor, value, color="#d9d9d9", lw=0.85, zorder=1)
        ax.scatter(value, i, s=32, color=METHOD_COLORS[method], edgecolor="white", lw=0.55, zorder=3)
        ax.text(value + (x_max - x_floor) * 0.024, i, f"{value:.3f}", va="center", ha="left", fontsize=7.2)

    ax.set_yticks(y)
    ax.set_yticklabels(ordered["Method"].astype(str) if show_ylabels else [])
    ax.invert_yaxis()
    ax.set_title(metric, pad=6, fontsize=8.4)
    ax.tick_params(axis="y", length=0, pad=3)
    ax.tick_params(axis="x", length=2.5, colors="#666666")
    ax.locator_params(axis="x", nbins=3)
    for spine in ["top", "right", "left"]:
        ax.spines[spine].set_visible(False)
    ax.spines["bottom"].set_color("#cfcfcf")
    ax.grid(False)


def make_figure_2(metrics: pd.DataFrame, outdir: Path) -> None:
    metrics = metrics.copy()
    metrics["Method"] = metrics["Method"].astype(str)
    metrics_main = metrics.loc[metrics["Method"].isin(MAIN_METHODS)].copy()

    fig = plt.figure(figsize=(7.35, 5.20))
    gs = fig.add_gridspec(3, 3, height_ratios=[1.0, 0.16, 1.34], hspace=0.28, wspace=0.34)
    metric_axes = [fig.add_subplot(gs[0, i]) for i in range(3)]
    ax_legend = fig.add_subplot(gs[1, :])
    trade_ax = fig.add_subplot(gs[1, :])
    trade_ax.remove()
    trade_ax = fig.add_subplot(gs[2, :])

    for i, metric in enumerate(CORE_METRICS_MAIN):
        floor = 0.60 if metric == "Batch mixing score" else 0.0
        _plot_lollipop_metric(
            metric_axes[i],
            metrics_main,
            metric,
            show_ylabels=(i == 0),
            x_floor=floor,
            method_order=MAIN_METHODS,
        )

    metric_axes[0].text(-0.20, 1.08, "a", transform=metric_axes[0].transAxes, fontsize=9.8, fontweight="bold", ha="left", va="bottom")

    ax_legend.axis("off")
    main_handles = [
        Line2D([0], [0], marker="o", linestyle="", markersize=5.6, markerfacecolor=METHOD_COLORS[m], markeredgecolor="none", label=m)
        for m in MAIN_METHODS
    ]
    ax_legend.legend(
        handles=main_handles,
        ncol=4,
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.40),
        columnspacing=1.3,
        handletextpad=0.45,
        borderaxespad=0,
        fontsize=8.2,
    )

    x = metrics["Batch mixing score"].to_numpy(float)
    y = metrics["Cell-type ASW"].to_numpy(float)
    methods = metrics["Method"].tolist()
    colors = [METHOD_COLORS[m] for m in methods]
    x_med, y_med = float(np.median(x)), float(np.median(y))
    xlim = (min(0.64, x.min() - 0.02), max(1.01, x.max() + 0.025))
    ylim = (min(0.055, y.min() - 0.015), y.max() + 0.035)

    trade_ax.add_patch(Rectangle((x_med, y_med), xlim[1] - x_med, ylim[1] - y_med, facecolor="#faf6ee", edgecolor="none", zorder=0))
    trade_ax.axvline(x_med, color="#d6d6d6", lw=0.8, zorder=1)
    trade_ax.axhline(y_med, color="#d6d6d6", lw=0.8, zorder=1)
    trade_ax.scatter(x, y, s=118, c=colors, alpha=0.95, edgecolors="white", linewidths=0.8, zorder=3)

    label_offsets = {
        "Uncorrected": (-0.011, -0.007),
        "Harmony": (0.010, 0.011),
        "Scanorama": (-0.020, -0.012),
        "Seurat RPCA": (-0.019, 0.010),
        "CONCORD": (0.010, -0.005),
    }
    texts = []
    for method, xv, yv in zip(methods, x, y):
        label = "CONCORD\n(exploratory)" if method == "CONCORD" else method
        dx, dy = label_offsets.get(method, (0.01, 0.01))
        texts.append(trade_ax.text(xv + dx, yv + dy, label, fontsize=7.5, ha="left", va="center", zorder=4))
    if adjust_text is not None:
        adjust_text(texts, ax=trade_ax, expand_points=(1.08, 1.18), force_text=0.22)
    trade_ax.set_xlim(xlim)
    trade_ax.set_ylim(ylim)
    trade_ax.set_xlabel("Batch mixing score")
    trade_ax.set_ylabel("Cell-type ASW")
    for spine in ["top", "right"]:
        trade_ax.spines[spine].set_visible(False)
    trade_ax.spines["left"].set_color("#cfcfcf")
    trade_ax.spines["bottom"].set_color("#cfcfcf")
    trade_ax.tick_params(length=2.5, colors="#555555")
    trade_ax.text(-0.055, 1.04, "b", transform=trade_ax.transAxes, fontsize=9.8, fontweight="bold", ha="left", va="bottom")

    save_figure(fig, outdir, "Figure_2_benchmark_synthesis_v6")
    plt.close(fig)


def make_supplementary_figure_s1(adata: ad.AnnData, outdir: Path, celltype_key: str) -> None:
    validate_obsm_keys(adata, MAIN_METHODS)
    cell_cats = list(adata.obs[celltype_key].astype(str).value_counts().index)
    cell_palette = get_category_palette(cell_cats, "celltype")
    coords_lookup = {method: np.asarray(adata.obsm[UMAP_KEYS[method]]) for method in MAIN_METHODS}
    fig = plt.figure(figsize=(7.45, 5.15))
    gs = fig.add_gridspec(2, 3, width_ratios=[1, 1, 1.02], hspace=0.14, wspace=0.14)
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
    ]
    ax_leg = fig.add_subplot(gs[:, 2])

    for ax, method in zip(axes, MAIN_METHODS):
        coords = coords_lookup[method]
        plot_reference_umap_panel(
            ax,
            coords,
            adata.obs[celltype_key],
            cell_cats,
            cell_palette,
            title=method,
            point_size=2.0,
            alpha=0.80,
        )
    axes[0].text(-0.17, 1.04, "a", transform=axes[0].transAxes, fontsize=10.6, fontweight="bold", ha="left", va="top")

    ax_leg.axis("off")
    ax_leg.legend(
        handles=legend_handles(cell_palette, marker_size=4.7),
        title="Cell type",
        ncol=1,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(0.0, 0.56),
        borderaxespad=0,
        handletextpad=0.44,
        columnspacing=0.82,
        labelspacing=0.56,
        fontsize=8.0,
        title_fontsize=8.8,
    )
    save_figure(fig, outdir, "Supplementary_Figure_S1_full_celltype_atlas_v6")
    plt.close(fig)


def make_supplementary_figure_s2(adata: ad.AnnData, metrics: pd.DataFrame, outdir: Path, batch_key: str, celltype_key: str) -> None:
    _ = metrics
    validate_obsm_keys(adata, ["CONCORD"])
    batch_cats = sorted(pd.unique(adata.obs[batch_key].astype(str)))
    cell_cats = list(adata.obs[celltype_key].astype(str).value_counts().index)
    batch_palette = get_category_palette(batch_cats, "batch")
    cell_palette = get_category_palette(cell_cats, "celltype")
    coords = np.asarray(adata.obsm[UMAP_KEYS["CONCORD"]])
    fig = plt.figure(figsize=(8.7, 4.85), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.38], hspace=0.02, wspace=0.10)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_leg_batch = fig.add_subplot(gs[1, 0])
    ax_leg_cell = fig.add_subplot(gs[1, 1])

    plot_reference_umap_panel(
        ax_a,
        coords,
        adata.obs[batch_key],
        batch_cats,
        batch_palette,
        title="Batch",
        panel_label="a",
        point_size=2.0,
        alpha=0.80,
        title_fontsize=13.0,
        panel_fontsize=13.0,
    )
    plot_reference_umap_panel(
        ax_b,
        coords,
        adata.obs[celltype_key],
        cell_cats,
        cell_palette,
        title="Cell type",
        panel_label="b",
        point_size=2.0,
        alpha=0.80,
        title_fontsize=13.0,
        panel_fontsize=13.0,
    )

    ax_leg_batch.axis("off")
    ax_leg_batch.legend(
        handles=legend_handles(batch_palette, marker_size=6.0, replace_underscores=False),
        title="Batch",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.60),
        borderaxespad=0,
        handletextpad=0.50,
        labelspacing=0.52,
        columnspacing=1.25,
        ncol=4,
        title_fontsize=11.0,
        fontsize=10.0,
    )

    ax_leg_cell.axis("off")
    ax_leg_cell.legend(
        handles=legend_handles(cell_palette, marker_size=5.6),
        title="Cell type",
        frameon=False,
        loc="center",
        bbox_to_anchor=(0.5, 0.58),
        borderaxespad=0,
        handletextpad=0.46,
        labelspacing=0.44,
        columnspacing=1.02,
        ncol=2,
        title_fontsize=11.0,
        fontsize=9.2,
    )
    save_figure(fig, outdir, "Supplementary_Figure_S2_CONCORD_v6")
    plt.close(fig)


def make_supplementary_table_s1(metrics: pd.DataFrame, outdir: Path) -> None:
    """Render the core metrics table reproducibly from the existing benchmark metrics CSV."""
    outdir.mkdir(parents=True, exist_ok=True)
    csv_path = outdir / "Supplementary_Table_S1_core_metrics_v6.csv"
    metrics.copy().to_csv(csv_path, index=False)

    fig, ax = plt.subplots(figsize=(7.1, 2.55))
    ax.axis("off")

    display = metrics.copy()
    for col in ["Batch ASW", "Batch mixing score", "Cell-type ASW", "Leiden ARI", "Leiden NMI"]:
        display[col] = display[col].map(lambda x: f"{float(x):.3f}")
    display["Dimensions"] = display["Dimensions"].astype(int).astype(str)
    display["Distance"] = display["Distance"].astype(str)

    table = ax.table(
        cellText=display.values,
        colLabels=display.columns,
        loc="center",
        cellLoc="center",
        colLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(7.5)
    table.scale(1.0, 1.40)

    best_cols = ["Batch mixing score", "Cell-type ASW", "Leiden ARI", "Leiden NMI"]
    best_rows = {col: int(metrics[col].astype(float).idxmax()) for col in best_cols}
    col_idx = {name: i for i, name in enumerate(display.columns)}

    for (r, c), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd")
        if r == 0:
            cell.set_facecolor("#efefef")
            cell.get_text().set_fontweight("bold")
            cell.set_linewidth(0.55)
        else:
            cell.set_facecolor("#fafafa" if r % 2 == 0 else "white")
            cell.set_linewidth(0.35)
            if c == 0:
                cell.get_text().set_ha("left")

    for col in best_cols:
        table[best_rows[col] + 1, col_idx[col]].get_text().set_fontweight("bold")

    ax.text(
        0.0,
        -0.12,
        "Core benchmark metrics from results/benchmark_metrics.csv. Higher values are better for batch mixing score, Cell-type ASW, Leiden ARI and Leiden NMI.",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=6.9,
        color="#333333",
    )
    save_figure(fig, outdir, "Supplementary_Table_S1_core_metrics_v6", formats=("png", "pdf"))
    plt.close(fig)


def make_dataset_composition_tables(adata: ad.AnnData, outdir: Path, batch_key: str, celltype_key: str) -> None:
    """Save descriptive n tables for report supplement."""
    outdir.mkdir(parents=True, exist_ok=True)
    obs = adata.obs[[batch_key, celltype_key]].astype(str).copy()
    batch_counts = obs[batch_key].value_counts().rename_axis("Batch").reset_index(name="n_cells")
    batch_counts["percent"] = batch_counts["n_cells"] / batch_counts["n_cells"].sum() * 100
    cell_counts = obs[celltype_key].value_counts().rename_axis("Cell type").reset_index(name="n_cells")
    cell_counts["percent"] = cell_counts["n_cells"] / cell_counts["n_cells"].sum() * 100
    cross = pd.crosstab(obs[celltype_key], obs[batch_key]).sort_index()
    batch_counts.to_csv(outdir / "Supplementary_Table_S2A_batch_composition_v6.csv", index=False)
    cell_counts.to_csv(outdir / "Supplementary_Table_S2B_celltype_composition_v6.csv", index=False)
    cross.to_csv(outdir / "Supplementary_Table_S2C_batch_by_celltype_counts_v6.csv")

    # A compact PDF overview with batch and cell-type counts; full cross-tab remains CSV.
    fig = plt.figure(figsize=(6.85, 3.45))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.78, 1.22], wspace=0.25)
    ax1, ax2 = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])
    for ax, df, title in [(ax1, batch_counts, "Batch composition"), (ax2, cell_counts.head(16), "Cell-type composition")]:
        ax.axis("off")
        temp = df.copy()
        temp["percent"] = temp["percent"].map(lambda x: f"{x:.1f}")
        tbl = ax.table(cellText=temp.values, colLabels=temp.columns, loc="center", cellLoc="center", colLoc="center")
        tbl.auto_set_font_size(False)
        tbl.set_fontsize(6.6)
        tbl.scale(1.0, 1.12)
        for (r, c), cell in tbl.get_celld().items():
            cell.set_edgecolor("#dddddd")
            cell.set_linewidth(0.35)
            if r == 0:
                cell.set_facecolor("#f0f0f0")
                cell.get_text().set_fontweight("bold")
            elif c == 0:
                cell.get_text().set_ha("left")
        ax.set_title(title, fontsize=8.3, pad=2)
    save_figure(fig, outdir, "Supplementary_Table_S2_dataset_composition_v6", formats=("png", "pdf"))
    plt.close(fig)


def _friendly_scib_name(name: str) -> str:
    mapping = {
        "X_uncorrected": "Uncorrected",
        "X_harmony": "Harmony",
        "X_scanorama": "Scanorama",
        "X_seurat_rpca": "Seurat RPCA",
        "X_concord": "CONCORD",
    }
    return mapping.get(name, name)


def _text_color_for_fill(fill_color, dark: str = "white", light: str = "#222222", threshold: float = 0.52) -> str:
    """Use white text on dark fills so labels stay readable in the scIB-style table."""
    r, g, b, _ = mpl.colors.to_rgba(fill_color)
    luminance = 0.2126 * r + 0.7152 * g + 0.0722 * b
    return dark if luminance < threshold else light


def make_supplementary_figure_s3_from_csv(scib_csv: str, outdir: Path) -> None:
    """Redraw a scaled scIB-style benchmark table from an existing CSV only."""
    csv_path = Path(scib_csv)
    if not csv_path.exists():
        raise FileNotFoundError(f"scIB CSV not found: {csv_path}")

    raw = pd.read_csv(csv_path)
    if raw.empty:
        raise ValueError(f"scIB CSV is empty: {csv_path}")

    metric_type_row = raw.loc[raw["Embedding"].astype(str) == "Metric Type"]
    metric_types = metric_type_row.iloc[0].to_dict() if not metric_type_row.empty else {}
    table_df = raw.loc[raw["Embedding"].astype(str) != "Metric Type"].copy()
    table_df["Embedding"] = table_df["Embedding"].map(_friendly_scib_name)

    numeric_cols = [c for c in table_df.columns if c != "Embedding"]
    for col in numeric_cols:
        table_df[col] = pd.to_numeric(table_df[col], errors="coerce")
    table_df = table_df.replace([np.inf, -np.inf], np.nan)

    preferred_order = ["Uncorrected", "Harmony", "Scanorama", "Seurat RPCA", "CONCORD"]
    table_df["Embedding"] = pd.Categorical(table_df["Embedding"], categories=preferred_order, ordered=True)
    table_df = table_df.sort_values("Embedding").reset_index(drop=True)
    table_df["Embedding"] = table_df["Embedding"].astype(str)

    metric_cols = [
        "Isolated labels",
        "KMeans NMI",
        "KMeans ARI",
        "Silhouette label",
        "cLISI",
        "BRAS",
        "iLISI",
        "KBET",
        "Graph connectivity",
        "PCR comparison",
    ]
    aggregate_cols = ["Batch correction", "Bio conservation", "Total"]

    fig = plt.figure(figsize=(12.2, 3.9))
    ax = fig.add_subplot(111)
    ax.set_xlim(0, 18.3)
    n_rows = len(table_df)
    ax.set_ylim(-0.7, n_rows + 1.35)
    ax.axis("off")

    left_margin = 2.25
    metric_step = 1.12
    metric_x = {name: left_margin + i * metric_step for i, name in enumerate(metric_cols)}
    agg_x = {
        "Batch correction": left_margin + len(metric_cols) * metric_step + 1.05,
        "Bio conservation": left_margin + len(metric_cols) * metric_step + 2.35,
        "Total": left_margin + len(metric_cols) * metric_step + 3.65,
    }

    header_band_y0 = n_rows + 0.58
    group_specs = [
        ("Bio conservation", metric_x["Isolated labels"] - 0.56, metric_x["cLISI"] + 0.56, "#dcecf8"),
        ("Batch correction", metric_x["BRAS"] - 0.56, metric_x["PCR comparison"] + 0.56, "#f8edd8"),
        ("Aggregate score", agg_x["Batch correction"] - 0.62, agg_x["Total"] + 0.62, "#eeeeee"),
    ]
    for label, x0, x1, color in group_specs:
        ax.add_patch(Rectangle((x0, header_band_y0), x1 - x0, 0.40, facecolor=color, edgecolor="none"))
        ax.text((x0 + x1) / 2, header_band_y0 + 0.20, label, ha="center", va="center", fontsize=8.6, color="#444444")

    header_y = n_rows + 0.10
    ax.text(0.15, header_y, "Method", ha="left", va="bottom", fontsize=8.8, fontweight="bold")
    for name in metric_cols:
        header_map = {
            "Isolated labels": "Isolated\nlabels",
            "KMeans NMI": "KMeans\nNMI",
            "KMeans ARI": "KMeans\nARI",
            "Silhouette label": "Silhouette\nlabel",
            "Graph connectivity": "Graph\nconnectivity",
            "PCR comparison": "PCR\ncomparison",
        }
        header = header_map.get(name, name)
        ax.text(metric_x[name], header_y, header, ha="center", va="bottom", fontsize=7.3, color="#333333", linespacing=0.95)
    for name in aggregate_cols:
        header_map = {
            "Batch correction": "Batch\ncorrection",
            "Bio conservation": "Bio\nconservation",
        }
        header = header_map.get(name, name)
        ax.text(agg_x[name], header_y, header, ha="center", va="bottom", fontsize=7.3, color="#333333", linespacing=0.95)

    cmap = plt.get_cmap("PRGn")
    agg_cmap = plt.get_cmap("YlGnBu")

    y_positions = np.arange(n_rows - 0.5, -0.5, -1.0)
    for row_idx, (_, row) in enumerate(table_df.iterrows()):
        y = y_positions[row_idx]
        method = str(row["Embedding"])
        ax.hlines(y - 0.5, xmin=0.0, xmax=18.1, color="#888888", lw=0.9, linestyles=(0, (1.2, 3.2)))
        ax.text(0.10, y, method, ha="left", va="center", fontsize=8.2, fontweight="bold" if method == "Harmony" else "normal", color="#222222")

        for name in metric_cols:
            value = float(row[name])
            x = metric_x[name]
            color = cmap(value)
            radius = 0.23 if metric_types.get(name, "") != "Aggregate score" else 0.20
            circ = plt.Circle((x, y), radius=radius, facecolor=color, edgecolor="none", alpha=0.95)
            ax.add_patch(circ)
            txt_color = _text_color_for_fill(color)
            ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=6.8, color=txt_color)

        for name in aggregate_cols:
            value = float(row[name])
            x = agg_x[name]
            color = agg_cmap(value)
            rect = Rectangle((x - 0.36, y - 0.18), 0.72, 0.36, facecolor=color, edgecolor="white", linewidth=0.4)
            ax.add_patch(rect)
            txt_color = _text_color_for_fill(color, light="#0f172a")
            ax.text(x, y, f"{value:.2f}", ha="center", va="center", fontsize=6.6, color=txt_color)

    ax.hlines(n_rows - 0.05, xmin=0.0, xmax=18.1, color="#777777", lw=1.0)
    ax.hlines(-0.02, xmin=0.0, xmax=18.1, color="#777777", lw=1.0)
    ax.vlines(agg_x["Batch correction"] - 0.68, ymin=-0.02, ymax=n_rows - 0.05, color="#777777", lw=0.8)

    outdir.mkdir(parents=True, exist_ok=True)
    out_csv = outdir / "Supplementary_Figure_S3_scib_results_table_scaled_v6.csv"
    if csv_path.resolve() != out_csv.resolve():
        shutil.copyfile(csv_path, out_csv)
    save_figure(fig, outdir, "Supplementary_Figure_S3_scib_results_table_scaled_v6", formats=("png", "pdf"))
    plt.close(fig)


def make_supplementary_figure_s4_from_flat(summary_flat_csv: str, outdir: Path) -> None:
    """Redraw the resampling sensitivity figure from an existing flat summary CSV only."""
    df = pd.read_csv(summary_flat_csv)
    required = ["Method"]
    for metric in CORE_METRICS_ALL:
        required.extend([f"{metric}__mean", f"{metric}__ci95"])
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Flat resampling summary is missing required columns: {missing}")

    df["Method"] = pd.Categorical(df["Method"].astype(str), categories=METHOD_ORDER, ordered=True)
    df = df.sort_values("Method").reset_index(drop=True)
    plot_metrics = CORE_METRICS_ALL

    fig, axes = plt.subplots(2, 2, figsize=(7.90, 5.95), sharey=True)
    axes = axes.ravel()
    panel_labels = ["a", "b", "c", "d"]
    for idx, (ax, metric) in enumerate(zip(axes, plot_metrics)):
        temp = df[["Method", f"{metric}__mean", f"{metric}__ci95"]].copy()
        temp.columns = ["Method", "mean", "ci95"]
        temp["Method"] = pd.Categorical(temp["Method"].astype(str), categories=METHOD_ORDER, ordered=True)
        temp = temp.sort_values("Method").reset_index(drop=True)
        temp["Method"] = temp["Method"].astype(str)
        y = np.arange(len(temp) - 1, -1, -1)
        x_floor = 0.60 if metric == "Batch mixing score" else 0.0
        x_max = float((temp["mean"] + temp["ci95"]).max()) + (0.03 if metric == "Batch mixing score" else 0.04)

        for i, row in temp.iterrows():
            yy = y[i]
            method = str(row["Method"])
            mean = float(row["mean"])
            ci95 = float(row["ci95"])
            ax.hlines(yy, x_floor, mean, color="#dddddd", lw=0.9, zorder=1)
            ax.errorbar(
                mean,
                yy,
                xerr=ci95,
                fmt="none",
                ecolor="#8f8f8f",
                elinewidth=0.85,
                capsize=2,
                zorder=2,
            )
            ax.scatter(mean, yy, s=31, color=METHOD_COLORS[method], edgecolor="white", linewidth=0.55, zorder=3)
            span = x_max - x_floor
            text_dx = span * 0.030
            if mean > x_max - span * 0.13:
                text_x = mean - text_dx
                text_ha = "right"
            else:
                text_x = mean + text_dx
                text_ha = "left"
            ax.text(
                text_x,
                yy,
                f"{mean:.3f}",
                va="center",
                ha=text_ha,
                fontsize=7.6,
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.72, pad=0.15),
            )

        ax.set_xlim(x_floor, x_max)
        ax.set_title(metric, fontsize=9.4, pad=6)
        ax.set_yticks(y)
        ax.set_yticklabels(temp["Method"].astype(str))
        ax.tick_params(axis="y", length=0, pad=3, labelsize=8.4)
        ax.tick_params(axis="x", length=2.5, colors="#666666", labelsize=8.2)
        ax.locator_params(axis="x", nbins=3)
        for spine in ["top", "right", "left"]:
            ax.spines[spine].set_visible(False)
        ax.spines["bottom"].set_color("#cfcfcf")
        ax.grid(False)
        ax.text(-0.16, 1.03, panel_labels[idx], transform=ax.transAxes, fontsize=10.0, fontweight="bold", ha="left", va="bottom")

    fig.text(0.5, 0.015, "Mean +/- 95% CI across 20 stratified subsamples", ha="center", va="bottom", fontsize=7.6, color="#333333")
    fig.tight_layout(rect=(0.0, 0.06, 1.0, 1.0))
    save_figure(fig, outdir, "Supplementary_Figure_S4_resampling_sensitivity_highdim", formats=("png", "pdf"))
    plt.close(fig)


def print_captions(made_umaps: bool) -> None:
    print("\nSuggested captions\n" + "=" * 18)
    if made_umaps:
        print("Figure 1. Visual atlas of pancreas integration results. UMAPs show batch structure (top row) and cell-type organization (bottom row) for the uncorrected data and the three main correction methods.")
        print("Supplementary Figure S1. Full cell-type atlas for the four main embeddings, with the complete legend shown separately.")
        print("Supplementary Figure S2. CONCORD shown separately as an exploratory embedding because it uses a 100-dimensional cosine representation rather than the 20-dimensional Euclidean embeddings used in the main comparison.")
        print("Supplementary Table S2. Dataset composition after filtering, including cell numbers per batch and per annotated cell type.")
    print("Figure 2. Benchmark synthesis across batch removal and biological conservation metrics. Lollipop panels show the three core report metrics, and the trade-off map summarises batch mixing against cell-type conservation.")
    print("Supplementary Table S1. Core benchmark metrics rendered reproducibly from the benchmark metrics CSV.")
    print("Supplementary Figure S3. Unscaled scIB-style benchmark table redrawn from the existing CSV only.")
    print("Supplementary Figure S4. Resampling sensitivity shown descriptively as mean +/- 95% CI across 20 stratified subsamples.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--h5ad", default="results/polished_figure_input.h5ad", help="Path to integrated AnnData with UMAP coordinates.")
    parser.add_argument("--metrics", default="results/benchmark_metrics.csv", help="Path to core metrics CSV.")
    parser.add_argument("--outdir", default="figures_polished_v6", help="Output directory.")
    parser.add_argument("--batch-key", default="batch", help="Batch column in adata.obs.")
    parser.add_argument("--celltype-key", default="cell_type", help="Cell-type column in adata.obs.")
    parser.add_argument(
        "--scib-unscaled-csv",
        default="figures_polished_v6/Supplementary_Figure_S3_scib_results_table_unscaled_v6.csv",
        help="Existing unscaled scIB results CSV used to redraw Supplementary Figure S3.",
    )
    parser.add_argument(
        "--resampling-summary-flat",
        default="results/highdim_resampling/resampling_metrics_summary_highdim_flat.csv",
        help="Existing flat summary CSV used to redraw Supplementary Figure S4.",
    )
    parser.add_argument(
        "--resampling-outdir",
        default="results/highdim_resampling",
        help="Output directory for the plot-only Supplementary Figure S4 render.",
    )
    args = parser.parse_args()

    set_plot_style()
    outdir = Path(args.outdir)
    adata, metrics, batch_key, celltype_key = load_data(args.h5ad, args.metrics, args.batch_key, args.celltype_key)
    export_main_metric_tables(metrics, args.metrics)

    make_figure_2(metrics, outdir)
    make_supplementary_table_s1(metrics, outdir)
    if args.scib_unscaled_csv and Path(args.scib_unscaled_csv).exists():
        make_supplementary_figure_s3_from_csv(args.scib_unscaled_csv, outdir)
    else:
        print("No existing unscaled scIB CSV provided; skipped Supplementary Figure S3 redraw.")
    if args.resampling_summary_flat and Path(args.resampling_summary_flat).exists():
        make_supplementary_figure_s4_from_flat(args.resampling_summary_flat, Path(args.resampling_outdir))
    else:
        print("No existing flat resampling summary provided; skipped Supplementary Figure S4 redraw.")

    made_umaps = False
    if adata is not None:
        make_figure_1(adata, outdir, batch_key, celltype_key)
        make_supplementary_figure_s1(adata, outdir, celltype_key)
        make_supplementary_figure_s2(adata, metrics, outdir, batch_key, celltype_key)
        make_dataset_composition_tables(adata, outdir, batch_key, celltype_key)
        made_umaps = True
    else:
        print("No AnnData file provided; UMAP and dataset-composition outputs were skipped.")

    print_captions(made_umaps)
    print(f"\nOutputs written to: {outdir.resolve()}")


if __name__ == "__main__":
    main()
