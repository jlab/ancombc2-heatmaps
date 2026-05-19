from __future__ import annotations

import os
from typing import Sequence

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .config import Subset
from .data import ANCOMData


def _format_percent(x: float) -> str:
    if pd.isna(x):
        return ""

    pct = x * 100
    if pct == 0:
        return "0%"
    if pct >= 10:
        return f"{pct:.0f}%"
    if pct >= 0.01:
        return f"{pct:.2f}%"
    return "<0.01%"


def _format_lfc(x: float, decimals: int) -> str:
    if pd.isna(x):
        return ""
    return f"{x:.{decimals}f}"


def _unique_labels(labels: Sequence[str]) -> list[str]:
    counts = {}
    out = []

    for label in labels:
        if label not in counts:
            counts[label] = 1
            out.append(label)
        else:
            counts[label] += 1
            out.append(f"{label} ({counts[label]})")

    return out


def _is_highlight(label: str, patterns: Sequence[str]) -> bool:
    lower = label.lower()
    return any(str(p).lower() in lower for p in patterns)


class HeatmapPlotter:
    def __init__(self, data: ANCOMData):
        self.data = data
        self.config = data.config

    def plot(
        self,
        subset: Subset,
        save_png: bool | None = None,
        save_pdf: bool | None = None,
        show: bool | None = None,
    ):
        cfg = self.config
        plot_cfg = cfg.plot

        if save_png is None:
            save_png = plot_cfg.save_png
        if save_pdf is None:
            save_pdf = plot_cfg.save_pdf
        if show is None:
            show = plot_cfg.show

        mean_rel, group_counts = self.data.mean_relative_abundance(subset)
        if mean_rel.empty:
            print(f"[NO HEATMAP] No relative abundance data: {subset.label}")
            return None

        lfc, sig, effect_col = self.data.ancom_matrices(subset)
        if lfc.empty:
            print(f"[NO HEATMAP] No ANCOM LFC data: {subset.label}")
            return None

        present_tps = [tp for tp in cfg.timepoints if tp in lfc.columns]
        lfc = lfc[present_tps]
        sig = sig[present_tps].astype("boolean").fillna(False).astype(bool)

        keep = sig.sum(axis=1) >= cfg.min_sig_cells_per_taxon
        lfc = lfc.loc[keep].copy()
        sig = sig.loc[keep].copy()

        if lfc.empty:
            print(f"[NO HEATMAP] No taxa after significance filter: {subset.label}")
            return None

        common = lfc.index.intersection(sig.index).intersection(mean_rel.index)
        lfc = lfc.loc[common]
        sig = sig.loc[common]
        mean_rel = mean_rel.loc[common]

        plot_df = lfc.mask(~sig)

        if cfg.remove_empty_rows_after_masking:
            plot_df = plot_df.loc[plot_df.notna().any(axis=1)]

        if plot_df.empty:
            print(f"[NO HEATMAP] No plottable significant cells: {subset.label}")
            return None

        mean_rel = mean_rel.reindex(plot_df.index)
        sig = sig.reindex(plot_df.index).fillna(False)
        lfc = lfc.reindex(plot_df.index)

        order = pd.DataFrame({
            "n_sig": sig.sum(axis=1),
            "mean_abs_lfc": lfc.abs().mean(axis=1),
        }).sort_values(["n_sig", "mean_abs_lfc"], ascending=[False, False]).index

        plot_df = plot_df.loc[order]
        mean_rel = mean_rel.loc[order]
        sig = sig.loc[order]
        lfc = lfc.loc[order]

        values = plot_df.to_numpy(dtype=float)
        if values.size == 0 or np.isnan(values).all():
            print(f"[NO HEATMAP] No numeric values: {subset.label}")
            return None

        max_abs = np.nanmax(np.abs(values))
        if not np.isfinite(max_abs) or max_abs == 0:
            max_abs = 1.0

        n_rows, n_cols = plot_df.shape

        row_sig_ra_mean = mean_rel.where(sig).mean(axis=1, skipna=True).fillna(0)
        scaled = (row_sig_ra_mean / max(row_sig_ra_mean.max(), 0.0001)).clip(0, 1)
        row_heights = 0.25 * (1 + 5 * scaled.to_numpy())

        y_edges = np.concatenate([[0], np.cumsum(row_heights)])
        y_centers = (y_edges[:-1] + y_edges[1:]) / 2
        x_edges = np.arange(n_cols + 1)
        x_centers = np.arange(n_cols) + 0.5

        fig_h = min(max(row_heights.sum() + 2.2, 4), 100)
        fig_w = min(max(n_cols * 0.8 + 7.5, 9), 16)

        fig, ax = plt.subplots(figsize=(fig_w, fig_h))

        cmap = mpl.colormaps[plot_cfg.heatmap_cmap].copy()
        cmap.set_bad("white")

        mesh = ax.pcolormesh(
            x_edges,
            y_edges,
            np.ma.masked_invalid(values),
            cmap=cmap,
            vmin=-max_abs,
            vmax=max_abs,
            shading="flat",
            edgecolors="lightgray",
            linewidth=0.35,
        )

        ax.set_ylim(y_edges[-1], 0)

        if plot_cfg.cell_text_mode != "none":
            white_threshold = 0.72 * max_abs

            for i, taxon in enumerate(plot_df.index):
                for j, tp in enumerate(plot_df.columns):
                    value = plot_df.loc[taxon, tp]
                    if pd.isna(value):
                        continue

                    if plot_cfg.cell_text_mode == "relative_abundance":
                        text = _format_percent(mean_rel.loc[taxon, tp])
                    elif plot_cfg.cell_text_mode == "lfc":
                        text = _format_lfc(value, plot_cfg.lfc_decimals)
                    else:
                        text = ""

                    color = "white" if abs(value) >= white_threshold else "black"
                    ax.text(x_centers[j], y_centers[i], text, ha="center", va="center", fontsize=8, color=color)

        if cfg.split_after_timepoint is not None and cfg.split_after_timepoint in plot_df.columns:
            split_pos = plot_df.columns.get_loc(cfg.split_after_timepoint) + 1
            ax.axvline(split_pos, color="black", linewidth=1.2)

        positive, negative = self.data._positive_negative_labels(effect_col)

        ax.set_title(
            f"ANCOM-BC2 log fold change ({positive} vs {negative})\n"
            f"{subset.title} | cells: {plot_cfg.cell_text_mode} | "
            f"taxa with ≥{cfg.min_sig_cells_per_taxon} significant cell(s)",
            fontsize=13,
            pad=12,
        )

        ax.set_xlabel("Timepoint")
        ax.set_ylabel("Taxon")

        ax.set_xticks(x_centers)
        ax.set_xticklabels(plot_df.columns, rotation=45, ha="right", fontsize=8)

        labels = _unique_labels([str(x) for x in plot_df.index])
        ax.set_yticks(y_centers)
        ax.set_yticklabels(labels, fontsize=7)

        for tick_label in ax.get_yticklabels():
            if _is_highlight(tick_label.get_text(), plot_cfg.highlight_taxa):
                tick_label.set_fontweight("bold")

        top_ax = ax.secondary_xaxis("top")
        top_ax.set_xticks(x_centers)

        top_labels = []
        for tp in plot_df.columns:
            counts = group_counts.get(tp, {})
            pos_n = counts.get(positive, "?")

            if negative == "reference":
                neg_n = sum(int(v) for k, v in counts.items() if k != positive)
            else:
                neg_n = counts.get(negative, "?")

            top_labels.append(f"{positive}={pos_n} | {negative}={neg_n}")

        top_ax.set_xticklabels(top_labels, rotation=45, ha="left", fontsize=8)
        top_ax.tick_params(axis="x", pad=6)

        cbar_ax = fig.add_axes([0.96, 0.60, 0.015, 0.18])
        cbar = fig.colorbar(mesh, cax=cbar_ax)
        cbar.set_label(
            f"ANCOM-BC2 log fold change\nred = higher in {positive}\nblue = higher in {negative}",
            fontsize=9,
        )
        cbar.ax.tick_params(labelsize=8)

        plt.subplots_adjust(left=0.52, right=0.95, top=0.82, bottom=0.10)

        out_base = (
            f"ancombc2_lfc_heatmap_"
            f"{cfg.variable_name}_"
            f"{subset.label}_"
            f"min{cfg.min_sig_cells_per_taxon}sig"
        )

        if save_png:
            out_png = os.path.join(cfg.output_dir, out_base + ".png")
            plt.savefig(out_png, dpi=300, bbox_inches="tight")
            print(f"Saved: {out_png}")

        if save_pdf:
            out_pdf = os.path.join(cfg.output_dir, out_base + ".pdf")
            plt.savefig(out_pdf, bbox_inches="tight")
            print(f"Saved: {out_pdf}")

        if show:
            plt.show()
        else:
            plt.close(fig)

        return fig, ax
