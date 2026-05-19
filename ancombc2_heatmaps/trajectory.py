from __future__ import annotations

from typing import Optional, Sequence

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter

from .config import Subset
from .data import ANCOMData


def clean_float_formatter(x, pos):
    return f"{x:.4f}".rstrip("0").rstrip(".")


def sig_to_label(is_sig: bool) -> str:
    return "*" if is_sig else "ns"


class TrajectoryPlotter:
    def __init__(self, data: ANCOMData):
        self.data = data
        self.config = data.config
        sns.set_theme(
            style=self.config.plot.sns_style,
            context=self.config.plot.sns_context,
        )

    def _apply_baseline(self, df: pd.DataFrame) -> pd.DataFrame:
        df = df.copy()
        df["tp_plot"] = df["tp"]

        if self.config.plot.merge_baselines:
            df["tp_plot"] = df["tp_plot"].replace({
                old: self.config.plot.merged_baseline_value
                for old in self.config.plot.baseline_values
            })

        return df

    def _comparison_levels(self, levels: Optional[Sequence[str]]) -> list[str]:
        if levels is not None:
            return list(levels)

        meta = self.data.metadata()
        return sorted(meta[self.config.group_col].dropna().unique())

    def _line_ylim(self, df: pd.DataFrame) -> tuple[float, float]:
        if df.empty:
            return (0, 1)

        upper = []
        for _, group_df in df.groupby(["tp_plot", "group"]):
            values = pd.to_numeric(group_df["abundance"], errors="coerce").dropna().to_numpy()
            if len(values) == 0:
                continue

            center = np.median(values) if self.config.plot.estimator == "median" else np.mean(values)

            if self.config.plot.error_style == "ci":
                rng = np.random.default_rng(42)
                func = np.median if self.config.plot.estimator == "median" else np.mean
                boots = [func(rng.choice(values, size=len(values), replace=True)) for _ in range(1000)]
                hi = np.percentile(boots, 97.5)
            else:
                hi = np.percentile(values, 75)

            upper.append(max(center, hi))

        if not upper:
            return (0, 1)

        ymax = max(upper)
        if not np.isfinite(ymax) or ymax <= 0:
            return (0, 1)

        return (0, ymax * 1.08)

    @staticmethod
    def _boxplot_whisker_bounds(values: pd.Series) -> tuple[float, float]:
        values = pd.to_numeric(values, errors="coerce").dropna()
        if values.empty:
            return np.nan, np.nan

        q1 = values.quantile(0.25)
        q3 = values.quantile(0.75)
        iqr = q3 - q1

        lower = q1 - 1.5 * iqr
        upper = q3 + 1.5 * iqr

        visible = values[(values >= lower) & (values <= upper)]
        if visible.empty:
            return float(values.min()), float(values.max())

        return float(visible.min()), float(visible.max())

    def _boxplot_ylim(self, df: pd.DataFrame, comparison_levels: Sequence[str]) -> tuple[float, float]:
        lows = []
        highs = []

        for tp_plot in sorted(df["tp_plot"].dropna().unique()):
            for group in comparison_levels:
                values = df.loc[
                    (df["tp_plot"] == tp_plot) & (df["group"] == group),
                    "abundance",
                ]

                low, high = self._boxplot_whisker_bounds(values)
                if not np.isnan(low):
                    lows.append(low)
                if not np.isnan(high):
                    highs.append(high)

        if lows and highs:
            y_min, y_max = min(lows), max(highs)
        else:
            y_min, y_max = float(df["abundance"].min()), float(df["abundance"].max())

        if y_min == y_max:
            pad = max(abs(y_max) * 0.1, 0.001)
            return max(0, y_min - pad), y_max + pad

        y_range = y_max - y_min
        return max(0, y_min - y_range * 0.08), y_max + y_range * 0.22

    def _add_significance_labels(self, ax, df: pd.DataFrame, sig_map: dict, x_is_categorical: bool):
        if not self.config.plot.show_significance or sig_map is None:
            return

        y0, y1 = ax.get_ylim()
        y_range = y1 - y0
        offset = y_range * 0.035

        tp_vals = sorted(df["tp_plot"].dropna().unique())

        for i, tp_plot in enumerate(tp_vals):
            g = df[df["tp_plot"] == tp_plot]
            if g.empty:
                continue

            original_tps = sorted(pd.to_numeric(g["tp"], errors="coerce").dropna().unique())

            if tp_plot == self.config.plot.merged_baseline_value and self.config.plot.merge_baselines:
                is_sig = any(sig_map.get(tp, False) for tp in self.config.plot.baseline_values)
            else:
                is_sig = any(sig_map.get(tp, False) for tp in original_tps)

            visible_upper = pd.to_numeric(g["abundance"], errors="coerce").quantile(0.75)
            if pd.isna(visible_upper):
                continue

            x = i if x_is_categorical else tp_plot
            y = min(visible_upper + offset, y1 - y_range * 0.04)

            ax.text(
                x,
                y,
                sig_to_label(bool(is_sig)),
                ha="center",
                va="bottom",
                fontsize=12,
                fontweight="bold",
                color="black",
            )

    def _clean_legend(self, ax):
        handles, labels = ax.get_legend_handles_labels()

        clean_handles = []
        clean_labels = []
        seen = set()

        for handle, label in zip(handles, labels):
            if label not in seen and not str(label).endswith("trend"):
                clean_handles.append(handle)
                clean_labels.append(label)
                seen.add(label)

        if self.config.plot.show_significance:
            star = Line2D(
                [0],
                [0],
                color="black",
                lw=0,
                label="* ANCOM significant\nns : not significant",
            )
            clean_handles.append(star)
            clean_labels.append(star.get_label())

        ax.legend(
            handles=clean_handles,
            labels=clean_labels,
            loc="upper left",
            fontsize=10,
            frameon=True,
            framealpha=0.9,
        )

    def line(
        self,
        taxon_query: str,
        subset: Subset,
        comparison_levels: Optional[Sequence[str]] = None,
    ):
        comparison_levels = self._comparison_levels(comparison_levels)

        df = self.data.taxon_timeseries(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

        if df.empty:
            print(f"[NO TRAJECTORY] No data: {taxon_query} — {subset.title}")
            return None

        df = self._apply_baseline(df)

        fig, ax = plt.subplots(figsize=self.config.plot.figsize)

        if self.config.plot.show_individual_lines:
            sns.lineplot(
                data=df,
                x="tp_plot",
                y="abundance",
                hue="group",
                style="group",
                units="subject",
                estimator=None,
                alpha=0.35,
                linewidth=1,
                legend=False,
                ax=ax,
            )

        errorbar = ("ci", 95) if self.config.plot.error_style == "ci" else ("pi", 50)

        sns.lineplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            style="group",
            estimator=self.config.plot.estimator,
            errorbar=errorbar,
            marker="o",
            linewidth=2.5,
            hue_order=comparison_levels,
            style_order=comparison_levels,
            ax=ax,
        )

        estimator_label = "median" if self.config.plot.estimator == "median" else "mean"
        error_label = "± 95% CI" if self.config.plot.error_style == "ci" else "± IQR"

        ax.set_title(f"{taxon_query} — {subset.title}\n({estimator_label} {error_label})", pad=15)
        ax.set_ylabel(self.config.plot.y_label)
        ax.set_xlabel("")
        ax.yaxis.set_major_formatter(FuncFormatter(clean_float_formatter))

        tp_vals = sorted(df["tp_plot"].dropna().unique())
        tp_labels = [self.config.timepoint_label_map.get(tp, str(tp)) for tp in tp_vals]

        ax.set_xticks(tp_vals)
        ax.set_xticklabels(tp_labels, rotation=35, ha="right")

        if self.config.plot.y_lim == "auto":
            ax.set_ylim(self._line_ylim(df))
        elif isinstance(self.config.plot.y_lim, tuple):
            ax.set_ylim(self.config.plot.y_lim)

        sig_map = self.data.significance_map(taxon_query, subset)
        self._add_significance_labels(ax, df, sig_map, x_is_categorical=False)
        self._clean_legend(ax)

        plt.tight_layout()
        plt.show()

        return fig, ax

    def boxplot(
        self,
        taxon_query: str,
        subset: Subset,
        comparison_levels: Optional[Sequence[str]] = None,
        show_trend: bool = True,
        trend_order: int = 2,
        jitter: float = 0.18,
    ):
        comparison_levels = self._comparison_levels(comparison_levels)

        df = self.data.taxon_timeseries(
            taxon_query=taxon_query,
            subset=subset,
            comparison_levels=comparison_levels,
        )

        if df.empty:
            print(f"[NO BOXPLOT] No data: {taxon_query} — {subset.title}")
            return None

        df = self._apply_baseline(df)
        df["tp_plot"] = pd.to_numeric(df["tp_plot"], errors="coerce")
        df["abundance"] = pd.to_numeric(df["abundance"], errors="coerce")
        df = df.dropna(subset=["tp_plot", "abundance", "group"])

        if df.empty:
            print(f"[NO BOXPLOT] No usable numeric data: {taxon_query} — {subset.title}")
            return None

        fig, ax = plt.subplots(figsize=self.config.plot.figsize)

        palette = dict(
            zip(
                comparison_levels,
                sns.color_palette(n_colors=len(comparison_levels)),
            )
        )

        sns.boxplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            hue_order=comparison_levels,
            palette=palette,
            dodge=True,
            showfliers=False,
            width=0.65,
            linewidth=1.2,
            ax=ax,
        )

        sns.stripplot(
            data=df,
            x="tp_plot",
            y="abundance",
            hue="group",
            hue_order=comparison_levels,
            palette=palette,
            dodge=True,
            jitter=jitter,
            alpha=0.65,
            size=3.5,
            linewidth=0.3,
            edgecolor="black",
            ax=ax,
        )

        tp_vals = sorted(df["tp_plot"].dropna().unique())
        x_positions = {tp: i for i, tp in enumerate(tp_vals)}

        if show_trend:
            grouped = df.groupby(["tp_plot", "group"], as_index=False)["abundance"].median()

            for group in comparison_levels:
                g = grouped[grouped["group"] == group].copy()
                if g.empty:
                    continue

                g["x_pos"] = g["tp_plot"].map(x_positions)
                g = g.dropna(subset=["x_pos", "abundance"])

                if len(g) < 2:
                    continue

                x = g["x_pos"].to_numpy(dtype=float)
                y = g["abundance"].to_numpy(dtype=float)

                order = min(trend_order, len(g) - 1)

                try:
                    coeffs = np.polyfit(x, y, deg=order)
                    poly = np.poly1d(coeffs)
                    x_smooth = np.linspace(x.min(), x.max(), 200)
                    y_smooth = poly(x_smooth)

                    ax.plot(
                        x_smooth,
                        y_smooth,
                        color=palette[group],
                        linewidth=2.2,
                        alpha=0.9,
                        label=f"{group} trend",
                    )
                except Exception:
                    ax.plot(
                        x,
                        y,
                        color=palette[group],
                        linewidth=2.2,
                        alpha=0.9,
                        label=f"{group} trend",
                    )

        ax.set_title(f"{taxon_query} — {subset.title}\n(boxplots + sample points + approximate trend)", pad=15)
        ax.set_ylabel(self.config.plot.y_label)
        ax.set_xlabel("")
        ax.yaxis.set_major_formatter(FuncFormatter(clean_float_formatter))

        tp_labels = [self.config.timepoint_label_map.get(tp, str(tp)) for tp in tp_vals]
        ax.set_xticks(range(len(tp_vals)))
        ax.set_xticklabels(tp_labels, rotation=35, ha="right")

        if self.config.plot.y_lim == "auto":
            ax.set_ylim(self._boxplot_ylim(df, comparison_levels))
        elif isinstance(self.config.plot.y_lim, tuple):
            ax.set_ylim(self.config.plot.y_lim)

        sig_map = self.data.significance_map(taxon_query, subset)
        self._add_significance_labels(ax, df, sig_map, x_is_categorical=True)
        self._clean_legend(ax)

        plt.tight_layout()
        plt.show()

        return fig, ax
