from __future__ import annotations

import os
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd
import qiime2
from biom import Table

from .config import ANCOMConfig, Subset
from .taxonomy import list_available_queries, match_taxa, normalize_taxon_label


def read_table_auto(path: str) -> pd.DataFrame:
    if str(path).endswith(".csv"):
        return pd.read_csv(path)
    return pd.read_csv(path, sep="\t")


def load_qza_table_as_df(path: str) -> pd.DataFrame:
    artifact = qiime2.Artifact.load(path)
    table = artifact.view(Table)
    df = table.to_dataframe(dense=True)
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)
    return df


def read_ancom_jsonl(path: str) -> pd.DataFrame:
    """
    Read a QIIME2 ANCOM-BC2 JSONL export table.

    QIIME2 JSONL exports can contain one metadata/header row before the real data.
    This function removes rows without a valid taxon entry.
    """
    df = pd.read_json(path, lines=True)

    if "taxon" in df.columns:
        df = df[df["taxon"].notna()].copy()

    return df


class ANCOMData:
    def __init__(self, config: ANCOMConfig):
        self.config = config
        Path(config.output_dir).mkdir(parents=True, exist_ok=True)
        self._metadata_cache: Optional[pd.DataFrame] = None
        self._rel_cache: Dict[Tuple[str, str], pd.DataFrame] = {}

    # ------------------------------------------------------------------
    # Metadata
    # ------------------------------------------------------------------
    def metadata(self) -> pd.DataFrame:
        if self._metadata_cache is None:
            self._metadata_cache = self._load_metadata()
        return self._metadata_cache.copy()

    def _load_metadata(self) -> pd.DataFrame:
        cfg = self.config

        if not os.path.exists(cfg.metadata_path):
            raise FileNotFoundError(f"Metadata file not found: {cfg.metadata_path}")

        meta = read_table_auto(cfg.metadata_path)

        required = [cfg.sample_col, cfg.timepoint_col, cfg.group_col]
        if cfg.subject_col is not None:
            required.append(cfg.subject_col)

        missing = [col for col in required if col not in meta.columns]
        if missing:
            raise ValueError(f"Missing metadata columns: {missing}")

        meta = meta.copy()

        for col in meta.columns:
            if meta[col].dtype == object:
                meta[col] = meta[col].astype(str).str.strip()

        meta[cfg.sample_col] = meta[cfg.sample_col].astype(str)
        meta[cfg.timepoint_col] = meta[cfg.timepoint_col].astype(str)

        if cfg.timepoint_map:
            meta[cfg.timepoint_col] = meta[cfg.timepoint_col].map(
                lambda x: cfg.timepoint_map.get(x, x)
            )

        if cfg.timepoints:
            meta = meta[meta[cfg.timepoint_col].isin(cfg.timepoints)].copy()

        for col, allowed in cfg.allowed_values.items():
            if col not in meta.columns:
                raise ValueError(f"allowed_values references unknown column: {col}")
            meta = meta[meta[col].isin(allowed)].copy()

        meta["_tp_num"] = meta[cfg.timepoint_col].map(cfg.timepoint_numeric_map)

        return meta

    def filter_metadata(self, subset: Optional[Subset]) -> pd.DataFrame:
        meta = self.metadata()

        if subset is None:
            return meta

        out = meta.copy()
        for col, value in subset.filters.items():
            if col not in out.columns:
                raise ValueError(f"Subset filter references unknown column: {col}")

            if isinstance(value, (list, tuple, set)):
                out = out[out[col].isin(list(value))].copy()
            else:
                out = out[out[col] == value].copy()

        return out

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    def table_path(self, timepoint: str, subset: Subset) -> str:
        cfg = self.config

        if cfg.table_paths is not None:
            if timepoint not in cfg.table_paths:
                raise ValueError(f"No table path configured for timepoint: {timepoint}")
            return cfg.table_paths[timepoint]

        rel = cfg.table_template.format(
            timepoint=timepoint,
            subset_label=subset.label,
            variable_name=cfg.variable_name,
        )
        return os.path.join(cfg.table_base, rel)

    def ancom_path(self, timepoint: str, subset: Subset) -> str:
        cfg = self.config

        if cfg.ancom_paths is not None:
            if timepoint not in cfg.ancom_paths:
                raise ValueError(f"No ANCOM path configured for timepoint: {timepoint}")
            return cfg.ancom_paths[timepoint]

        rel = cfg.ancom_template.format(
            timepoint=timepoint,
            subset_label=subset.label,
            variable_name=cfg.variable_name,
        )
        return os.path.join(cfg.ancom_base, rel)

    def explain_paths(self, subset: Subset) -> pd.DataFrame:
        rows = []
        for tp in self.config.timepoints:
            table = self.table_path(tp, subset)
            ancom = self.ancom_path(tp, subset)
            rows.append({
                "timepoint": tp,
                "table_path": table,
                "table_exists": os.path.exists(table),
                "ancom_path": ancom,
                "ancom_exists": os.path.isdir(ancom),
            })
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # QZA tables / relative abundance
    # ------------------------------------------------------------------
    def relative_abundance(self, timepoint: str, subset: Subset) -> pd.DataFrame:
        cache_key = (timepoint, subset.label)

        if cache_key in self._rel_cache:
            return self._rel_cache[cache_key].copy()

        path = self.table_path(timepoint, subset)
        if not os.path.exists(path):
            raise FileNotFoundError(f"QZA table not found: {path}")

        table = load_qza_table_as_df(path)
        table.index = [normalize_taxon_label(x) for x in table.index]
        table = table.groupby(table.index).sum()

        col_sums = table.sum(axis=0)
        col_sums[col_sums == 0] = np.nan

        rel = table.div(col_sums, axis=1)
        self._rel_cache[cache_key] = rel

        return rel.copy()

    def mean_relative_abundance(
        self,
        subset: Subset,
    ) -> Tuple[pd.DataFrame, Dict[str, Dict[str, int]]]:
        cfg = self.config
        meta = self.filter_metadata(subset)

        mean_by_tp = {}
        group_counts = {}

        for tp in cfg.timepoints:
            try:
                rel = self.relative_abundance(tp, subset)
            except FileNotFoundError:
                continue

            tp_meta = meta[meta[cfg.timepoint_col] == tp].copy()
            samples = [s for s in tp_meta[cfg.sample_col] if s in rel.columns]

            if not samples:
                continue

            mean_by_tp[tp] = rel[samples].mean(axis=1, skipna=True)

            used_meta = tp_meta[tp_meta[cfg.sample_col].isin(samples)]
            group_counts[tp] = (
                used_meta[cfg.group_col]
                .value_counts(dropna=False)
                .astype(int)
                .to_dict()
            )

        return pd.DataFrame(mean_by_tp), group_counts

    # ------------------------------------------------------------------
    # ANCOM exports
    # ------------------------------------------------------------------
    def detect_effect_column(self, df: pd.DataFrame) -> str:
        cfg = self.config

        if cfg.effect_column is not None:
            if cfg.effect_column not in df.columns:
                raise ValueError(
                    f"Configured effect_column '{cfg.effect_column}' not found. "
                    f"Available columns: {list(df.columns)}"
                )
            return cfg.effect_column

        prefix = f"{cfg.variable_name}::"
        candidates = [col for col in df.columns if col.startswith(prefix)]

        if len(candidates) == 1:
            return candidates[0]

        if len(candidates) > 1:
            raise ValueError(
                f"Multiple effect columns found: {candidates}. "
                f"Set effect_column explicitly in ANCOMConfig."
            )

        raise ValueError(
            f"No effect column found with prefix '{prefix}'. "
            f"Available columns: {list(df.columns)}"
        )

    def _positive_negative_labels(self, effect_col: Optional[str]) -> Tuple[str, str]:
        cfg = self.config

        if cfg.positive_class is not None:
            positive = cfg.positive_class
        elif effect_col and "::" in effect_col:
            positive = effect_col.split("::", 1)[1]
        else:
            positive = cfg.variable_name or cfg.group_col

        if cfg.negative_class is not None:
            negative = cfg.negative_class
        else:
            allowed = cfg.allowed_values.get(cfg.group_col)
            if allowed is not None:
                others = [x for x in allowed if x != positive]
                negative = others[0] if len(others) == 1 else "reference"
            else:
                negative = "reference"

        return positive, negative

    def _significance_filename(self) -> str:
        cfg = self.config

        if cfg.significance_metric == "q":
            return cfg.q_filename

        if cfg.significance_metric == "p":
            return cfg.p_filename

        raise ValueError(
            "significance_metric must be either 'q' or 'p'. "
            f"Got: {cfg.significance_metric}"
        )

    def _require_diff_for_significance(self) -> bool:
        cfg = self.config

        if cfg.require_diff_for_significance is not None:
            return bool(cfg.require_diff_for_significance)

        # Default behavior:
        # - q-values: require ANCOM-BC2 diff == True
        # - p-values: only use unadjusted p < cutoff
        if cfg.significance_metric == "q":
            return True

        if cfg.significance_metric == "p":
            return False

        raise ValueError(
            "significance_metric must be either 'q' or 'p'. "
            f"Got: {cfg.significance_metric}"
        )

    def read_ancom_timepoint(
        self,
        timepoint: str,
        subset: Subset,
    ) -> Optional[pd.DataFrame]:
        cfg = self.config
        export_dir = self.ancom_path(timepoint, subset)

        if not os.path.isdir(export_dir):
            return None

        lfc_path = os.path.join(export_dir, cfg.lfc_filename)
        sig_path = os.path.join(export_dir, self._significance_filename())
        diff_path = os.path.join(export_dir, cfg.diff_filename)

        if not (
            os.path.exists(lfc_path)
            and os.path.exists(sig_path)
            and os.path.exists(diff_path)
        ):
            return None

        lfc = read_ancom_jsonl(lfc_path)
        sig = read_ancom_jsonl(sig_path)
        diff = read_ancom_jsonl(diff_path)

        tax_col_lfc = "taxon" if "taxon" in lfc.columns else lfc.columns[0]
        tax_col_sig = "taxon" if "taxon" in sig.columns else sig.columns[0]
        tax_col_diff = "taxon" if "taxon" in diff.columns else diff.columns[0]

        effect_col = self.detect_effect_column(lfc)

        if effect_col not in sig.columns or effect_col not in diff.columns:
            raise ValueError(
                f"Effect column '{effect_col}' missing in "
                f"{cfg.significance_metric}/diff files for {export_dir}"
            )

        sig_value_name = cfg.significance_metric

        out = (
            lfc[[tax_col_lfc, effect_col]]
            .rename(columns={tax_col_lfc: "taxon_raw", effect_col: "lfc_raw"})
            .merge(
                sig[[tax_col_sig, effect_col]].rename(
                    columns={
                        tax_col_sig: "taxon_raw",
                        effect_col: sig_value_name,
                    }
                ),
                on="taxon_raw",
                how="inner",
            )
            .merge(
                diff[[tax_col_diff, effect_col]].rename(
                    columns={
                        tax_col_diff: "taxon_raw",
                        effect_col: "diff",
                    }
                ),
                on="taxon_raw",
                how="inner",
            )
        )

        out["taxon"] = out["taxon_raw"].apply(normalize_taxon_label)
        out["lfc"] = pd.to_numeric(out["lfc_raw"], errors="coerce")

        if cfg.invert_sign:
            out["lfc"] = -out["lfc"]

        out[sig_value_name] = pd.to_numeric(out[sig_value_name], errors="coerce")
        out["significance_value"] = out[sig_value_name]

        # Keep a q column for backward compatibility with older code.
        # If significance_metric == "p", this q column contains the selected
        # significance value used by the package, not an FDR-adjusted q-value.
        out["q"] = out["significance_value"]

        out["diff"] = out["diff"].astype(str).isin(["True", "true", "1"])

        if self._require_diff_for_significance():
            out["significant"] = (
                (out["significance_value"] < cfg.q_cutoff)
                & out["diff"]
            )
        else:
            out["significant"] = out["significance_value"] < cfg.q_cutoff

        out["timepoint"] = timepoint
        out["effect_column"] = effect_col
        out["significance_metric"] = cfg.significance_metric

        return out

    def ancom_matrices(
        self,
        subset: Subset,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, Optional[str]]:
        lfc_by_tp = {}
        sig_by_tp = {}
        effect_col = None

        for tp in self.config.timepoints:
            result = self.read_ancom_timepoint(tp, subset)
            if result is None or result.empty:
                continue

            effect_col = str(result["effect_column"].iloc[0])

            lfc_by_tp[tp] = result.groupby("taxon")["lfc"].mean()
            sig_by_tp[tp] = result.groupby("taxon")["significant"].any()

        return pd.DataFrame(lfc_by_tp), pd.DataFrame(sig_by_tp), effect_col

    def significance_map(self, taxon_query: str, subset: Subset) -> Dict[float, bool]:
        sig_map = {}

        for tp in self.config.timepoints:
            result = self.read_ancom_timepoint(tp, subset)
            tp_num = self.config.timepoint_numeric_map.get(tp)

            if result is None or result.empty or tp_num is None:
                continue

            sig_by_taxon = result.groupby("taxon")["significant"].any()
            matches = match_taxa(sig_by_taxon.index, taxon_query)

            sig_map[tp_num] = bool(sig_by_taxon.loc[matches].any()) if matches else False

        return sig_map

    # ------------------------------------------------------------------
    # Trajectory data
    # ------------------------------------------------------------------
    def taxon_timeseries(
        self,
        taxon_query: str,
        subset: Subset,
        comparison_levels: Optional[Sequence[str]] = None,
    ) -> pd.DataFrame:
        cfg = self.config
        meta = self.filter_metadata(subset)

        if comparison_levels is not None:
            meta = meta[meta[cfg.group_col].isin(comparison_levels)].copy()

        rows = []

        for tp in cfg.timepoints:
            try:
                rel = self.relative_abundance(tp, subset)
            except FileNotFoundError:
                continue

            matches = match_taxa(rel.index, taxon_query)
            if not matches:
                continue

            abundance = rel.loc[matches].sum(axis=0)

            tp_meta = meta[meta[cfg.timepoint_col] == tp].copy()
            tp_meta = tp_meta[tp_meta[cfg.sample_col].isin(abundance.index)]

            for _, row in tp_meta.iterrows():
                rows.append({
                    "sample": row[cfg.sample_col],
                    "subject": row[cfg.subject_col] if cfg.subject_col else row[cfg.sample_col],
                    "group": row[cfg.group_col],
                    "timepoint": tp,
                    "tp": cfg.timepoint_numeric_map.get(tp),
                    "abundance": abundance.get(row[cfg.sample_col], np.nan),
                    "taxon_query": taxon_query,
                })

        df = pd.DataFrame(rows)
        if df.empty:
            return pd.DataFrame(
                columns=[
                    "sample",
                    "subject",
                    "group",
                    "timepoint",
                    "tp",
                    "abundance",
                    "taxon_query",
                ]
            )

        return df.dropna(subset=["tp", "abundance", "group"])

    def available_taxa(
        self,
        timepoint: Optional[str] = None,
        subset: Optional[Subset] = None,
    ) -> Dict[str, List[str]]:
        if subset is None:
            subset = Subset(label="all", title="all", filters={})

        if timepoint is None:
            timepoint = self.config.timepoints[0]

        rel = self.relative_abundance(timepoint, subset)
        return list_available_queries(list(rel.index))