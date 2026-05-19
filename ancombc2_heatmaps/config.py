from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple, Union


@dataclass
class Subset:
    label: str
    title: Optional[str] = None
    filters: Dict[str, Union[str, Sequence[str]]] = field(default_factory=dict)

    def __post_init__(self):
        if self.title is None:
            self.title = self.label


@dataclass
class PlotConfig:
    figsize: Tuple[float, float] = (12, 8)
    y_label: str = "relative abundance"

    estimator: str = "mean"          # "mean" or "median"
    error_style: str = "iqr"         # "iqr" or "ci"
    show_individual_lines: bool = False
    show_significance: bool = True

    merge_baselines: bool = False
    baseline_values: Tuple[Union[int, float], ...] = (-7, -4, -1)
    merged_baseline_value: Union[int, float] = 0

    y_lim: Union[str, Tuple[float, float]] = "auto"
    sns_style: str = "whitegrid"
    sns_context: str = "talk"

    heatmap_cmap: str = "RdBu_r"
    cell_text_mode: str = "relative_abundance"  # "relative_abundance", "lfc", "none"
    lfc_decimals: int = 2
    highlight_taxa: List[str] = field(default_factory=list)

    save_png: bool = True
    save_pdf: bool = True
    show: bool = True


@dataclass
class ANCOMConfig:
    # Files
    metadata_path: str
    output_dir: str

    # Metadata columns
    sample_col: str
    timepoint_col: str
    group_col: str
    subject_col: Optional[str] = None

    # Timepoints
    timepoints: List[str] = field(default_factory=list)
    timepoint_map: Dict[str, str] = field(default_factory=dict)
    timepoint_numeric_map: Dict[str, Union[int, float]] = field(default_factory=dict)
    timepoint_label_map: Dict[Union[int, float], str] = field(default_factory=dict)

    # Optional metadata filtering
    allowed_values: Dict[str, List[str]] = field(default_factory=dict)

    # Path mode 1: template paths
    table_base: str = ""
    ancom_base: str = ""
    table_template: str = "{timepoint}/table_{timepoint}_{subset_label}.qza"
    ancom_template: str = "{timepoint}/table_{timepoint}_{subset_label}_ANCOMB_exported"

    # Path mode 2: explicit fixed paths
    table_paths: Optional[Dict[str, str]] = None
    ancom_paths: Optional[Dict[str, str]] = None

    # Export filenames
    lfc_filename: str = "lfc.jsonl"
    q_filename: str = "q.jsonl"
    diff_filename: str = "diff.jsonl"

    # ANCOM effect handling
    variable_name: Optional[str] = None
    effect_column: Optional[str] = None
    positive_class: Optional[str] = None
    negative_class: Optional[str] = None
    invert_sign: bool = False
    q_cutoff: float = 0.05

    # Default taxa for workflow.trajectories() / workflow.boxplots()
    taxa: List[str] = field(default_factory=list)

    # Heatmap filtering
    min_sig_cells_per_taxon: int = 1
    split_after_timepoint: Optional[str] = None
    remove_empty_rows_after_masking: bool = True

    # Plot settings
    plot: PlotConfig = field(default_factory=PlotConfig)

    def __post_init__(self):
        if self.variable_name is None:
            self.variable_name = self.group_col

        if not self.timepoints:
            raise ValueError("ANCOMConfig.timepoints must contain at least one timepoint.")

        if not self.timepoint_numeric_map:
            self.timepoint_numeric_map = {tp: i for i, tp in enumerate(self.timepoints)}

        if not self.timepoint_label_map:
            self.timepoint_label_map = {
                self.timepoint_numeric_map[tp]: tp
                for tp in self.timepoints
                if tp in self.timepoint_numeric_map
            }
