[!FIXME]
# Heatmap packgage 

## Table of Contents

- [Summary](#summary)
- [Installation](#installation)
- [Input Data](#input-data)
- [Heatmap Workflow](#heatmap-workflow)
- [Trajectory Workflow](#trajectory-workflow)
- [Boxplot Workflow](#boxplot-workflow)
- [ANCOM-BC2](#ancom-bc2) # if required


---

## Summary

blablabla

---

## Installation

The installation of the package can be done from GitHub via git in the console:

[!FIXME] # maybe a change after making it public
```
pip install git+ssh://git@github.com/jlab/ancombc2-heatmaps.git  
```



Please be aware that you have to activate the correct conda environment with an Qiime2 Version that is not too old.<br>
I recommend "qiime2-amplicon-2026.1".

---

To see your available environments use:

```
conda env list
```

and to activate the correct one use:

```
conda activate qiime2-amplicon-2026.1
```

if you don't have a working Qiime2 Environment you can find a good YAML here:

https://github.com/qiime2/distributions

```
wget https://raw.githubusercontent.com/qiime2/distributions/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml \
  -O qiime2-amplicon-2026.1.yml
```

and install it with

```
conda env create \
  -n qiime2-amplicon-2026.1 \
  -f ~/qiime2-amplicon-2026.1.yml
```


---

After installation you can check the version with

```
qiime --version
```

and import it for example in your notebook with

```
import ancombc2_heatmaps as ah
```


You also have to import numpy and Path from pathlib.

```
from pathlib import Path
import pandas as pd
```
---

## Input Data

## Required input files

After importing the package, the plotting functions require three main types of input files:

1. A metadata table
2. QIIME2 feature tables
3. Exported ANCOM-BC2 result files

---

## 1. Metadata table

The metadata file must be a tab-separated text file (`.tsv` or `.txt`).

Example:

```text
#SampleID	time_point	description_of_treatment	mice_model	sex
sample_1	baseline1	sham	WT	female
sample_2	baseline1	irradiated	WT	male
sample_3	day1	sham	Apc	female
```

The required columns depend on the configuration.

For heatmaps, the metadata must contain at least:

```text
sample_col
timepoint_col
comparison_col
```

Example:

```python
metadata=ah.MetadataConfig(
    sample_col="#SampleID",
    timepoint_col="time_point",
    comparison_col="description_of_treatment",
)
```

This means the metadata table must contain the columns:

```text
#SampleID
time_point
description_of_treatment
```

For trajectory plots, the metadata must contain at least:

```text
sample_col
timepoint_col
mouse_col
comparison_col
```

Example:

```python
metadata=ah.TrajectoryMetadataConfig(
    sample_col="#SampleID",
    timepoint_col="time_point",
    mouse_col="mouse_id",
    comparison_col="description_of_treatment",
)
```

---

## 2. QIIME2 feature tables

The package expects QIIME2 feature tables in `.qza` format.

These tables are used to calculate relative abundances for the heatmap cell annotations and trajectory plots.

Example file structure:

```text
collapsed_tables/
├── genus/
│   └── by_treatment/
│       ├── baseline1/
│       │   └── table_baseline1_all_genus.qza
│       ├── day1/
│       │   └── table_day1_all_genus.qza
│       └── day3/
│           └── table_day3_all_genus.qza
```

The exact filenames are defined in the configuration:

```python
paths=ah.PathConfig(
    base_table_dir="path/to/collapsed_tables/genus/by_treatment",
    table_template="{timepoint}/table_{timepoint}_{subset_label}.qza",
)
```

The `.qza` files must contain feature tables that can be loaded as BIOM tables through QIIME2.

---

## 3. Exported ANCOM-BC2 result files

The package expects exported ANCOM-BC2 result folders.

Each exported folder must contain:

```text
lfc.jsonl
q.jsonl
diff.jsonl
```

Example:

```text
ancombc2_results/
├── genus/
│   └── by_treatment/
│       ├── baseline1/
│       │   └── table_baseline1_all_genus_ANCOMBC2_exported/
│       │       ├── lfc.jsonl
│       │       ├── q.jsonl
│       │       └── diff.jsonl
│       ├── day1/
│       │   └── table_day1_all_genus_ANCOMBC2_exported/
│       │       ├── lfc.jsonl
│       │       ├── q.jsonl
│       │       └── diff.jsonl
```

The files must contain the ANCOM-BC2 output tables in JSON lines format.

More information about the required ANCOM-BC2 result files is provided in the [ANCOM-BC2](#ancom-bc2) section.

---

## Summary of required files

| File type | Format | Used for |
|---|---|---|
| Metadata table | `.tsv` or `.txt` | Sample information, groups, timepoints |
| Feature table | `.qza` | Relative abundance calculation |
| ANCOM-BC2 LFC table | `lfc.jsonl` | Log fold changes |
| ANCOM-BC2 q-value table | `q.jsonl` | Adjusted p-values |
| ANCOM-BC2 diff table | `diff.jsonl` | Significance status |

---

## Minimal required directory structure

```text
project/
├── metadata.tsv
├── collapsed_tables/
│   └── genus/
│       └── by_treatment/
│           └── baseline1/
│               └── table_baseline1_all_genus.qza
└── ancombc2_results/
    └── genus/
        └── by_treatment/
            └── baseline1/
                └── table_baseline1_all_genus_ANCOMBC2_exported/
                    ├── lfc.jsonl
                    ├── q.jsonl
                    └── diff.jsonl
```


---


## Heatmap Workflow


This section explains how to create ANCOM-BC2 heatmaps with the package.

The heatmap workflow is used to visualize exported ANCOM-BC2 results across multiple timepoints.  
It combines:

- ANCOM-BC2 log fold changes
- ANCOM-BC2 significance information
- mean relative abundance from QIIME2 feature tables
- metadata-based subgroup filtering

---

## Basic idea

The heatmap needs three types of input files:

| Input | Description |
|---|---|
| QIIME2 `.qza` feature tables | Used to calculate mean relative abundance |
| exported ANCOM-BC2 result folders | Contain `lfc.jsonl`, `q.jsonl`, and `diff.jsonl` |
| metadata file | Used for sample grouping, timepoints, and subset filtering |

The final output is a heatmap showing significant ANCOM-BC2 log fold changes across timepoints.

---

## 1. Import the package

```python
import ancombc2_heatmaps as ah
```

---

## 2. Define the timepoints

First, define the timepoints that should appear in the heatmap.

```python
timepoints = [
    "baseline1",
    "baseline2",
    "baseline3",
    "day1",
    "day3",
    "day7",
    "day14",
]
```

The order in this list is also the order used in the heatmap.

---

## 3. Create the heatmap configuration

The main configuration object is `HeatmapConfig`.

It contains several smaller configuration blocks:

| Config class | Purpose |
|---|---|
| `MetadataConfig` | Defines metadata columns and timepoint handling |
| `ComparisonConfig` | Defines the ANCOM-BC2 comparison variable |
| `PathConfig` | Defines input and output paths |
| `PlotTextConfig` | Customizes titles and labels |
| `HeatmapStyleConfig` | Customizes figure size, fonts, colors, and layout |
| `TaxonomyConfig` | Customizes taxon label formatting |

A minimal example:

```python
heatmap_config = ah.HeatmapConfig(
    metadata=ah.MetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        comparison_col="description_of_treatment",
        timepoints=timepoints,
        timepoint_map={
            "baseline_1": "baseline1",
            "baseline_2": "baseline2",
            "baseline_3": "baseline3",
            "day_1_post": "day1",
            "day_3_post": "day3",
            "day_7_post": "day7",
            "day_14_post": "day14",
        },
        allowed_values={
            "description_of_treatment": ["sham", "irradiated"],
        },
    ),

    comparison=ah.ComparisonConfig(
        variable_name="description_of_treatment",
    ),

    paths=ah.PathConfig(
        base_table_dir="/path/to/qza_tables",
        base_ancom_dir="/path/to/exported_ancom_results",
        metadata_path="/path/to/metadata.tsv",
        output_dir="/path/to/output",

        table_template="{timepoint}/table_{timepoint}_genus_ANCOM.qza",
        ancom_template="{timepoint}_treat_ANCOMB_exported",
    ),

    split_after_timepoint="baseline3",
)
```

---

## 4. MetadataConfig options

`MetadataConfig` tells the package how to read and interpret the metadata.

```python
metadata=ah.MetadataConfig(
    sample_col="sample_name",
    timepoint_col="time_point",
    comparison_col="description_of_treatment",
    timepoints=timepoints,
    timepoint_map={...},
    allowed_values={...},
)
```

### Options

| Option | Description |
|---|---|
| `sample_col` | Column containing sample IDs |
| `timepoint_col` | Column containing timepoint information |
| `comparison_col` | Column used for the main comparison |
| `timepoints` | Timepoints to include and their plotting order |
| `timepoint_map` | Optional renaming of metadata timepoints |
| `allowed_values` | Optional filtering of metadata values |

Example:

```python
allowed_values={
    "description_of_treatment": ["sham", "irradiated"],
}
```

This means only samples with treatment `sham` or `irradiated` are used.

This is usually only necessary if the metadata contains additional groups that should be excluded.

---

## 5. ComparisonConfig options

`ComparisonConfig` controls which ANCOM-BC2 effect column is used.

```python
comparison=ah.ComparisonConfig(
    variable_name="description_of_treatment",
)
```

The package searches for ANCOM-BC2 columns starting with:

```text
description_of_treatment::
```

For example:

```text
description_of_treatment::sham
```

### Options

| Option | Description |
|---|---|
| `variable_name` | Metadata variable used in ANCOM-BC2 |
| `positive_class` | Optional custom label for the positive group |
| `negative_class` | Optional custom label for the reference group |
| `effect_column` | Manually specify the exact ANCOM-BC2 effect column |
| `invert_sign` | Manually invert the log fold change sign |

Example with manual effect column:

(In general you don't need to modify this manually.)

```python
comparison=ah.ComparisonConfig(
    variable_name="description_of_treatment",
    effect_column="description_of_treatment::sham",
)
```

Example with custom labels:

```python
comparison=ah.ComparisonConfig(
    variable_name="description_of_treatment",
    positive_class="sham",
    negative_class="irradiated",
)
```

Example with sign inversion:

```python
comparison=ah.ComparisonConfig(
    variable_name="description_of_treatment",
    invert_sign=True,
)
```

Use `invert_sign=True` only if you explicitly want to reverse the ANCOM-BC2 log fold change direction.

---

## 6. PathConfig options

`PathConfig` defines where the package finds the input files and where it saves the output.

```python
paths=ah.PathConfig(
    base_table_dir="/path/to/qza_tables",
    base_ancom_dir="/path/to/exported_ancom_results",
    metadata_path="/path/to/metadata.tsv",
    output_dir="/path/to/output",

    table_template="{timepoint}/table_{timepoint}_{subset_label}.qza",
    ancom_template="{timepoint}/table_{timepoint}_{subset_label}_ANCOMB_exported",
)
```

### Required paths

| Option | Description |
|---|---|
| `base_table_dir` | Base folder containing QIIME2 `.qza` feature tables |
| `base_ancom_dir` | Base folder containing exported ANCOM-BC2 result folders |
| `metadata_path` | Path to metadata file |
| `output_dir` | Folder where heatmaps are saved |

### Template options

| Option | Description |
|---|---|
| `table_template` | Relative path pattern for `.qza` tables |
| `ancom_template` | Relative path pattern for ANCOM-BC2 export folders |
| `lfc_filename` | Filename for log fold changes, default: `lfc.jsonl` |
| `q_filename` | Filename for q-values, default: `q.jsonl` |
| `diff_filename` | Filename for significance flags, default: `diff.jsonl` |

The templates can use:

```text
{timepoint}
{subset_label}
```

Example:

```python
table_template="{timepoint}/table_{timepoint}_{subset_label}.qza"
```

For:

```python
timepoint = "day1"
subset.label = "WT"
```

this becomes:

```text
day1/table_day1_WT.qza
```

---

## 7. General HeatmapConfig options

These options are set directly inside `HeatmapConfig`.

```python
heatmap_config = ah.HeatmapConfig(
    ...,
    q_cutoff=0.05,
    min_sig_cells_per_taxon=1,
    split_after_timepoint="baseline3",
    remove_empty_rows_after_masking=True,
    cell_text_mode="relative_abundance",
    lfc_decimals=2,
)
```

### Options

| Option | Description |
|---|---|
| `q_cutoff` | q-value threshold for significance |
| `min_sig_cells_per_taxon` | Minimum number of significant cells required to keep a taxon |
| `split_after_timepoint` | Draws a vertical line after a selected timepoint |
| `remove_empty_rows_after_masking` | Removes taxa with no significant values after masking |
| `cell_text_mode` | Controls text inside heatmap cells |
| `lfc_decimals` | Number of decimals for LFC labels |

### Cell text modes

```python
cell_text_mode="relative_abundance"
```

Options:

| Value | Meaning |
|---|---|
| `"relative_abundance"` | Show mean relative abundance inside significant cells |
| `"lfc"` | Show log fold change values inside significant cells |
| `"none"` | Show no text inside cells |

Example:

```python
heatmap_config = ah.HeatmapConfig(
    ...,
    cell_text_mode="lfc",
    lfc_decimals=3,
)
```

---

## 8. Styling options

You can customize the appearance with `HeatmapStyleConfig`.

```python
style=ah.HeatmapStyleConfig(
    title_fontsize=13,
    axis_label_fontsize=11,
    xtick_fontsize=8,
    ytick_fontsize=7,
    celltext_fontsize=8,

    base_cell_height=0.25,
    base_cell_width=0.78,
    height_scale=6.0,

    heatmap_cmap="RdBu_r",
    missing_color="white",

    highlight_taxa=[
        "Akkermansia",
        "Muribaculaceae",
    ],
)
```

### Common style options

| Option | Description |
|---|---|
| `title_fontsize` | Title font size |
| `axis_label_fontsize` | Axis label font size |
| `xtick_fontsize` | X-axis label font size |
| `ytick_fontsize` | Y-axis taxon label font size |
| `celltext_fontsize` | Text size inside cells |
| `base_cell_height` | Base height of each heatmap row |
| `base_cell_width` | Base width of each heatmap column |
| `height_scale` | Strength of row-height scaling |
| `rowheight_ra_min` | Minimum relative abundance for row-height scaling |
| `rowheight_ra_max` | Maximum relative abundance for row-height scaling |
| `figure_min_height` | Minimum figure height |
| `figure_max_height` | Maximum figure height |
| `figure_min_width` | Minimum figure width |
| `figure_max_width` | Maximum figure width |
| `heatmap_cmap` | Matplotlib colormap |
| `missing_color` | Color for non-significant or missing cells |
| `edge_color` | Cell border color |
| `split_line_color` | Color of the vertical split line |
| `highlight_taxa` | Taxa names that should be highlighted in bold |

---

## 9. Custom title and labels

You can customize text with `PlotTextConfig`.

```python
text=ah.PlotTextConfig(
    title_template=(
        "ANCOM-BC2 log fold change ({positive_label} vs {negative_label})\n"
        "{subset_title}"
    ),
    x_label="Timepoint",
    y_label="Taxon",
    colorbar_template=(
        "ANCOM-BC2 log fold change\n"
        "red = higher in {positive_label}\n"
        "blue = higher in {negative_label}"
    ),
    top_axis_count_template="{positive_label}={positive_n} | {negative_label}={negative_n}",
)
```

Available placeholders include:

| Placeholder | Meaning |
|---|---|
| `{positive_label}` | Positive comparison group |
| `{negative_label}` | Reference group |
| `{subset_title}` | Title from the selected subset |
| `{min_sig_cells}` | Minimum significant cells per taxon |
| `{cell_text_description}` | Description of cell text mode |
| `{positive_n}` | Number of samples in positive group |
| `{negative_n}` | Number of samples in negative group |

---

## 10. Create the plotter

After defining the config, create the plotter:

```python
plotter = ah.ANCOMBC2HeatmapPlotter(
    heatmap_config
)
```

The plotter will use the config to find the required files and generate heatmaps.

---

## 11. Define a subset

A subset defines which samples are used.
Please be aware that you only use subsets you have a individually calculatet ancom export for.

```python
subset = ah.SubsetSpec(
    label="WT",
    title="WT mice",
    filters={
        "mice_model": "WT",
    },
)
```

### Subset options

| Option | Description |
|---|---|
| `label` | Used in file paths and output filenames |
| `title` | Used in the heatmap title |
| `filters` | Metadata filters applied before plotting |

Example with multiple filters:

```python
subset = ah.SubsetSpec(
    label="WT_female",
    title="WT female mice",
    filters={
        "mice_model": "WT",
        "sex": "female",
    },
)
```

Example without filtering:

```python
subset = ah.SubsetSpec(
    label="",
    title="All samples",
    filters={},
)
```

---

## 12. Generate one heatmap

```python
plotter.plot_subset(
    meta_df=plotter.load_metadata(),
    subset=subset,
)
```

This will:

1. Load and filter metadata
2. Load QIIME2 feature tables
3. Calculate mean relative abundance
4. Load exported ANCOM-BC2 results
5. Mask non-significant cells
6. Generate the heatmap
7. Save the output files

---

## 13. Control saving and display

```python
plotter.plot_subset(
    meta_df=plotter.load_metadata(),
    subset=subset,
    save_png=True,
    save_pdf=True,
    show=True,
)
```

### Options

| Option | Description |
|---|---|
| `save_png` | Save heatmap as PNG |
| `save_pdf` | Save heatmap as PDF |
| `show` | Display the plot directly |

For scripts or batch jobs, use:

```python
plotter.plot_subset(
    meta_df=plotter.load_metadata(),
    subset=subset,
    save_png=True,
    save_pdf=True,
    show=False,
)
```

---

## 14. Generate multiple heatmaps

To generate multiple heatmaps in one run:

```python
subsets = [
    ah.SubsetSpec(
        label="WT",
        title="WT mice",
        filters={
            "mice_model": "WT",
        },
    ),

    ah.SubsetSpec(
        label="Apc",
        title="Apc mice",
        filters={
            "mice_model": "Apc",
        },
    ),

    ah.SubsetSpec(
        label="WT_female",
        title="WT female mice",
        filters={
            "mice_model": "WT",
            "sex": "female",
        },
    ),
]

plotter.plot_all_subsets(
    subsets=subsets,
    save_png=True,
    save_pdf=True,
    show=False,
)
```

---



## 15. Full minimal example

```python
import ancombc2_heatmaps as ah


timepoints = [
    "baseline1",
    "baseline2",
    "baseline3",
    "day1",
    "day3",
    "day7",
    "day14",
]


heatmap_config = ah.HeatmapConfig(
    metadata=ah.MetadataConfig(
        sample_col="sample_name",
        timepoint_col="time_point",
        comparison_col="description_of_treatment",
        timepoints=timepoints,
        timepoint_map={
            "baseline_1": "baseline1",
            "baseline_2": "baseline2",
            "baseline_3": "baseline3",
            "day_1_post": "day1",
            "day_3_post": "day3",
            "day_7_post": "day7",
            "day_14_post": "day14",
        },
        allowed_values={
            "description_of_treatment": ["sham", "irradiated"],
        },
    ),

    comparison=ah.ComparisonConfig(
        variable_name="description_of_treatment",
    ),

    paths=ah.PathConfig(
        base_table_dir="/path/to/qza_tables",
        base_ancom_dir="/path/to/exported_ancom_results",
        metadata_path="/path/to/metadata.tsv",
        output_dir="/path/to/output",

        table_template="{timepoint}/table_{timepoint}_genus_ANCOM.qza",
        ancom_template="{timepoint}_treat_ANCOMB_exported",
    ),

    q_cutoff=0.05,
    min_sig_cells_per_taxon=1,
    split_after_timepoint="baseline3",
    cell_text_mode="relative_abundance",
)


plotter = ah.ANCOMBC2HeatmapPlotter(
    heatmap_config
)


subset = ah.SubsetSpec(
    label="",
    title="All samples",
    filters={},
)


plotter.plot_subset(
    meta_df=plotter.load_metadata(),
    subset=subset,
    save_png=True,
    save_pdf=True,
    show=True,
)
```

---

## 17. Interpreting the heatmap

The heatmap shows ANCOM-BC2 log fold changes for significant taxa across timepoints.

General interpretation:

| Visual element | Meaning |
|---|---|
| Red cells | Higher abundance in the positive comparison group |
| Blue cells | Higher abundance in the reference group |
| White cells | Not significant or missing |
| Cell numbers | Depending on `cell_text_mode`, relative abundance or LFC |
| Row height | Scaled by mean relative abundance across significant cells |
| Top labels | Sample counts per group and timepoint |
| Vertical split line | Optional separation between baseline and later timepoints |

The exact direction depends on the ANCOM-BC2 effect column.  
For example, if the effect column is:

```text
description_of_treatment::sham
```

then positive log fold changes correspond to higher abundance in `sham`, unless `invert_sign=True` is used.



# Trajectory Workflow
# Trajectory Workflow

The trajectory workflow is used to visualize taxon abundance changes across timepoints.

Trajectory plots combine:

- QIIME2 feature tables (`.qza`)
- metadata
- optional ANCOM-BC2 significance information

The workflow supports:

- mean or median trajectories
- interquartile range or confidence intervals
- optional individual mouse trajectories
- optional ANCOM significance annotations
- optional baseline merging

---

# 1. Import the package

```python
import ancombc2_heatmaps as ah
```

---

# 2. Define timepoints

Trajectory plots require ordered and numeric timepoints.

Example:

```python
timepoint_order = [
    "baseline1",
    "baseline2",
    "baseline3",
    "day1",
    "day3",
    "day7",
    "day14",
]
```

---

# 3. Create the trajectory configuration

Trajectory plots use `TrajectoryConfig`.

It consists of:

| Config class | Purpose |
|---|---|
| `TrajectoryMetadataConfig` | Metadata handling and timepoint mapping |
| `TrajectoryPathConfig` | Input paths |
| `TrajectoryPlotConfig` | Plot styling and statistics |

---

# 4. Metadata configuration

Example:

```python
trajectory_config = ah.TrajectoryConfig(

    metadata=ah.TrajectoryMetadataConfig(

        sample_col="sample_name",

        timepoint_col="time_point",

        mouse_col="mouse_id",

        comparison_col="description_of_treatment",

        genotype_col="mice_model",

        treatment_col="description_of_treatment",

        timepoint_order=[
            "baseline1",
            "baseline2",
            "baseline3",
            "day1",
            "day3",
            "day7",
            "day14",
        ],

        timepoint_numeric_map={
            "baseline1": -7,
            "baseline2": -4,
            "baseline3": -1,
            "day1": 1,
            "day3": 3,
            "day7": 7,
            "day14": 14,
        },

        timepoint_label_map={
            -7: "baseline1",
            -4: "baseline2",
            -1: "baseline3",
            1: "day1",
            3: "day3",
            7: "day7",
            14: "day14",
        },

        allowed_values={
            "description_of_treatment": [
                "sham",
                "irradiated",
            ],
        },
    ),
```

---

# 5. MetadataConfig options

| Option | Description |
|---|---|
| `sample_col` | Metadata column containing sample IDs |
| `timepoint_col` | Metadata column containing timepoints |
| `mouse_col` | Metadata column identifying individual animals |
| `comparison_col` | Main comparison variable |
| `genotype_col` | Optional genotype column |
| `treatment_col` | Optional treatment column |
| `timepoint_order` | Defines loading and plotting order |
| `timepoint_numeric_map` | Converts metadata timepoints into numeric positions |
| `timepoint_label_map` | Controls displayed x-axis labels |
| `allowed_values` | Optional metadata filtering |

---

# 6. Numeric timepoint mapping

Trajectory plots require numeric timepoints.

Example:

```python
timepoint_numeric_map={
    "baseline1": -7,
    "baseline2": -4,
    "baseline3": -1,
    "day1": 1,
    "day3": 3,
    "day7": 7,
    "day14": 14,
}
```

The numeric values are used for:

- trajectory positioning
- temporal ordering
- baseline merging
- significance alignment

The exact numbers do not need to represent real days, but should preserve the temporal order.

---

# 7. Path configuration

Example:

```python
paths=ah.TrajectoryPathConfig(

    metadata_path="/path/to/metadata.tsv",

    table_base="/path/to/qza_tables",

    ancom_base="/path/to/exported_ancom_results",

    table_template=(
        "{timepoint}/table_{timepoint}_{subset_label}.qza"
    ),

    ancom_template=(
        "{timepoint}/table_{timepoint}_{subset_label}_{variable_name}_ANCOMB_exported"
    ),
)
```

---

# 8. PathConfig options

| Option | Description |
|---|---|
| `metadata_path` | Metadata file |
| `table_base` | Base folder containing `.qza` feature tables |
| `ancom_base` | Base folder containing exported ANCOM-BC2 results |
| `table_template` | Relative path pattern for `.qza` tables |
| `ancom_template` | Relative path pattern for exported ANCOM folders |

Available template placeholders:

| Placeholder | Meaning |
|---|---|
| `{timepoint}` | Current timepoint |
| `{subset_label}` | Current subset label |
| `{variable_name}` | Comparison variable |

---

# 9. Plot configuration

Example:

```python
plot=ah.TrajectoryPlotConfig(

    estimator="mean",

    error_style="iqr",

    show_individual_lines=False,

    merge_baselines=False,

    y_lim="auto_fix",

    show_significance=True,

    q_cutoff=0.05,

    figsize=(12, 8),

    line_styles={
        "sham": "-",
        "irradiated": "--",
    },

    sns_style="whitegrid",

    sns_context="talk",

    y_label="relative abundance",
)
```

---

# 10. PlotConfig options

## Estimator

Controls the central statistic.

```python
estimator="mean"
```

Options:

| Value | Meaning |
|---|---|
| `"mean"` | Plot group means |
| `"median"` | Plot group medians |

---

## Error intervals

```python
error_style="iqr"
```

Options:

| Value | Meaning |
|---|---|
| `"iqr"` | Interquartile range |
| `"ci"` | Bootstrap 95% confidence interval |

---

## Individual mouse trajectories

```python
show_individual_lines=True
```

Displays one trajectory per mouse in the background.

Useful for:

- variability inspection
- outlier detection
- trajectory consistency

---

## Baseline merging

```python
merge_baselines=True
```

Merges:

```text
baseline1
baseline2
baseline3
```

into one baseline position.

Useful when baseline replicates should be visualized together.

---

## Y-axis scaling

```python
y_lim="auto_fix"
```

Options:

| Value | Meaning |
|---|---|
| `"auto_fix"` | Automatically determine y-axis limits |
| `(0, 0.05)` | Manual y-axis limits |

Example:

```python
y_lim=(0, 0.05)
```

---

## ANCOM significance

```python
show_significance=True
```

Displays:

```text
*
ns
```

above timepoints.

The significance threshold is controlled by:

```python
q_cutoff=0.05
```

---

## Figure size

```python
figsize=(12, 8)
```

Controls matplotlib figure dimensions.

---

## Line styles

```python
line_styles={
    "sham": "-",
    "irradiated": "--",
}
```

Allows custom line styles for groups.

---



# 11. Create the plotter

```python
plotter = ah.TaxonTrajectoryPlotter(
    trajectory_config
)
```

---

# 12. Create a subset

Example:

```python
subset = ah.SubsetSpec(
    label="WT",
    title="WT mice",
    filters={
        "mice_model": "WT",
    },
)
```

Example with multiple filters:

```python
subset = ah.SubsetSpec(
    label="WT_female",
    title="WT female mice",
    filters={
        "mice_model": "WT",
        "sex": "female",
    },
)
```

---

# 13. Plot a trajectory

Example:

```python
plotter.plot_taxon(
    taxon_query="g_Akkermansia",

    subset=subset,

    comparison_levels=[
        "sham",
        "irradiated",
    ],
)
```

---

# 14. Taxon queries

The trajectory workflow supports several query formats.

## Genus

```python
taxon_query="g_Akkermansia"
```

## Family

```python
taxon_query="f_Akkermansiaceae"
```

## Phylum

```python
taxon_query="p_Verrucomicrobiota"
```

## Exact normalized label

```python
taxon_query=(
    "p_Verrucomicrobiota; "
    "f_Akkermansiaceae; "
    "g_Akkermansia"
)
```

---

# 15. Comparison levels

```python
comparison_levels=[
    "sham",
    "irradiated",
]
```

Defines which groups should be plotted.

If omitted, all available groups are plotted automatically.


---

# 16. Full minimal example

```python
import ancombc2_heatmaps as ah


trajectory_config = ah.TrajectoryConfig(

    metadata=ah.TrajectoryMetadataConfig(

        sample_col="sample_name",

        timepoint_col="time_point",

        mouse_col="mouse_id",

        comparison_col="description_of_treatment",

        timepoint_order=[
            "baseline1",
            "baseline2",
            "baseline3",
            "day1",
            "day3",
            "day7",
            "day14",
        ],

        timepoint_numeric_map={
            "baseline1": -7,
            "baseline2": -4,
            "baseline3": -1,
            "day1": 1,
            "day3": 3,
            "day7": 7,
            "day14": 14,
        },

        timepoint_label_map={
            -7: "baseline1",
            -4: "baseline2",
            -1: "baseline3",
            1: "day1",
            3: "day3",
            7: "day7",
            14: "day14",
        },
    ),

    paths=ah.TrajectoryPathConfig(
        metadata_path="/path/to/metadata.tsv",

        table_base="/path/to/qza_tables",

        ancom_base="/path/to/ancom_exports",
    ),

    plot=ah.TrajectoryPlotConfig(
        estimator="mean",
        error_style="iqr",
        show_significance=True,
        merge_baselines=False,
    ),
)


plotter = ah.TaxonTrajectoryPlotter(
    trajectory_config
)


subset = ah.SubsetSpec(
    label="WT",

    title="WT mice",

    filters={
        "mice_model": "WT",
    },
)


plotter.plot_taxon(
    taxon_query="g_Akkermansia",

    subset=subset,

    comparison_levels=[
        "sham",
        "irradiated",
    ],
)
```

---

# 17. Trajectory interpretation

Trajectory plots visualize abundance changes across timepoints.

General interpretation:

| Visual element | Meaning |
|---|---|
| Colored lines | Group trajectories |
| Points | Group mean or median |
| Error intervals | IQR or confidence interval |
| Background lines | Optional individual mouse trajectories |
| `*` | ANCOM significant |
| `ns` | Not significant |

The exact interpretation depends on:

- selected estimator
- selected error interval
- ANCOM-BC2 comparison direction
- subset definition

```

## ANCOM-BC2