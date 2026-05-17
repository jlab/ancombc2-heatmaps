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
.

---

After installation you can check the version with

```
qiime --version
```

and import it for example in your notebook with

```
import ancombc2_heatmaps as ah
```
.

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

The important columns are:

```text
taxon
<effect_column>
```

For example:

```text
taxon
description_of_treatment::sham
```

or:

```text
taxon
sex::female
```

The effect column is defined in the comparison configuration:

```python
comparison=ah.ComparisonConfig(
    variable_name="description_of_treatment",
    positive_class="irradiated",
    negative_class="sham",
    effect_column="description_of_treatment::sham",
    invert_sign=True,
)
```

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

The heatmap workflow consists of four main steps:

1. Create a heatmap configuration
2. Create a heatmap plotter
3. Define subsets
4. Generate the heatmap

---

### Minimal Heatmap Example

```python
import ancombc2_heatmaps as ah

heatmap_config = ah.HeatmapConfig(

    metadata=ah.MetadataConfig(

        sample_col="#SampleID",

        timepoint_col="time_point",

        comparison_col="description_of_treatment",

        timepoints=[
            "baseline1",
            "day1",
            "day3",
            "day14",
        ],
    ),

    comparison=ah.ComparisonConfig(

        variable_name="description_of_treatment",

        positive_class="irradiated",

        negative_class="sham",
    ),

    paths=ah.PathConfig(

        base_table_dir="collapsed_tables/genus/by_treatment",

        base_ancom_dir="ancombc2_results/genus/by_treatment",

        metadata_path="metadata.tsv",

        output_dir="plots",
    ),
)
```

---

### Create Heatmap Plotter

```python
plotter = ah.ANCOMBC2HeatmapPlotter(
    heatmap_config
)
```

---

### Define a Subset

Subsets are used for subgroup analyses.

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

Additional examples:

```python
filters={
    "sex": "female"
}
```

```python
filters={
    "mice_model": "WT",
    "sex": "male",
}
```

```python
filters={
    "description_of_treatment": [
        "sham",
        "irradiated",
    ]
}
```

---

### Generate Heatmap

```python
plotter.plot_subset(

    meta_df=plotter.load_metadata(),

    subset=subset,
)
```

---

## Heatmap Features

### Relative abundance annotations

```python
cell_text_mode="relative_abundance"
```

Displays mean relative abundance values inside heatmap cells.

---

### Log fold change annotations

```python
cell_text_mode="lfc"
```

Displays ANCOM-BC2 log fold changes inside heatmap cells.

---

### Disable cell text

```python
cell_text_mode="none"
```

---

### Significance filtering

Only keep taxa with a minimum number of significant cells.

```python
min_sig_cells_per_taxon=2
```

---

### Highlight taxa

```python
style=ah.HeatmapStyleConfig(

    highlight_taxa=[
        "Akkermansia",
        "Muribaculaceae",
    ]
)
```

Highlighted taxa are shown in bold.

---

### Split baseline and treatment timepoints

```python
split_after_timepoint="baseline3"
```

Useful for separating baseline and post-treatment phases.

---

### Heatmap styling

Example:

```python
style=ah.HeatmapStyleConfig(

    figure_min_width=10,

    figure_max_width=16,

    base_cell_height=0.25,

    height_scale=6.0,

    heatmap_cmap="RdBu_r",
)
```

The package automatically rescales:

- row heights
- figure height
- figure width

based on abundance and number of taxa.

---

## Trajectory Workflow

Trajectory plots visualize abundance changes over time.

Supported:

- mean trajectories
- median trajectories
- IQR visualization
- confidence intervals
- significance annotations
- individual sample trajectories

---

### Minimal Trajectory Example

```python
trajectory_config = ah.TrajectoryConfig(

    metadata=ah.TrajectoryMetadataConfig(

        sample_col="#SampleID",

        timepoint_col="time_point",

        mouse_col="mouse_id",

        comparison_col="description_of_treatment",

        timepoint_order=[
            "baseline1",
            "day1",
            "day3",
            "day14",
        ],

        timepoint_numeric_map={
            "baseline1": -1,
            "day1": 1,
            "day3": 3,
            "day14": 14,
        },

        timepoint_label_map={
            -1: "baseline",
            1: "day1",
            3: "day3",
            14: "day14",
        },
    ),

    paths=ah.TrajectoryPathConfig(

        metadata_path="metadata.tsv",

        table_base="collapsed_tables/genus/by_treatment",

        ancom_base="ancombc2_results/genus/by_treatment",
    ),
)
```

---

### Create Trajectory Plotter

```python
plotter = ah.TaxonTrajectoryPlotter(
    trajectory_config
)
```

---

### Generate Trajectory Plot

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

## Trajectory Features

### Mean trajectories

```python
estimator="mean"
```

---

### Median trajectories

```python
estimator="median"
```

---

### IQR visualization

```python
error_style="iqr"
```

---

### Bootstrap confidence intervals

```python
error_style="ci"
```

---

### Individual sample trajectories

```python
show_individual_lines=True
```

---

### Merge baselines

```python
merge_baselines=True
```

Merges multiple baseline measurements into one baseline group.

---

### Significance labels

```python
show_significance=True
```

Displays:

```text
*  = significant
ns = not significant
```

---

## Boxplot Workflow

The package also supports boxplot-based trajectories.

Features:

- boxplots
- sample points
- trend lines
- significance labels
- whisker-based scaling

---

### Generate Boxplot Trajectory

```python
workflow.plot_boxplot_trajectory(

    subset=subset,

    taxon_query="g_Akkermansia",

    comparison_levels=[
        "sham",
        "irradiated",
    ],

    show_trend=True,

    trend_order=3,
)
```

---

### Trend Lines

Trend lines are polynomial approximations generated using:

```python
numpy.polyfit()
```

The polynomial order can be adjusted:

```python
trend_order=2
```

```python
trend_order=4
```

Higher values create more flexible curves.

---

### Workflow API

The package provides a workflow wrapper which combines:

- heatmaps
- trajectory plots
- boxplot trajectories

into one interface.

---

### Create Workflow

```python
workflow = ah.PlotWorkflow(

    heatmap_config=heatmap_config,

    trajectory_config=trajectory_config,
)
```

---

### Heatmap

```python
workflow.plot_heatmap(
    subset=subset,
)
```

---

### Trajectory

```python
workflow.plot_trajectory(

    subset=subset,

    taxon_query="g_Akkermansia",
)
```

---

### Boxplot trajectory

```python
workflow.plot_boxplot_trajectory(

    subset=subset,

    taxon_query="g_Akkermansia",
)
```

---

### Combined plotting

```python
workflow.plot_heatmap_and_trajectory(

    subset=subset,

    taxon_query="g_Akkermansia",
)
```

```python
workflow.plot_heatmap_and_boxplot_trajectory(

    subset=subset,

    taxon_query="g_Akkermansia",
)
```

---

### Taxon Queries

The package supports shorthand taxon queries.

Examples:

| Rank | Example |
|---|---|
| Phylum | `p_Bacteroidota` |
| Family | `f_Muribaculaceae` |
| Genus | `g_Akkermansia` |

---

### Exact labels

Full exact labels are also supported.

```python
taxon_query="p__Bacteroidota; f__Muribaculaceae; g__Muribaculum"
```

---

### Automatic Query Discovery

Available taxa queries can be listed automatically.

```python
plotter.list_available_queries(
    qza_fp
)
```

Returns:

```python
{
    "family_queries": [...],
    "phylum_queries": [...],
    "genus_queries": [...],
}
```

---

### Output

Plots can:

- be displayed interactively
- be saved as PNG
- be saved as PDF

Example:

```python
plotter.plot_subset(

    meta_df=plotter.load_metadata(),

    subset=subset,

    save_png=True,

    save_pdf=True,

    show=True,
)
```








