# ANCOM-BC2 Heatmaps and Taxon Trajectories

Reusable plotting utilities for ANCOM-BC2 log fold change heatmaps, taxon trajectory plots and boxplot trajectory plots from QIIME2 feature tables and exported ANCOM-BC2 results.

## Table of Contents

- [Summary](#summary)
- [Installation](#installation)
- [Input Data](#input-data)
- [Quick Start](#quick-start)
- [Configuration](#configuration)
- [Path Handling](#path-handling)
- [Heatmap Workflow](#heatmap-workflow)
- [Trajectory Workflow](#trajectory-workflow)
- [Boxplot Workflow](#boxplot-workflow)
- [Subsets](#subsets)
- [Taxon Queries](#taxon-queries)
- [ANCOM-BC2 Direction](#ancom-bc2-direction)
- [Troubleshooting](#troubleshooting)


---

## Summary

This package helps visualize ANCOM-BC2 results across multiple timepoints.

It supports:

- ANCOM-BC2 log fold change heatmaps
- taxon trajectory plots over time
- boxplot trajectory plots with sample points
- optional ANCOM significance labels for trajectory plots
- template-based or fully explicit file paths
- metadata-based subsets, for example WT only, Apc only, male only, female only, or all samples

The package is designed for QIIME2-based microbiome workflows where each timepoint has:

1. a QIIME2 feature table (`.qza`)
2. an exported ANCOM-BC2 result directory containing:
   - `lfc.jsonl`
   - `q.jsonl`
   - `diff.jsonl`

---

## Installation

Install the package from GitHub:

```bash
pip install git+ssh://git@github.com/jlab/ancombc2-heatmaps.git
```


The package requires QIIME2 because `.qza` feature tables are loaded through the QIIME2 API.

Recommended environment:

```bash
conda activate qiime2-amplicon-2026.1
```

To see available conda environments:

```bash
conda env list
```

To check your QIIME2 version:

```bash
qiime --version
```

If you do not have a working QIIME2 environment, you can install one from the official QIIME2 distribution YAMLs:

```bash
wget https://raw.githubusercontent.com/qiime2/distributions/dev/2026.1/amplicon/released/qiime2-amplicon-ubuntu-latest-conda.yml \
  -O qiime2-amplicon-2026.1.yml

conda env create \
  -n qiime2-amplicon-2026.1 \
  -f qiime2-amplicon-2026.1.yml
```

Import the package in Python or in a notebook:

```python
import ancombc2_heatmaps as ah

print(ah.__version__)
```

---

## Input Data

The package expects three main input types:

| File type | Format | Used for |
|---|---|---|
| Metadata table | `.tsv`, `.txt` or `.csv` | sample information, groups, timepoints and subsets |
| QIIME2 feature tables | `.qza` | relative abundance calculation |
| ANCOM-BC2 export folders | folder with `lfc.jsonl`, `q.jsonl`, `diff.jsonl` | log fold changes and significance |

---

## 1. Metadata table

The metadata file should contain at least:

```text
sample column
timepoint column
group column
```

Example metadata:

```text
sample_name	time_point	description_of_treatment	mice_model	sex
sample_1	baseline_1	sham	WT	female
sample_2	baseline_1	irradiated	WT	male
sample_3	day_1_post	sham	Apc	female
sample_4	day_1_post	irradiated	Apc	male
```

In the config, these columns are specified with:

```python
sample_col="sample_name"
timepoint_col="time_point"
group_col="description_of_treatment"
```

A subject column can be used for individual lines in trajectory plots:

```python
subject_col="host_taxon_id"
```

If there is no subject or host ID column, use:

```python
subject_col=None
```

Then the package uses the sample ID as fallback.

---

## 2. QIIME2 feature tables

Feature tables must be QIIME2 `.qza` tables that can be loaded as BIOM tables.

The tables are used to calculate relative abundances for:

- heatmap cell annotations
- trajectory plots
- boxplot trajectory plots

Example:

```text
heatmaps_genus_by_timepoint/
├── table_baseline1_genus_ANCOM.qza
├── table_baseline2_genus_ANCOM.qza
├── table_baseline3_genus_ANCOM.qza
├── table_day1_genus_ANCOM.qza
├── table_day3_genus_ANCOM.qza
├── table_day7_genus_ANCOM.qza
└── table_day14_genus_ANCOM.qza
```

This structure can be represented with:

```python
table_base="/path/to/heatmaps_genus_by_timepoint"
table_template="table_{timepoint}_genus_ANCOM.qza"
```

---

## 3. Exported ANCOM-BC2 result files

Each exported ANCOM-BC2 result folder must contain:

```text
lfc.jsonl
q.jsonl
diff.jsonl
```

Example:

```text
real_ANCOMB_BC2/
├── baseline1_treat_ANCOMB_exported/
│   ├── lfc.jsonl
│   ├── q.jsonl
│   └── diff.jsonl
├── baseline2_treat_ANCOMB_exported/
│   ├── lfc.jsonl
│   ├── q.jsonl
│   └── diff.jsonl
└── day1_treat_ANCOMB_exported/
    ├── lfc.jsonl
    ├── q.jsonl
    └── diff.jsonl
```

This structure can be represented with:

```python
ancom_base="/path/to/real_ANCOMB_BC2"
ancom_template="{timepoint}_treat_ANCOMB_exported"
```

---

## Quick Start

This example reproduces the simplified version 2.0 workflow.

```python
import ancombc2_heatmaps as ah

timepoints = [
    "baseline1", "baseline2", "baseline3",
    "day1", "day3", "day7", "day14",
]

timepoint_map = {
    "baseline_1": "baseline1",
    "baseline_2": "baseline2",
    "baseline_3": "baseline3",
    "day_1_post": "day1",
    "day_3_post": "day3",
    "day_7_post": "day7",
    "day_14_post": "day14",
}

config = ah.ANCOMConfig(
    metadata_path="/path/to/metadata.tsv",
    output_dir="/path/to/output",

    table_base="/path/to/heatmaps_genus_by_timepoint",
    ancom_base="/path/to/heatmaps_genus_by_timepoint/real_ANCOMB_BC2",

    sample_col="sample_name",
    timepoint_col="time_point",
    group_col="description_of_treatment",
    subject_col=None,

    timepoints=timepoints,
    timepoint_map=timepoint_map,

    allowed_values={
        "description_of_treatment": ["sham", "irradiated"],
    },

    table_template="table_{timepoint}_genus_ANCOM.qza",
    ancom_template="{timepoint}_treat_ANCOMB_exported",

    variable_name="description_of_treatment",

    split_after_timepoint="baseline3",

    taxa=[
        "g_Akkermansia",
    ],
)

workflow = ah.ANCOMWorkflow(config)

subset = ah.Subset(
    label="",
    title="All samples",
    filters={},
)
```

Before plotting, check whether all paths are resolved correctly:

```python
workflow.check(subset)
```

Then create plots:

```python
workflow.heatmap(subset)

workflow.trajectory("g_Akkermansia", subset)
workflow.boxplot("g_Akkermansia", subset)
```

Because `g_Akkermansia` is already listed in `config.taxa`, you can also use:

```python
workflow.trajectories(subset)
workflow.boxplots(subset)
```

---

## Configuration

The package uses one central config class:

```python
ah.ANCOMConfig(...)
```

Important fields:

| Field | Meaning |
|---|---|
| `metadata_path` | path to metadata file |
| `output_dir` | directory for saved figures |
| `sample_col` | metadata column containing sample IDs |
| `timepoint_col` | metadata column containing timepoints |
| `group_col` | metadata column used for comparison groups |
| `subject_col` | optional mouse/subject ID column for individual lines |
| `timepoints` | ordered list of timepoints to plot |
| `timepoint_map` | optional mapping from metadata names to plot names |
| `allowed_values` | optional metadata filtering before plotting |
| `table_base` | base directory for QIIME2 tables |
| `table_template` | filename or relative path template for QIIME2 tables |
| `ancom_base` | base directory for ANCOM-BC2 export folders |
| `ancom_template` | folder template for ANCOM-BC2 exports |
| `variable_name` | ANCOM-BC2 variable name / effect prefix |
| `taxa` | default taxa for `workflow.trajectories()` and `workflow.boxplots()` |

Additional plotting options are grouped in `ah.PlotConfig`.

Example:

```python
config = ah.ANCOMConfig(
    ...,
    plot=ah.PlotConfig(
        save_png=False,
        save_pdf=False,
        show=True,
        estimator="mean",
        error_style="iqr",
        show_significance=True,
        merge_baselines=False,
    ),
)
```

---

## Path Handling

The package supports two ways to define file paths.

---

### 1. Template mode

Use this when your files follow a regular naming pattern.

```python
config = ah.ANCOMConfig(
    ...,
    table_base="/path/to/tables",
    ancom_base="/path/to/ancom_exports",

    table_template="{timepoint}/table_{timepoint}_{subset_label}.qza",
    ancom_template="{timepoint}/table_{timepoint}_{subset_label}_ANCOMB_exported",
)
```

For example, for `timepoint="day1"` and `subset.label="WT_genus_ANCOM"`, this resolves to:

```text
/path/to/tables/day1/table_day1_WT_genus_ANCOM.qza
/path/to/ancom_exports/day1/table_day1_WT_genus_ANCOM_ANCOMB_exported
```

Templates can use:

```text
{timepoint}
{subset_label}
{variable_name}
```

---

### 2. Fixed path mode

Use this when your filenames are irregular or when you want to specify every file manually.

```python
config = ah.ANCOMConfig(
    ...,
    table_paths={
        "baseline1": "/path/to/baseline1_table.qza",
        "day1": "/path/to/day1_table.qza",
        "day3": "/path/to/day3_table.qza",
    },
    ancom_paths={
        "baseline1": "/path/to/baseline1_ANCOMB_exported",
        "day1": "/path/to/day1_ANCOMB_exported",
        "day3": "/path/to/day3_ANCOMB_exported",
    },
)
```

If `table_paths` is given, it overrides `table_base` and `table_template`.

If `ancom_paths` is given, it overrides `ancom_base` and `ancom_template`.

---

## Heatmap Workflow

Create one heatmap:

```python
workflow.heatmap(subset)
```

Disable saving:

```python
workflow.heatmap(
    subset,
    save_png=False,
    save_pdf=False,
    show=True,
)
```

Create heatmaps for multiple subsets:

```python
subsets = [
    ah.Subset(label="", title="All samples", filters={}),
    ah.Subset(label="WT_genus_ANCOM", title="WT", filters={"mice_model": "WT"}),
    ah.Subset(label="Apc_genus_ANCOM", title="Apc", filters={"mice_model": "Apc"}),
]

workflow.heatmaps(subsets)
```

The heatmap shows ANCOM-BC2 log fold changes across timepoints.

By default:

- only significant cells are colored
- taxa are filtered by `min_sig_cells_per_taxon`
- cell text shows mean relative abundance
- row height is scaled by relative abundance
- the top axis shows group sample counts per timepoint

Relevant config options:

```python
config = ah.ANCOMConfig(
    ...,
    min_sig_cells_per_taxon=1,
    split_after_timepoint="baseline3",
    remove_empty_rows_after_masking=True,
    plot=ah.PlotConfig(
        cell_text_mode="relative_abundance",  # "relative_abundance", "lfc", or "none"
        heatmap_cmap="RdBu_r",
        highlight_taxa=["Akkermansia", "Muribaculaceae"],
    ),
)
```

---

## Trajectory Workflow

Plot one taxon:

```python
workflow.trajectory("g_Akkermansia", subset)
```

Plot all taxa listed in `config.taxa`:

```python
workflow.trajectories(subset)
```

Temporarily plot a different list of taxa:

```python
workflow.trajectories(
    subset,
    taxa=[
        "g_Akkermansia",
        "f_Muribaculaceae",
        "p_Verrucomicrobiota",
    ],
)
```

Choose comparison levels manually:

```python
workflow.trajectory(
    "g_Akkermansia",
    subset,
    comparison_levels=["sham", "irradiated"],
)
```

Relevant plot settings:

```python
config = ah.ANCOMConfig(
    ...,
    plot=ah.PlotConfig(
        estimator="mean",            # "mean" or "median"
        error_style="iqr",           # "iqr" or "ci"
        show_individual_lines=False,
        show_significance=True,
        y_label="relative abundance",
    ),
)
```

If `show_significance=True`, labels are added above the trajectory:

```text
*  = ANCOM significant
ns = not significant
```

---

## Boxplot Workflow

Plot one taxon as grouped boxplots over time:

```python
workflow.boxplot("g_Akkermansia", subset)
```

Plot all taxa listed in `config.taxa`:

```python
workflow.boxplots(subset)
```

Disable the trend line:

```python
workflow.boxplot(
    "g_Akkermansia",
    subset,
    show_trend=False,
)
```

Change the polynomial degree of the trend line:

```python
workflow.boxplot(
    "g_Akkermansia",
    subset,
    show_trend=True,
    trend_order=3,
)
```

Typical values:

```text
trend_order=1  linear trend
trend_order=2  quadratic trend
trend_order=3  cubic trend
```

For multiple taxa:

```python
workflow.boxplots(
    subset,
    show_trend=True,
    trend_order=2,
)
```

The boxplot workflow shows:

- boxplots per timepoint and group
- individual sample points
- optional approximate polynomial trend line
- optional ANCOM significance labels

---

## Subsets

Subsets are defined with:

```python
ah.Subset(...)
```

A subset has:

| Field | Meaning |
|---|---|
| `label` | used in file templates |
| `title` | shown in plot titles |
| `filters` | metadata filters applied before plotting |

All samples:

```python
subset_all = ah.Subset(
    label="",
    title="All samples",
    filters={},
)
```

WT only:

```python
subset_wt = ah.Subset(
    label="WT_genus_ANCOM",
    title="WT",
    filters={"mice_model": "WT"},
)
```

Female WT only:

```python
subset_wt_female = ah.Subset(
    label="WT_female_genus_ANCOM",
    title="WT female",
    filters={
        "mice_model": "WT",
        "sex": "female",
    },
)
```

Multiple accepted values are also possible:

```python
subset_baseline_groups = ah.Subset(
    label="WT_genus_ANCOM",
    title="WT, both sexes",
    filters={
        "mice_model": "WT",
        "sex": ["female", "male"],
    },
)
```

Important: `subset.label` is also used to build file paths when the template contains `{subset_label}`.

---

## Taxon Queries

Taxa can be selected with short query strings.

Examples:

```python
"g_Akkermansia"
"f_Akkermansiaceae"
"p_Verrucomicrobiota"
```

Supported prefixes:

| Prefix | Rank |
|---|---|
| `g_` | genus |
| `f_` | family |
| `p_` | phylum |
| `o_` | order |
| `c_` | class |
| `k_` | kingdom |
| `d_` | domain |

You can also use the normalized full taxon label:

```python
"p__Verrucomicrobiota; f__Akkermansiaceae; g__Akkermansia"
```


---

## ANCOM-BC2 Direction

The package uses the ANCOM-BC2 log fold change direction directly from the selected effect column.

Example effect column:

```text
description_of_treatment::sham
```

This means:

```text
positive LFC = higher in sham
negative LFC = lower in sham / higher in the reference group
```

The heatmap colorbar labels are inferred from the effect column and metadata where possible.

You can set labels manually:

```python
config = ah.ANCOMConfig(
    ...,
    positive_class="sham",
    negative_class="irradiated",
)
```

If you explicitly want to reverse the LFC direction:

```python
config = ah.ANCOMConfig(
    ...,
    invert_sign=True,
)
```

Use `invert_sign=True` only if you are sure that you want to flip the ANCOM-BC2 LFC values.

---

## Troubleshooting

### 1. Check resolved paths

Before plotting a new dataset, run:

```python
workflow.check(subset)
```

This returns a table with:

```text
timepoint
table_path
table_exists
ancom_path
ancom_exists
```

If `table_exists` or `ancom_exists` is `False`, your base paths or templates do not match your actual file structure.

---

### 2. No relative abundance data

Message:

```text
[NO HEATMAP] No relative abundance data
```

Common causes:

- QIIME2 table paths are wrong
- metadata sample IDs do not match table sample IDs
- timepoint names in metadata do not match `config.timepoints`
- subset filters remove all samples

Check metadata and table overlap:

```python
meta = workflow.data.filter_metadata(subset)
rel = workflow.data.relative_abundance("baseline1", subset)

metadata_samples = set(meta.loc[meta["time_point"] == "baseline1", "sample_name"])
table_samples = set(rel.columns)

len(metadata_samples & table_samples)
```

If the overlap is `0`, the sample names do not match.

---

### 3. Missing metadata columns

Message:

```text
ValueError: Missing metadata columns: [...]
```

The column name in your config does not exist in the metadata file.

---

### 4. No LFC data found

Message:

```text
[NO HEATMAP] No ANCOM LFC data
```

Common causes:

- ANCOM export paths are wrong
- `lfc.jsonl`, `q.jsonl` or `diff.jsonl` is missing
- `variable_name` does not match the effect column prefix

Check one export folder manually:

```bash
ls /path/to/ANCOMB_exported
```

Expected:

```text
lfc.jsonl
q.jsonl
diff.jsonl
```

If the ANCOM-BC2 column is for example:

```text
description_of_treatment::sham
```

then use:

```python
variable_name="description_of_treatment"
```

or specify exactly:

```python
effect_column="description_of_treatment::sham"
```

---


## Minimal complete example

```python
import ancombc2_heatmaps as ah

timepoints = [
    "baseline1", "baseline2", "baseline3",
    "day1", "day3", "day7", "day14",
]

config = ah.ANCOMConfig(
    metadata_path="/path/to/metadata.tsv",
    output_dir="/path/to/output",

    table_base="/path/to/tables",
    ancom_base="/path/to/ancom",

    sample_col="sample_name",
    timepoint_col="time_point",
    group_col="description_of_treatment",
    subject_col=None,

    timepoints=timepoints,

    table_template="table_{timepoint}_genus_ANCOM.qza",
    ancom_template="{timepoint}_treat_ANCOMB_exported",

    variable_name="description_of_treatment",

    taxa=["g_Akkermansia"],

    plot=ah.PlotConfig(
        save_png=False,
        save_pdf=False,
        show=True,
    ),
)

workflow = ah.ANCOMWorkflow(config)

subset = ah.Subset(
    label="",
    title="All samples",
    filters={},
)

workflow.check(subset)

workflow.heatmap(subset)
workflow.trajectories(subset)
workflow.boxplots(subset)
```
