STEFAN
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
- [Significance Filtering](#significance-filtering)
- [Trajectory Workflow](#trajectory-workflow)
- [Boxplot Trajectory Workflow](#boxplot-trajectory-workflow)
- [Subsets](#subsets)
- [Taxon Queries](#taxon-queries)
- [ANCOM-BC2 Direction](#ancom-bc2-direction)
- [Troubleshooting](#troubleshooting)
- [Minimal complete example](#minimal-complete-example)

---

## Summary

This package helps visualize ANCOM-BC2 results across multiple timepoints.

It supports:

- ANCOM-BC2 log fold change heatmaps
- taxon trajectory plots over time
- boxplot trajectory plots over time

The package is designed for QIIME2-based microbiome workflows where each timepoint has:

1. a QIIME2 feature table (`.qza`)
2. an exported ANCOM-BC2 result directory containing:
   - `lfc.jsonl`
   - `p.jsonl`
   - `q.jsonl`
   - `diff.jsonl`

---

<img width="1127" height="1127" alt="grafik" src="https://github.com/user-attachments/assets/866a5592-3334-408a-9873-31eee1e7df3d" />

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
| ANCOM-BC2 export folders | folder with `lfc.jsonl`, `p.jsonl`, `q.jsonl`, `diff.jsonl` | log fold changes and significance |

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

To learn more about setting up the correct paths you should look at [Path Handling](#path-handling).

---

## 3. Exported ANCOM-BC2 result files

Each exported ANCOM-BC2 result folder must contain:

```text
lfc.jsonl
p.jsonl
q.jsonl
diff.jsonl
```

The package uses:

| File | Used for |
|---|---|
| `lfc.jsonl` | ANCOM-BC2 log fold changes |
| `p.jsonl` | optional unadjusted p-value filtering |
| `q.jsonl` | FDR-adjusted q-value filtering |
| `diff.jsonl` | optional ANCOM-BC2 differential-abundance flag |

Example:

```text
real_ANCOMB_BC2/
├── baseline1_treat_ANCOMB_exported/
│   ├── lfc.jsonl
│   ├── p.jsonl
│   ├── q.jsonl
│   └── diff.jsonl
├── baseline2_treat_ANCOMB_exported/
│   ├── lfc.jsonl
│   ├── p.jsonl
│   ├── q.jsonl
│   └── diff.jsonl
└── day1_treat_ANCOMB_exported/
    ├── lfc.jsonl
    ├── p.jsonl
    ├── q.jsonl
    └── diff.jsonl
```

This structure can be represented with:

```python
ancom_base="/path/to/real_ANCOMB_BC2"
ancom_template="{timepoint}_treat_ANCOMB_exported"
```

To learn more about setting up the correct paths you should look at [Path Handling](#path-handling).

---

## Quick Start

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

    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=False,

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
| `timepoint_numeric_map` | optional mapping from timepoint labels to numeric values |
| `timepoint_label_map` | optional mapping from numeric timepoints to x-axis labels |
| `allowed_values` | optional metadata filtering before plotting |
| `table_base` | base directory for QIIME2 tables |
| `table_template` | filename or relative path template for QIIME2 tables |
| `ancom_base` | base directory for ANCOM-BC2 export folders |
| `ancom_template` | folder template for ANCOM-BC2 exports |
| `table_paths` | optional explicit paths for QIIME2 feature tables |
| `ancom_paths` | optional explicit paths for ANCOM-BC2 export folders |
| `lfc_filename` | filename of the LFC export, default: `lfc.jsonl` |
| `p_filename` | filename of the unadjusted p-value export, default: `p.jsonl` |
| `q_filename` | filename of the FDR-adjusted q-value export, default: `q.jsonl` |
| `diff_filename` | filename of the ANCOM-BC2 differential-abundance flag export, default: `diff.jsonl` |
| `variable_name` | ANCOM-BC2 variable name / effect prefix |
| `effect_column` | exact ANCOM-BC2 effect column if automatic detection is ambiguous |
| `positive_class` | manual label for the positive LFC direction |
| `negative_class` | manual label for the negative/reference LFC direction |
| `invert_sign` | reverses the LFC sign if set to `True` |
| `significance_metric` | significance file used for heatmap filtering: `"q"` or `"p"` |
| `q_cutoff` | cutoff used for the selected significance metric, default: `0.05` |
| `require_diff_for_significance` | whether `diff == True` is additionally required |
| `min_sig_cells_per_taxon` | minimum number of significant cells required to keep a taxon in the heatmap |
| `split_after_timepoint` | optional vertical separator after a selected timepoint |
| `remove_empty_rows_after_masking` | remove taxa without plottable significant cells |
| `taxa` | default taxa for `workflow.trajectories()` and `workflow.boxplots()` |

Additional plotting options are grouped in `ah.PlotConfig`.

Example:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=False,
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

- only cells passing the selected significance filter are colored
- taxa are filtered by `min_sig_cells_per_taxon`
- cell text shows mean relative abundance
- row height is scaled by relative abundance
- the top axis shows group sample counts per timepoint

Relevant config options:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=False,
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

### Vertical separator after a baseline

Use `split_after_timepoint` to draw a vertical line after a selected timepoint:

```python
config = ah.ANCOMConfig(
    ...,
    split_after_timepoint="baseline3",
)
```

For example, if your heatmap timepoints are:

```python
timepoints = ["baseline1", "day 1", "day 3", "day 7"]
```

you can draw the separator after `baseline1` with:

```python
split_after_timepoint="baseline1"
```

---

## Significance Filtering

Heatmap cells are filtered using ANCOM-BC2 significance values.

The package supports two significance metrics:

| Option | File used | Meaning |
|---|---|---|
| `significance_metric="q"` | `q.jsonl` | FDR-adjusted q-values |
| `significance_metric="p"` | `p.jsonl` | unadjusted p-values |

The cutoff is set with:

```python
q_cutoff=0.05
```

The name `q_cutoff` is kept for backwards compatibility. It is used as the cutoff for whichever significance metric is selected.

---

### Recommended q-value filtering

For final interpretation, FDR-adjusted q-values are recommended:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=False,
)
```

This displays heatmap cells where:

```text
q < 0.05
```

The q-values are read from:

```text
q.jsonl
```

---

### Exploratory p-value filtering

For exploratory plots, unadjusted p-values can be used:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="p",
    q_cutoff=0.05,
    require_diff_for_significance=False,
)
```

This displays heatmap cells where:

```text
p < 0.05
```

The p-values are read from:

```text
p.jsonl
```

When using:

```python
significance_metric="p"
```

the package raises a warning because unadjusted p-values are not corrected for multiple testing and may increase the number of false positives.

Use this option for exploratory visualization only.

---

### Optional use of `diff.jsonl`

ANCOM-BC2 also exports `diff.jsonl`, which contains a `True`/`False` flag for differential abundance.

If you want a stricter filter, require both the selected significance value and `diff == True`:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=True,
)
```

This displays heatmap cells where:

```text
q < 0.05 and diff == True
```

For unadjusted p-values, this stricter version would be:

```python
config = ah.ANCOMConfig(
    ...,
    significance_metric="p",
    q_cutoff=0.05,
    require_diff_for_significance=True,
)
```

which displays cells where:

```text
p < 0.05 and diff == True
```

For final interpretation, FDR-adjusted q-values are recommended. Unadjusted p-values should be clearly marked as exploratory.

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

The significance labels use the same significance filtering settings as the heatmap, for example:

```python
significance_metric="q"
q_cutoff=0.05
require_diff_for_significance=False
```

---

<img width="1166" height="766" alt="grafik" src="https://github.com/user-attachments/assets/6c68b63e-d9e5-478b-be2d-ca24a6573b17" />

---

## Boxplot Trajectory Workflow

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

The significance labels use the same significance filtering settings as the heatmap.

---

<img width="1166" height="766" alt="grafik" src="https://github.com/user-attachments/assets/f431cb66-ea8f-4290-b878-0bf5ffcbe279" />

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
- `lfc.jsonl`, `p.jsonl`, `q.jsonl` or `diff.jsonl` is missing
- `variable_name` does not match the effect column prefix

Check one export folder manually:

```bash
ls /path/to/ANCOMB_exported
```

Expected:

```text
lfc.jsonl
p.jsonl
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

### 5. No taxa after significance filter

Message:

```text
[NO HEATMAP] No taxa after significance filter
```

Common causes:

- no taxon passes the selected significance filter
- `significance_metric` is set to `"q"` but no q-values are below the cutoff
- `require_diff_for_significance=True` removes taxa where `diff == False`
- `effect_column` points to the wrong ANCOM-BC2 effect

For a less strict q-value filter, use:

```python
significance_metric="q"
q_cutoff=0.05
require_diff_for_significance=False
```

For exploratory unadjusted p-value plots, use:

```python
significance_metric="p"
q_cutoff=0.05
require_diff_for_significance=False
```

Unadjusted p-values are not corrected for multiple testing and should be used for exploratory visualization only.

---

### 6. Suppress seaborn FutureWarnings

Depending on the pandas/seaborn versions in your environment, seaborn may print `FutureWarning` messages during plotting, for example:

> [!WARNING]
> `FutureWarning: use_inf_as_na option is deprecated and will be removed in a future version.`

This warning comes from seaborn/pandas compatibility and does not usually affect the generated figures.

To hide this warning in a notebook, add this at the top of the notebook before creating plots:

```python
import warnings

warnings.filterwarnings(
    "ignore",
    category=FutureWarning,
    module="seaborn",
)
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

    significance_metric="q",
    q_cutoff=0.05,
    require_diff_for_significance=False,

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
