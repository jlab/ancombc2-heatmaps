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



## Trajectory Workflow


## Boxplot Workflow


## ANCOM-BC2