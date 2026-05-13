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
- [Examples](#examples)



## Summary

blablabla

## Installation

The installation of the package can be done from GitHub via git in the console:

[!FIXME] # maybe a change after making it public
```
pip install git+ssh://git@github.com/jlab/ancombc2-heatmaps.git  
```



Please be aware that you have to activate the correct conda environment with an Qiime2 Version that is not too old.
I recommend "qiime2-amplicon-2026.1".

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

After installation you can check the version with

```
!qiime --version
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


## Input Data















