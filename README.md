# DNAm Network Clock

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Apply the K=60 PCA network based epigenetic age clock from Carcedo et al. (2026) to your own methylation data, in Python, R, or MATLAB.

This repository is deliberately simple. There is no package to install. For each language there is a single file that defines a `predict_age` function, plus a short example. You copy the one file for your language into your project (or keep it next to the `weights` folder) and call it. The `weights` folder is shared by all three languages and is the single source of truth.

## What it does

`predict_age(betas)` takes a methylation beta value matrix, samples in rows and CpG IDs in columns, and returns one predicted age in years per sample.

Under the hood it groups the input CpGs into the 4,595 modules discovered by Infomap clustering on a CpG to CpG Spearman correlation network, aggregates each module to one value per sample (the median of the available CpGs for multi CpG modules, the single CpG's value for singleton modules), standardises and centres that samples by 4,595 matrix with the training time scaler, projects it onto the first 60 principal components, and applies a ridge regression linear model (alpha = 1000) to produce the age.

Training metrics, 1,917 whole blood samples, ages 18 to 94, 12 studies.

## Your input

A table with samples in rows and CpG IDs in columns, where the column names are Illumina cg identifiers (for example `cg00006081`). These are beta values, not M values.

Missing CpGs and NaN values are tolerated. A warning is emitted if more than 5% of the clock's CpGs are absent from your input. Modules with no available CpGs are filled with the training time mean for that module, which is the standard way to handle platform and coverage differences (for example EPIC versus 450K), but if a large fraction of the clock is missing the predictions will be biased toward the training set mean age (around 50 years). The clock was trained on whole blood methylation, performance on other tissues has not been characterised.

## Usage

### Python

Requires numpy and pandas (Python 3.9 or newer). Copy `python/dnam_network_clock.py` into your project, or run from a folder where the `weights` directory is reachable.

```python
import pandas as pd
from dnam_network_clock import predict_age

# Your methylation matrix, samples (rows) x CpGs (columns)
betas = pd.read_csv("my_methylation.csv", index_col=0)

ages = predict_age(betas)        # numpy array, length n_samples
print(ages[:5])

# If the weights live somewhere specific, point to them
ages = predict_age(betas, weights_dir="path/to/weights")

# Coverage diagnostics alongside the predictions
ages, diag = predict_age(betas, return_diagnostics=True)
print(f"{diag['cpg_coverage']*100:.1f}% of clock CpGs present")
```

### R

Uses only base R (tested on R 4.0 or newer). Source the one file, then call `predict_age`.

```r
source("path/to/dnam_network_clock.R")

# Your methylation matrix, rows = samples, columns = CpG IDs
betas <- read.csv("my_methylation.csv", row.names = 1, check.names = FALSE)

ages <- predict_age(betas)       # numeric vector, length n_samples
head(ages)

# Point to a specific weights folder if needed
ages <- predict_age(betas, weights_dir = "path/to/weights")
```

If the script cannot find the weights automatically, set `options(dnam_clock.weights_dir = "path/to/weights")` or pass `weights_dir=`.

### MATLAB

Uses only built in functions (tested on R2020b or newer). Put `predict_age.m` on your path. Unlike Python and R, you pass the beta matrix and the CpG ID vector separately.

```matlab
addpath('path/to/matlab')

T = readtable('my_methylation.csv', 'ReadRowNames', true, ...
              'VariableNamingRule', 'preserve');
cpgs  = T.Properties.VariableNames;
betas = table2array(T);

ages = predict_age(betas, cpgs);          % [n_samples x 1] column vector
disp(ages(1:5))

% Point to a specific weights folder if needed
ages = predict_age(betas, cpgs, 'path/to/weights');
```

## Checking it works

The `examples` folder holds a deterministic reference, `example_input.csv` (20 synthetic samples by 7,328 CpGs) and `expected_output.csv` (the predictions from the Python implementation, which produced the published results). Running `predict_age` on the input reproduces the expected output to roughly 1e-10. Differences larger than about 1e-6 indicate a real bug. Each language folder has a ready to run example that does this comparison for you:

```
python python/example.py
Rscript R/example.R
matlab -batch example          % or, from the MATLAB prompt, run: example
```

## Repository layout

```
dnam-network-clock/
  README.md
  LICENSE                          MIT
  weights/                         shared by all three languages
    module_definitions.csv         CpG, Module, Strategy (7,328 rows)
    module_params.csv              per module scaler and PCA mean (4,595 rows)
    pca_components_k60.csv          PC loadings (60 by 4,595)
    ridge_k60.csv                  ridge PC coefficients and intercept
    metadata.json                  K, module counts, training metrics
  examples/
    example_input.csv              reference input (20 samples by 7,328 CpGs)
    expected_output.csv            reference predictions
  python/
    dnam_network_clock.py          the function, copy this into your project
    example.py
  R/
    dnam_network_clock.R           the function, source this
    example.R
  matlab/
    predict_age.m                  the function, put this on your path
    example.m
```

## Citing the clock

If you use this clock in published work, please cite the paper:

> Carcedo, A., et al. (2026). [Manuscript title]. PLOS Computational Biology. https://doi.org/[paper DOI]

## License

All code in this repository is released under the MIT License, see [LICENSE](LICENSE). The trained weight files in `weights/` are released under CC BY 4.0, you may use, redistribute, and adapt them freely with attribution.
