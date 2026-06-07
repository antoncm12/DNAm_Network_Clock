# DNAm Network Clock — Apply

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

Apply the **K = 60 PCA network-based epigenetic age clock** from Carcedo *et al.* (2026) to your own methylation data — in **Python**, **R**, or **MATLAB**.

This is a *use-the-clock* repository. To see how the clock was built and trained, see the companion analysis repository [\[link\]](https://github.com/antoncm12/dnam-network-clock-analysis).

## What it does

`predict_age(betas)` takes a methylation beta-value matrix (samples × CpGs) and returns a vector of predicted ages in years, one per sample.

Under the hood it:

1. Groups the input CpGs into the **4,595 modules** discovered by Infomap clustering on a CpG–CpG Spearman correlation network.
2. Aggregates each module to one value per sample — **median** of the available CpGs for multi-CpG modules, or the single CpG's value for singleton modules.
3. Standardises and centres the resulting `samples × 4,595` matrix using the training-time scaler.
4. Projects onto the **first 60 principal components** of the training-set module matrix.
5. Applies a **ridge-regression linear model** (α = 1000) trained on chronological age to produce the predicted age.

**Training metrics** (1,917 whole-blood samples, ages 18–94, 12 studies, 80/20 split): test R² = 0.840, test MAE = 4.44 years.

## Requirements for your input

A table with **samples in rows** and **CpG IDs in columns**, where the column names are Illumina cg-identifiers (e.g. `cg00006081`). β-values, not M-values.

- Missing CpGs and `NaN` values are tolerated. A warning is emitted if more than 5% of the clock's CpGs are absent from your input.
- Modules with no available CpGs are filled with the training-time mean for that module. This is the standard way to handle platform / coverage differences (e.g. EPIC vs 450K) but if a large fraction of the clock is missing, predictions will be biased toward the training-set mean age (~50 years).
- The clock was trained on **whole-blood** methylation. Performance on other tissues has not been characterised.

## Installation

### Python

```bash
pip install git+https://github.com/antoncm12/DNAm_Network_Clock
```

That's it. The trained weights ship inside the package; no extra downloads.

For a local development install:

```bash
git clone https://github.com/antoncm12/DNAm_Network_Clock
cd DNAm_Network_Clock
pip install -e .
```

Requires Python ≥ 3.9, `numpy ≥ 1.20`, `pandas ≥ 1.3`.

### R

No package install — just clone the repository and `source()` the script:

```r
source("path/to/DNAm_Network_Clock/R/dnam_network_clock.R")
```

The script auto-locates the `weights/` directory relative to itself. Tested on R ≥ 4.0; uses only base R.

### MATLAB

Clone and add the `matlab/` folder to your path:

```matlab
addpath('path/to/DNAm_Network_Clock/matlab')
```

The functions auto-locate the `weights/` directory relative to their own location. Tested on MATLAB R2020b+; uses only built-in toolboxes.

## Usage

### Python

```python
import pandas as pd
from dnam_network_clock import predict_age

# Your methylation matrix: samples (rows) × CpGs (columns)
betas = pd.read_csv("my_methylation.csv", index_col=0)

ages = predict_age(betas)        # numpy array, length n_samples
print(ages[:5])                  # [42.1  67.8  28.3  55.4  71.9]
```

To get coverage diagnostics alongside the predictions:

```python
ages, diag = predict_age(betas, return_diagnostics=True)
print(f"{diag['cpg_coverage']*100:.1f}% of clock CpGs present")
print(f"{diag['n_modules_filled_with_training_mean']} modules imputed")
```

To inspect the clock metadata:

```python
from dnam_network_clock import clock_info
info = clock_info()
print(info["training_metrics"])   # {'r2_test': 0.840, 'mae_test': 4.44, ...}
```

### R

```r
source("R/dnam_network_clock.R")

# Your methylation matrix: rows = samples, columns = CpG IDs
betas <- read.csv("my_methylation.csv", row.names = 1, check.names = FALSE)

ages <- predict_age(betas)       # numeric vector, length n_samples
head(ages)
```

To inspect metadata (requires the optional `jsonlite` package for structured output):

```r
clock_info()
```

### MATLAB

```matlab
addpath('matlab')

% Your methylation data — pass the beta matrix and CpG ID vector separately
T = readtable('my_methylation.csv', 'ReadRowNames', true, ...
              'VariableNamingRule', 'preserve');
cpgs  = T.Properties.VariableNames;
betas = table2array(T);

ages = predict_age(betas, cpgs);    % [n_samples × 1] column vector
disp(ages(1:5))
```

## Verifying your installation

The repository contains a deterministic reference example you can use to confirm any of the three implementations matches the Python reference (which produced the published results):

```
examples/example_input.csv      20 samples × 7,328 CpGs, synthetic
examples/expected_output.csv    Reference predictions from the Python implementation
```

Running `predict_age()` on `example_input.csv` should reproduce `expected_output.csv` to ~1e-10 absolute tolerance. (Tiny BLAS-level differences across platforms are expected; differences larger than ~1e-6 indicate a real bug.)

Quick sanity script in Python:

```python
import pandas as pd, numpy as np
from dnam_network_clock import predict_age

X = pd.read_csv("examples/example_input.csv", index_col=0)
expected = pd.read_csv("examples/expected_output.csv")
got = predict_age(X)
diff = np.abs(got - expected["predicted_age"].values).max()
print(f"max |new − expected| = {diff:.2e}  {'PASS' if diff < 1e-6 else 'FAIL'}")
```

## Repository layout

```
DNAm_Network_Clock/
├── README.md
├── LICENSE                          (MIT)
├── pyproject.toml                   Python package config
├── weights/                         Single source of truth — used by all three languages
│   ├── module_definitions.csv       CpG → Module → Strategy  (7,328 rows)
│   ├── module_params.csv            Per-module scaler + PCA mean  (4,595 rows)
│   ├── pca_components_k60.csv       PC loadings  (60 × 4,595)
│   ├── ridge_k60.csv                Ridge PC coefs + intercept
│   └── metadata.json                K, module counts, training metrics
├── src/dnam_network_clock/          Python package
│   ├── __init__.py
│   ├── predict.py                   The full pipeline
│   └── weights/                     Copy of weights/ bundled with the package
├── R/
│   └── dnam_network_clock.R         source() and call predict_age()
├── matlab/
│   ├── predict_age.m                Main entry point
│   ├── load_clock_weights.m
│   └── build_module_matrix.m
└── examples/
    ├── example_input.csv            Reference input (20 samples × 7,328 CpGs)
    ├── expected_output.csv          Reference predictions
    └── example.{py,R,m}             Minimal usage scripts
```

The weight CSVs are the canonical, language-agnostic representation. The Python package also bundles a copy under `src/dnam_network_clock/weights/` so `pip install` works without external data.

## Citing the clock

If you use this clock in published work, please cite the original paper:

> Carcedo, A., et al. (2026). [Manuscript title]. *PLOS Computational Biology*. https://doi.org/[paper DOI]

And the archived code/weights:

> Carcedo, A. (2026). DNAm Network Clock — apply (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.antoncm12

## License

Code in this repository is licensed under the **MIT License** (see [LICENSE](LICENSE)).
The trained weight files in `weights/` are released under **CC BY 4.0** — you may use, redistribute, and adapt them freely with attribution.

## Issues & contact

Open a [GitHub issue](https://github.com/antoncm12/DNAm_Network_Clock/issues) for bug reports or feature requests, or email [antoncm12@antoncm12] for scientific questions.
