"""DNAm Network Clock, age prediction from methylation beta values (Python).

This is a single, self contained file. There is nothing to install. Either
copy it into your own project, or keep it next to the "weights" folder and
import it:

    from dnam_network_clock import predict_age
    ages = predict_age(betas)        # betas: pandas DataFrame, samples x CpGs

The clock is a K=60 principal component ridge regression on per module medians
of methylation beta values. Modules were discovered by Infomap clustering of a
CpG to CpG Spearman correlation network (threshold tau=0.70) over 7,328 CpGs
pre filtered by absolute Spearman rho with age >= 0.35, on a training set of
1,917 whole blood methylation samples spanning ages 18 to 94 across 12 studies.
See README.md for the full methodology, the trained model metrics, and the
paper reference.

Requires: numpy >= 1.20, pandas >= 1.3 (Python >= 3.9).
Licence: MIT.
"""
from __future__ import annotations

import json
import warnings
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd


# Weight loading, cached once per process.
_WEIGHTS_CACHE: dict = {}


def _locate_weights_dir() -> Path:
    """Find the "weights" folder shipped alongside this file.

    Searched in order: the current working directory, the folder this file
    lives in, and the parent of that folder (the repository layout, where the
    code sits in python/ and the weights sit in ../weights). If none is found,
    pass weights_dir= explicitly to predict_age().
    """
    here = Path(__file__).resolve().parent
    candidates = [
        Path.cwd() / "weights",
        here / "weights",
        here.parent / "weights",
    ]
    for c in candidates:
        if c.is_dir():
            return c
    raise FileNotFoundError(
        "Could not locate the 'weights' directory. Pass weights_dir= "
        "explicitly to predict_age(), e.g. predict_age(betas, weights_dir='/path/to/weights')."
    )


def _load_weights(weights_dir: Optional[Union[str, Path]] = None) -> dict:
    """Load and cache the K=60 clock weights from the CSV files in weights_dir."""
    key = str(weights_dir) if weights_dir is not None else "<auto>"
    if key in _WEIGHTS_CACHE:
        return _WEIGHTS_CACHE[key]

    d = Path(weights_dir) if weights_dir is not None else _locate_weights_dir()

    defs   = pd.read_csv(d / "module_definitions.csv")
    params = pd.read_csv(d / "module_params.csv")
    comps  = pd.read_csv(d / "pca_components_k60.csv")
    ridge  = pd.read_csv(d / "ridge_k60.csv")
    with open(d / "metadata.json") as f:
        meta = json.load(f)

    # Module column order is given by module_params.csv.
    module_order = params["Module"].to_numpy()
    assert len(module_order) == meta["n_modules"]

    # Per module lookup: Module -> (list of CpGs, strategy).
    grouped = defs.groupby("Module", sort=False)
    module_to_cpgs = {int(m): g["CpG"].tolist() for m, g in grouped}
    module_to_strategy = {int(m): g["Strategy"].iloc[0] for m, g in grouped}

    # PCA components: columns are module IDs in module_order, rows are PCs.
    pca_components = comps.drop(columns=["PC"]).to_numpy()
    assert pca_components.shape[0] == meta["K"]
    assert pca_components.shape[1] == meta["n_modules"]

    # Ridge: PC1..PCK rows, plus an intercept row.
    pc_rows = ridge[ridge["name"].str.startswith("PC")]
    pc_coefs = pc_rows["value"].to_numpy()
    intercept = float(ridge.loc[ridge["name"] == "intercept", "value"].iloc[0])
    assert pc_coefs.shape == (meta["K"],)

    weights = {
        "module_order":       module_order,
        "module_to_cpgs":     module_to_cpgs,
        "module_to_strategy": module_to_strategy,
        "scaler_mean":        params["scaler_mean"].to_numpy(),
        "scaler_scale":       params["scaler_scale"].to_numpy(),
        "pca_mean":           params["pca_mean"].to_numpy(),
        "pca_components":     pca_components,
        "pc_coefs":           pc_coefs,
        "intercept":          intercept,
        "metadata":           meta,
    }
    _WEIGHTS_CACHE[key] = weights
    return weights


# Per module aggregation: build samples x modules matrix from samples x CpGs.
def _build_module_matrix(betas: pd.DataFrame, weights: dict,
                         warn_missing_above: float = 0.05) -> "tuple[np.ndarray, dict]":
    """Aggregate the input beta matrix (samples x CpGs) into a samples x modules
    matrix using each module's stored strategy.

    Missing CpGs are tolerated. 'median' modules use the median of the
    available CpGs per sample, falling back to the training scaler_mean when a
    module has no input CpGs at all. 'singleton' modules use the single CpG's
    value, or scaler_mean if that CpG is missing.
    """
    n_samples = betas.shape[0]
    module_order = weights["module_order"]
    n_modules = len(module_order)
    scaler_mean = weights["scaler_mean"]

    X = np.empty((n_samples, n_modules), dtype=np.float64)
    available_cpgs = set(betas.columns.astype(str))

    n_missing_modules = 0
    n_singleton_fallback = 0
    n_partial_median = 0
    n_cpgs_used = 0
    n_cpgs_needed = 0

    for j, mid in enumerate(module_order):
        cpgs = weights["module_to_cpgs"][int(mid)]
        strategy = weights["module_to_strategy"][int(mid)]
        n_cpgs_needed += len(cpgs)
        present = [c for c in cpgs if c in available_cpgs]
        n_cpgs_used += len(present)

        if strategy == "singleton":
            if present:
                X[:, j] = betas[present[0]].to_numpy()
            else:
                X[:, j] = scaler_mean[j]
                n_singleton_fallback += 1
        else:  # 'median'
            if not present:
                X[:, j] = scaler_mean[j]
                n_missing_modules += 1
            else:
                if len(present) < len(cpgs):
                    n_partial_median += 1
                # Median across the available CpGs for each sample.
                X[:, j] = np.nanmedian(betas[present].to_numpy(), axis=1)
                # If every CpG of a module was NaN for a sample, fall back.
                nan_mask = np.isnan(X[:, j])
                if nan_mask.any():
                    X[nan_mask, j] = scaler_mean[j]

    coverage = n_cpgs_used / n_cpgs_needed if n_cpgs_needed else 0.0
    diagnostics = {
        "n_modules": n_modules,
        "n_cpgs_used": n_cpgs_used,
        "n_cpgs_needed": n_cpgs_needed,
        "cpg_coverage": coverage,
        "n_modules_filled_with_training_mean": n_missing_modules + n_singleton_fallback,
        "n_partial_median_modules": n_partial_median,
    }
    missing_frac = 1.0 - coverage
    if missing_frac > warn_missing_above:
        warnings.warn(
            f"DNAm Network Clock: {missing_frac*100:.1f}% of the clock's CpGs are "
            f"missing from the input ({n_cpgs_used:,}/{n_cpgs_needed:,} present). "
            f"{diagnostics['n_modules_filled_with_training_mean']:,} modules were "
            f"filled with the training mean. Predictions may be biased toward "
            f"the training set mean age.",
            UserWarning, stacklevel=2,
        )
    return X, diagnostics


# Public API.
def predict_age(
    betas: pd.DataFrame,
    weights_dir: Optional[Union[str, Path]] = None,
    return_diagnostics: bool = False,
) -> Union[np.ndarray, "tuple[np.ndarray, dict]"]:
    """Predict chronological age from a methylation beta value matrix.

    Parameters
    ----------
    betas
        Methylation beta values as a pandas DataFrame with samples as rows and
        CpG IDs (Illumina cg identifiers) as columns. Missing CpGs and within
        column NaN are tolerated, see Notes.
    weights_dir
        Optional path to a directory of weight CSVs. By default the "weights"
        folder shipped with this file is located automatically.
    return_diagnostics
        If True, also return a dict summarising CpG coverage and which modules
        were imputed from the training mean.

    Returns
    -------
    ndarray of shape (n_samples,)
        Predicted ages in years, in the same order as betas.index. If
        return_diagnostics=True, a 2 tuple (ages, diagnostics).

    Notes
    -----
    The clock has 4,595 modules covering 7,328 CpGs. Each module is aggregated
    to a single value per sample by taking the median of the available CpGs
    (multi CpG modules) or using the single CpG (singleton modules). A module
    with no input CpGs at all is filled with the training time mean for that
    module. A warning is emitted if more than 5% of the clock's CpGs are
    missing. The model is trained on whole blood methylation, performance on
    other tissues has not been characterised.

    Examples
    --------
    >>> import pandas as pd
    >>> from dnam_network_clock import predict_age
    >>> betas = pd.read_csv("my_methylation.csv", index_col=0)   # samples x CpGs
    >>> ages = predict_age(betas)
    >>> ages[:5]
    array([42.1, 67.8, 28.3, 55.4, 71.9])
    """
    if not isinstance(betas, pd.DataFrame):
        raise TypeError(
            "betas must be a pandas DataFrame with samples as rows and CpG "
            "IDs as columns. Got: " + type(betas).__name__
        )
    if betas.shape[0] == 0:
        raise ValueError("betas has 0 rows (samples).")

    weights = _load_weights(weights_dir)

    # 1. samples x modules matrix
    X, diagnostics = _build_module_matrix(betas, weights)

    # 2. standardise
    X = (X - weights["scaler_mean"]) / weights["scaler_scale"]

    # 3. centre for PCA, then project onto the top K PCs
    X = X - weights["pca_mean"]
    PC_scores = X @ weights["pca_components"].T          # (n_samples, 60)

    # 4. ridge
    ages = PC_scores @ weights["pc_coefs"] + weights["intercept"]

    if return_diagnostics:
        return ages, diagnostics
    return ages


def clock_info(weights_dir: Optional[Union[str, Path]] = None) -> dict:
    """Return the trained clock's metadata (K, module counts, training metrics)."""
    return dict(_load_weights(weights_dir)["metadata"])


__all__ = ["predict_age", "clock_info"]
