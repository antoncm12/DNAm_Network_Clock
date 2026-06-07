"""Minimal usage example for the DNAm Network Clock.

Run from the repository root:

    python examples/example.py
"""
from pathlib import Path
import numpy as np
import pandas as pd

from dnam_network_clock import predict_age, clock_info


HERE = Path(__file__).resolve().parent

# Inspect the bundled clock
info = clock_info()
print("Clock metadata:")
print(f"  K = {info['K']} PCs over {info['n_modules']:,} modules ({info['n_cpgs']:,} CpGs)")
print(f"  Trained R²(test) = {info['training_metrics']['r2_test']:.3f}, "
      f"MAE(test) = {info['training_metrics']['mae_test']:.2f} years")

# Apply to the reference example shipped with the repo
betas = pd.read_csv(HERE / "example_input.csv", index_col=0)
print(f"\nLoaded reference input: {betas.shape[0]} samples × {betas.shape[1]:,} CpGs")

ages, diag = predict_age(betas, return_diagnostics=True)
print(f"\nPredictions:")
print(f"  range  : {ages.min():.2f} – {ages.max():.2f} years")
print(f"  median : {np.median(ages):.2f} years")
print(f"  CpG coverage : {diag['cpg_coverage']*100:.1f}%")

# Check against the stored expected output
expected = pd.read_csv(HERE / "expected_output.csv")
max_diff = float(np.abs(ages - expected["predicted_age"].to_numpy()).max())
print(f"\n  max |predicted − expected| = {max_diff:.2e}  "
      f"{'✓ PASS' if max_diff < 1e-6 else '✗ FAIL'}")
