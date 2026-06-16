"""Example usage of the DNAm Network Clock (Python).

Run it from the repository root:

    python python/example.py

The only requirement is that the "weights" folder is reachable. When you run
from the repo root it is found automatically. If you copy dnam_network_clock.py
into your own project, either keep a "weights" folder in your working directory
or pass weights_dir= to predict_age().
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

# Make sure Python can import the single file module sitting next to this script.
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from dnam_network_clock import predict_age, clock_info

REPO = HERE.parent
WEIGHTS = REPO / "weights"
EXAMPLES = REPO / "examples"

# Look at what the clock is.
info = clock_info(weights_dir=WEIGHTS)
print("Clock metadata:")
print(f"  K = {info['K']} PCs over {info['n_modules']:,} modules ({info['n_cpgs']:,} CpGs)")
print(f"  Trained R2(test) = {info['training_metrics']['r2_test']:.3f}, "
      f"MAE(test) = {info['training_metrics']['mae_test']:.2f} years")

# ---- The two lines that matter for your own data ----
# betas must be samples (rows) x CpG IDs (columns).
betas = pd.read_csv(EXAMPLES / "example_input.csv", index_col=0)
ages = predict_age(betas, weights_dir=WEIGHTS)
# -----------------------------------------------------

print(f"\nLoaded reference input: {betas.shape[0]} samples x {betas.shape[1]:,} CpGs")
print("Predictions:")
print(f"  range  : {ages.min():.2f} to {ages.max():.2f} years")
print(f"  median : {np.median(ages):.2f} years")

# Optional: coverage diagnostics.
_, diag = predict_age(betas, weights_dir=WEIGHTS, return_diagnostics=True)
print(f"  CpG coverage : {diag['cpg_coverage']*100:.1f}%")

# Check against the stored reference predictions.
expected = pd.read_csv(EXAMPLES / "expected_output.csv")
max_diff = float(np.abs(ages - expected["predicted_age"].to_numpy()).max())
print(f"\n  max |predicted - expected| = {max_diff:.2e}  "
      f"{'PASS' if max_diff < 1e-6 else 'FAIL'}")
