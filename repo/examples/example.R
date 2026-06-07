# Minimal usage example for the DNAm Network Clock (R).
#
# Run from the repository root:
#
#     Rscript examples/example.R

source("R/dnam_network_clock.R")

# Load reference input shipped with the repo
betas <- read.csv("examples/example_input.csv",
                  row.names = 1, check.names = FALSE)
cat(sprintf("Loaded reference input: %d samples × %d CpGs\n",
            nrow(betas), ncol(betas)))

# Predict
ages <- predict_age(betas)
cat(sprintf("\nPredictions:\n"))
cat(sprintf("  range  : %.2f – %.2f years\n", min(ages), max(ages)))
cat(sprintf("  median : %.2f years\n", median(ages)))

# Verify against stored expected output (Python-generated)
expected <- read.csv("examples/expected_output.csv")
max_diff <- max(abs(ages - expected$predicted_age))
cat(sprintf("\n  max |R prediction − expected| = %.2e   %s\n",
            max_diff, if (max_diff < 1e-6) "✓ PASS" else "✗ FAIL"))
