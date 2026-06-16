# DNAm Network Clock, age prediction from methylation beta values (R)
#
# This is a single, self contained file that uses only base R. There is
# nothing to install. Source it into your own script:
#
#     source("path/to/dnam_network_clock.R")
#     ages <- predict_age(betas)
#
# where `betas` is a data.frame or matrix with samples as ROWS and CpG IDs
# (Illumina cg identifiers) as COLUMN names. Missing CpGs are tolerated.
#
# The clock is a K=60 principal component ridge regression on per module
# medians of methylation beta values. See README.md for the full methodology,
# the trained model metrics, and the paper reference.
#
# Author: Anton Carcedo
# Licence: MIT

# Record the directory this file lives in, at source time, so the weights
# folder can be found automatically later. Scans the call frames for the path
# that source() records, with a fallback for Rscript.
.dnam_this_script_dir <- function() {
  for (i in seq_len(sys.nframe())) {
    of <- sys.frame(i)$ofile
    if (!is.null(of)) return(normalizePath(dirname(of)))
  }
  args <- commandArgs(trailingOnly = FALSE)
  fa <- grep("^--file=", args, value = TRUE)
  if (length(fa) > 0) {
    return(normalizePath(dirname(sub("^--file=", "", fa[1]))))
  }
  NA_character_
}
.DNAM_SCRIPT_DIR <- .dnam_this_script_dir()


# Locate the "weights" folder.
#
# Tries, in order: getOption("dnam_clock.weights_dir"), a "weights" folder next
# to this file, a "weights" folder one level up (the repository layout, where
# this file sits in R/ and the weights sit in ../weights), and "weights" in the
# current working directory.
.dnam_locate_weights_dir <- function() {
  opt <- getOption("dnam_clock.weights_dir", NULL)
  if (!is.null(opt) && dir.exists(opt)) return(normalizePath(opt))

  candidates <- character(0)
  if (!is.na(.DNAM_SCRIPT_DIR)) {
    candidates <- c(candidates,
                    file.path(.DNAM_SCRIPT_DIR, "weights"),
                    file.path(dirname(.DNAM_SCRIPT_DIR), "weights"))
  }
  candidates <- c(candidates, file.path(getwd(), "weights"))

  for (p in candidates) if (dir.exists(p)) return(normalizePath(p))

  stop("DNAm Network Clock: could not locate the 'weights' directory. ",
       "Set options(dnam_clock.weights_dir = '/path/to/weights') ",
       "or pass weights_dir= to predict_age().")
}

# Lightweight string hash for cache keys (avoids a 'digest' package dependency).
.dnam_digest_str <- function(s) {
  chars <- strsplit(s, "")[[1]]
  paste0("h", as.character(sum(utf8ToInt(s) * seq_along(chars))))
}

# Load and cache the K=60 clock weights from the CSV files.
.dnam_load_weights <- function(weights_dir = NULL) {
  if (is.null(weights_dir)) weights_dir <- .dnam_locate_weights_dir()
  cache_key <- paste0(".dnam_weights_cache_", .dnam_digest_str(weights_dir))
  if (exists(cache_key, envir = .GlobalEnv)) {
    return(get(cache_key, envir = .GlobalEnv))
  }

  defs   <- read.csv(file.path(weights_dir, "module_definitions.csv"),
                     stringsAsFactors = FALSE)
  params <- read.csv(file.path(weights_dir, "module_params.csv"),
                     stringsAsFactors = FALSE)
  comps  <- read.csv(file.path(weights_dir, "pca_components_k60.csv"),
                     stringsAsFactors = FALSE, check.names = FALSE)
  ridge  <- read.csv(file.path(weights_dir, "ridge_k60.csv"),
                     stringsAsFactors = FALSE)

  module_order <- params$Module                   # canonical column order

  # Split CpGs by module, then reorder to module_order.
  defs_split <- split(defs$CpG, defs$Module)
  module_to_cpgs <- defs_split[as.character(module_order)]

  # Strategy, one per module (take the first row).
  strat_per_module <- tapply(defs$Strategy, defs$Module, function(x) x[1])
  module_to_strategy <- as.character(strat_per_module[as.character(module_order)])

  # PCA components: 60 rows x 4595 modules, drop the 'PC' label column.
  pca_components <- as.matrix(comps[, -1, drop = FALSE])

  # Ridge coefficients plus intercept.
  pc_mask    <- grepl("^PC", ridge$name)
  pc_coefs   <- ridge$value[pc_mask]
  intercept  <- ridge$value[ridge$name == "intercept"]

  stopifnot(length(module_order)     == 4595,
            length(module_to_cpgs)   == 4595,
            length(pc_coefs)         == 60,
            nrow(pca_components)      == 60,
            ncol(pca_components)      == 4595,
            length(intercept)         == 1)

  weights <- list(
    module_order       = module_order,
    module_to_cpgs     = module_to_cpgs,
    module_to_strategy = module_to_strategy,
    scaler_mean        = params$scaler_mean,
    scaler_scale       = params$scaler_scale,
    pca_mean           = params$pca_mean,
    pca_components     = pca_components,
    pc_coefs           = pc_coefs,
    intercept          = intercept,
    n_modules          = length(module_order),
    K                  = length(pc_coefs)
  )
  assign(cache_key, weights, envir = .GlobalEnv)
  weights
}


# Build the samples x modules matrix from samples x CpGs.
# Missing CpGs and NA values are tolerated.
.dnam_build_module_matrix <- function(betas, weights, warn_missing_above = 0.05) {
  if (is.data.frame(betas)) betas <- as.matrix(betas)
  if (!is.matrix(betas))
    stop("betas must be a matrix or data.frame (samples x CpGs).")
  if (is.null(colnames(betas)))
    stop("betas must have CpG IDs as column names.")

  n_samples <- nrow(betas)
  n_modules <- weights$n_modules
  X <- matrix(0.0, nrow = n_samples, ncol = n_modules)
  available <- colnames(betas)

  n_missing_modules     <- 0L
  n_singleton_fallback  <- 0L
  n_partial_median      <- 0L
  n_cpgs_used           <- 0L
  n_cpgs_needed         <- 0L

  for (j in seq_len(n_modules)) {
    cpgs     <- weights$module_to_cpgs[[j]]
    strategy <- weights$module_to_strategy[j]
    n_cpgs_needed <- n_cpgs_needed + length(cpgs)
    present  <- intersect(cpgs, available)
    n_cpgs_used <- n_cpgs_used + length(present)

    if (strategy == "singleton") {
      if (length(present) == 1L) {
        X[, j] <- betas[, present]
      } else {
        X[, j] <- weights$scaler_mean[j]
        n_singleton_fallback <- n_singleton_fallback + 1L
      }
    } else {  # 'median'
      if (length(present) == 0L) {
        X[, j] <- weights$scaler_mean[j]
        n_missing_modules <- n_missing_modules + 1L
      } else {
        if (length(present) < length(cpgs)) n_partial_median <- n_partial_median + 1L
        if (length(present) == 1L) {
          X[, j] <- betas[, present]
        } else {
          X[, j] <- apply(betas[, present, drop = FALSE], 1, median, na.rm = TRUE)
        }
        nan_rows <- is.na(X[, j])
        if (any(nan_rows)) X[nan_rows, j] <- weights$scaler_mean[j]
      }
    }
  }

  coverage <- if (n_cpgs_needed > 0) n_cpgs_used / n_cpgs_needed else 0
  missing_frac <- 1 - coverage
  if (missing_frac > warn_missing_above) {
    warning(sprintf(
      paste("DNAm Network Clock: %.1f%% of the clock's CpGs are missing from",
            "the input (%d/%d present). %d modules were filled with the",
            "training mean. Predictions may be biased toward the training",
            "set mean age."),
      missing_frac * 100, n_cpgs_used, n_cpgs_needed,
      n_missing_modules + n_singleton_fallback
    ), call. = FALSE)
  }

  list(
    X = X,
    diagnostics = list(
      n_modules = n_modules,
      n_cpgs_used = n_cpgs_used,
      n_cpgs_needed = n_cpgs_needed,
      cpg_coverage = coverage,
      n_modules_filled_with_training_mean = n_missing_modules + n_singleton_fallback,
      n_partial_median_modules = n_partial_median
    )
  )
}


# Predict chronological age from a methylation beta matrix.
#
# betas              A matrix or data.frame with samples as rows and CpG IDs
#                    (Illumina cg identifiers) as column names.
# weights_dir        Optional path to a directory of weight CSVs. If NULL, the
#                    "weights" folder is located automatically.
# return_diagnostics If TRUE, return a list(ages, diagnostics).
#
# Returns a numeric vector of predicted ages (years), one per sample in the row
# order of `betas`. If return_diagnostics=TRUE, a list with elements `ages` and
# `diagnostics`.
predict_age <- function(betas, weights_dir = NULL, return_diagnostics = FALSE) {
  if (is.null(rownames(betas)) && is.data.frame(betas))
    rownames(betas) <- seq_len(nrow(betas))

  weights <- .dnam_load_weights(weights_dir)

  built <- .dnam_build_module_matrix(betas, weights)
  X <- built$X

  # standardise, then centre, then project, then ridge.
  X <- sweep(X, 2, weights$scaler_mean,  "-")
  X <- sweep(X, 2, weights$scaler_scale, "/")
  X <- sweep(X, 2, weights$pca_mean,     "-")
  PC_scores <- X %*% t(weights$pca_components)               # n x 60
  ages <- as.numeric(PC_scores %*% weights$pc_coefs + weights$intercept)

  if (return_diagnostics) list(ages = ages, diagnostics = built$diagnostics) else ages
}

# Return the trained clock's metadata.
clock_info <- function(weights_dir = NULL) {
  if (is.null(weights_dir)) weights_dir <- .dnam_locate_weights_dir()
  meta_path <- file.path(weights_dir, "metadata.json")
  if (!file.exists(meta_path))
    stop("metadata.json not found in ", weights_dir)
  raw <- paste(readLines(meta_path, warn = FALSE), collapse = "\n")
  if (requireNamespace("jsonlite", quietly = TRUE)) {
    jsonlite::fromJSON(raw)
  } else {
    message("Install 'jsonlite' for a structured metadata object, ",
            "returning the raw JSON string.")
    raw
  }
}
