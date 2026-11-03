# Run from the repository root: Rscript tests/test_core.R
source("R/correspondence.R")
if (!requireNamespace("Matrix", quietly = TRUE)) stop("Install Matrix first.")
mk <- function(m, groups) list(expression = m,
  labels = stats::setNames(groups, colnames(m)), population_order = unique(groups))
reference <- mk(matrix(1:12, nrow = 3,
  dimnames = list(c("g1", "g2", "g3"), paste0("r", 1:4))), c("x", "x", "y", "y"))
query <- mk(matrix(21:29, nrow = 3,
  dimnames = list(c("g3", "g1", "g4"), paste0("q", 1:3))), c("u", "v", "v"))
aligned <- align_dataset_pair(reference, query)
stopifnot(identical(aligned$genes, c("g1", "g3")),
  identical(as.numeric(aligned$query_x[1, ]), c(22, 21)),
  identical(rownames(aligned$train_x), colnames(reference$expression)),
  identical(rownames(aligned$query_x), colnames(query$expression)))

# The direction determines the training cells and feature order, not a transpose operation.
reverse <- align_dataset_pair(query, reference)
stopifnot(identical(reverse$genes, c("g3", "g1")),
  identical(rownames(reverse$train_x), colnames(query$expression)))

# This asymmetric example distinguishes inverse-logit-of-mean from mean probability.
eta <- matrix(c(-4, 0, 1, 2, 2, -1), nrow = 3,
              dimnames = list(paste0("q", 1:3), c("x", "y")))
scores <- aggregate_correspondence(eta, c("u", "u", "v"), c("v", "u"))
stopifnot(identical(rownames(scores), c("x", "y")),
  identical(colnames(scores), c("v", "u")),
  abs(scores["x", "u"] - stats::plogis(-2)) < 1e-12,
  abs(scores["x", "u"] - mean(stats::plogis(c(-4, 0)))) > 0.1)

# Misaligned cell annotations must fail rather than silently relabel cells.
bad <- reference
bad$labels <- bad$labels[4:1]
stopifnot(inherits(tryCatch(check_dataset(bad), error = identity), "error"))
message("Core alignment and score-definition checks passed.")
