# Run from the repository root: Rscript example_ABC.R
# These data are generated entirely for demonstration; they are not experimental data.
source("R/correspondence.R")
if (!requireNamespace("Matrix", quietly = TRUE) || !requireNamespace("glmnet", quietly = TRUE))
  stop('Install dependencies first: install.packages(c("Matrix", "glmnet"))')

make_example_dataset <- function(name, cells_per_population, seed) {
  set.seed(seed)
  n_genes <- 80L
  labels <- rep(paste0("population_", 1:3), each = cells_per_population)
  cells <- paste0(name, "_cell_", seq_along(labels))
  rate <- matrix(0.5, n_genes, length(labels))
  for (k in 1:3) rate[((k - 1L) * 15L + 1L):(k * 15L), labels == paste0("population_", k)] <- 5
  # Independent count noise and varying library sizes in each synthetic dataset.
  rate <- sweep(rate, 2, runif(length(labels), 0.6, 1.4), "*")
  counts <- matrix(rpois(length(rate), lambda = as.vector(rate)), nrow = n_genes,
                   dimnames = list(paste0("gene_", seq_len(n_genes)), cells))
  expression <- log1p(sweep(counts, 2, colSums(counts), "/") * 10000)
  expression <- expression[sample(seq_len(n_genes)), , drop = FALSE]
  list(expression = Matrix::Matrix(expression, sparse = TRUE),
       labels = stats::setNames(labels, cells), population_order = paste0("population_", 1:3))
}

datasets <- list(A = make_example_dataset("A", 36L, 101L),
                 B = make_example_dataset("B", 40L, 102L),
                 C = make_example_dataset("C", 32L, 103L))
comparisons <- data.frame(reference = c("B", "B", "A"), query = c("A", "C", "C"))
args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else paste0("example_output_", format(Sys.time(), "%Y%m%d_%H%M%S"))
if (file.exists(output_dir)) stop("Use a new output folder: ", output_dir)
if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create: ", output_dir)
saveRDS(datasets, file.path(output_dir, "SYNTHETIC_example_datasets_ABC.rds"))
writeLines("SYNTHETIC DEMONSTRATION DATA AND OUTPUTS. No experimental data or results.",
           file.path(output_dir, "SYNTHETIC_DATA_NOTICE.txt"))
utils::write.csv(comparisons, file.path(output_dir, "comparisons.csv"), row.names = FALSE)
for (i in seq_len(nrow(comparisons))) {
  ref <- comparisons$reference[i]
  qry <- comparisons$query[i]
  key <- paste(ref, "to", qry, sep = "_")
  result <- fit_correspondence(datasets[[ref]], datasets[[qry]], ref, qry)
  stopifnot(identical(dim(result$scores), c(3L, 3L)), all(is.finite(result$scores)),
            all(result$scores >= 0 & result$scores <= 1))
  write_correspondence(result, file.path(output_dir, key))
  plot_correspondence(result, file.path(output_dir, paste0(key, "_SYNTHETIC.pdf")),
                       paste("SYNTHETIC:", ref, "reference ->", qry, "query"))
}
capture.output(utils::sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
message("Synthetic example completed: ", normalizePath(output_dir))
