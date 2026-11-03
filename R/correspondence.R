# Directed transcriptional correspondence between annotated datasets.
# This file defines functions; sourcing it does not fit models or write files.

check_dataset <- function(x, name = "dataset") {
  if (!is.list(x) || !all(c("expression", "labels", "population_order") %in% names(x)))
    stop(name, ": provide expression, labels and population_order.")
  m <- x$expression
  if (!(is.matrix(m) || inherits(m, "Matrix")) || !is.numeric(m) && !inherits(m, "Matrix"))
    stop(name, ": expression must be a numeric genes-by-cells matrix.")
  if (nrow(m) < 2L || ncol(m) < 2L || is.null(rownames(m)) || is.null(colnames(m)) ||
      anyNA(rownames(m)) || anyNA(colnames(m)) || any(!nzchar(rownames(m))) ||
      any(!nzchar(colnames(m))) || anyDuplicated(rownames(m)) || anyDuplicated(colnames(m)))
    stop(name, ": use unique, nonempty gene and cell identifiers.")
  if (inherits(m, "sparseMatrix") && !"x" %in% methods::slotNames(m))
    stop(name, ": expression must contain numeric values, not a pattern-only sparse matrix.")
  values <- if (inherits(m, "sparseMatrix")) m@x else as.vector(m)
  if (!is.numeric(values) || any(!is.finite(values)))
    stop(name, ": expression values must be finite numbers.")
  if (!is.character(x$labels) && !is.factor(x$labels))
    stop(name, ": labels must be a named character vector or factor.")
  if (!identical(names(x$labels), colnames(m)) || anyNA(x$labels) ||
      any(!nzchar(as.character(x$labels))))
    stop(name, ": label names and order must exactly match expression cell columns.")
  if (!is.character(x$population_order) || anyNA(x$population_order) ||
      anyDuplicated(x$population_order) || length(x$population_order) < 2L ||
      !setequal(x$population_order, as.character(x$labels)))
    stop(name, ": population_order must list every observed population exactly once.")
  invisible(x)
}

align_dataset_pair <- function(reference, query) {
  check_dataset(reference, "reference")
  check_dataset(query, "query")
  genes <- intersect(rownames(reference$expression), rownames(query$expression))
  if (length(genes) < 2L) stop("Fewer than two shared gene identifiers.")
  train_x <- Matrix::t(reference$expression[genes, , drop = FALSE])
  query_x <- Matrix::t(query$expression[genes, , drop = FALSE])
  stopifnot(identical(colnames(train_x), colnames(query_x)))
  list(genes = genes, train_x = train_x, query_x = query_x)
}

aggregate_correspondence <- function(eta, query_labels, query_order) {
  if (!is.matrix(eta) || !is.numeric(eta) || any(!is.finite(eta)) ||
      is.null(colnames(eta)) || length(query_labels) != nrow(eta) ||
      anyNA(query_labels) || anyDuplicated(query_order) ||
      !setequal(query_order, as.character(query_labels)))
    stop("Invalid linear predictors, labels or query population order.")
  scores <- matrix(NA_real_, ncol(eta), length(query_order),
                   dimnames = list(colnames(eta), query_order))
  for (q in query_order)
    scores[, q] <- stats::plogis(colMeans(eta[query_labels == q, , drop = FALSE]))
  scores
}

fit_correspondence <- function(reference, query, reference_name = "reference",
                               query_name = "query", seed = 2026L) {
  needed <- c("Matrix", "glmnet")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "))
  if (length(seed) != 1L || is.na(seed) || !is.finite(seed) || seed != as.integer(seed))
    stop("seed must be one integer.")
  aligned <- align_dataset_pair(reference, query)
  train_x <- aligned$train_x
  query_x <- aligned$query_x
  if (nrow(train_x) < 12L) stop("At least 12 reference cells are needed for 12-fold CV.")
  reference_order <- reference$population_order
  query_order <- query$population_order

  # Set a fixed RNG kind and reset once per directed pair. Restore the caller's RNG.
  old_kind <- RNGkind()
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  on.exit({
    do.call(RNGkind, as.list(old_kind))
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  set.seed(seed)
  rng_start <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  fit_cv <- if ("control" %in% names(formals(glmnet::glmnet)))
    function(...) glmnet::cv.glmnet(..., control = list(thresh = 1e-7, maxit = 100000L)) else
    function(...) glmnet::cv.glmnet(..., thresh = 1e-7, maxit = 100000L)
  eta <- matrix(NA_real_, nrow(query_x), length(reference_order),
                dimnames = list(rownames(query_x), reference_order))
  folds <- coefficients <- settings <- vector("list", length(reference_order))
  names(folds) <- names(coefficients) <- names(settings) <- reference_order
  lambda_ratio <- if (nrow(train_x) < ncol(train_x)) 0.01 else 0.0001
  for (g in reference_order) {
    message("Fitting ", reference_name, " -> ", query_name, ": ", g)
    y <- factor(reference$labels == g, levels = c(FALSE, TRUE))
    if (any(table(y) < 2L)) stop(g, ": at least two cells in each response class are required.")
    prevalence <- mean(y == "TRUE")
    offset <- stats::qlogis(prevalence)
    foldid <- sample(rep(seq_len(12L), length.out = nrow(train_x)))
    fit <- tryCatch(fit_cv(
      x = train_x, y = y, family = "binomial", offset = rep(offset, nrow(train_x)),
      intercept = FALSE, alpha = 0.99, standardize = TRUE,
      nfolds = 12L, foldid = foldid, type.measure = "class",
      nlambda = 100L, lambda.min.ratio = lambda_ratio, keep = FALSE, parallel = FALSE
    ), error = function(e) stop(reference_name, " -> ", query_name, " / ", g,
                               ": ", conditionMessage(e), call. = FALSE))
    prediction <- as.numeric(stats::predict(fit, newx = query_x, type = "link",
                         newoffset = rep(0, nrow(query_x)), s = "lambda.min"))
    if (any(!is.finite(prediction))) stop(g, ": non-finite predictions.")
    eta[, g] <- prediction
    folds[[g]] <- foldid
    beta <- as.matrix(stats::coef(fit, s = "lambda.min"))[, 1]
    beta <- beta[names(beta) != "(Intercept)"]
    nonzero <- beta[beta != 0]
    nonzero <- nonzero[order(-abs(nonzero))]
    coefficients[[g]] <- data.frame(
      reference_population = rep(g, length(nonzero)), gene = names(nonzero),
      coefficient = as.numeric(nonzero), rank_by_abs = seq_along(nonzero))
    settings[[g]] <- data.frame(
      reference_dataset = reference_name, query_dataset = query_name,
      reference_population = g, reference_cells = nrow(train_x),
      positive_cells = sum(y == "TRUE"), prevalence = prevalence,
      training_offset = offset, prediction_offset = 0, alpha = 0.99,
      intercept = FALSE, standardize = TRUE, nfolds = 12L,
      CV_measure = "misclassification error", lambda_min = fit$lambda.min,
      nlambda = 100L, lambda_min_ratio = lambda_ratio, thresh = 1e-7,
      maxit = 100000L, shared_genes = length(aligned$genes),
      nonzero_genes = length(nonzero), seed = seed)
    rm(fit)
  }
  scores <- aggregate_correspondence(eta, query$labels, query_order)
  list(reference_name = reference_name, query_name = query_name,
       scores = scores, per_cell_linear_predictors = eta, shared_genes = aligned$genes,
       fold_ids = data.frame(reference_cell_id = rownames(train_x),
                             as.data.frame(folds, check.names = FALSE), check.names = FALSE),
       nonzero_coefficients = do.call(rbind, coefficients),
       settings = do.call(rbind, settings),
       reference_labels = reference$labels, query_labels = query$labels,
       reference_population_order = reference_order, query_population_order = query_order,
       RNG_kind = RNGkind(), RNG_state_at_pair_start = rng_start,
       session_info = utils::sessionInfo())
}

write_correspondence <- function(result, output_dir) {
  if (file.exists(output_dir)) stop("Use a new output folder: ", output_dir)
  if (!dir.create(output_dir, recursive = TRUE)) stop("Cannot create: ", output_dir)
  write_table <- function(x, file) utils::write.csv(x, file.path(output_dir, file), row.names = FALSE)
  write_table(data.frame(reference_population = rownames(result$scores), result$scores,
                         check.names = FALSE), "scores.csv")
  write_table(result$fold_ids, "actual_CV_fold_ids.csv")
  write_table(result$settings, "model_settings.csv")
  write_table(result$nonzero_coefficients, "nonzero_coefficients.csv")
  write_table(data.frame(gene = result$shared_genes), "shared_genes.csv")
  saveRDS(result, file.path(output_dir, "analysis.rds"))
  capture.output(result$session_info, file = file.path(output_dir, "sessionInfo.txt"))
  invisible(output_dir)
}

plot_correspondence <- function(result, file, title = NULL) {
  if (file.exists(file)) stop("Use a new PDF path: ", file)
  if (is.null(title)) title <- paste(result$reference_name, "reference ->",
                                    result$query_name, "query")
  scores <- result$scores
  if (!is.matrix(scores) || any(!is.finite(scores)) || any(scores < 0 | scores > 1))
    stop("Invalid score matrix.")
  grDevices::pdf(file, width = 7, height = 6, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE)
  graphics::layout(matrix(c(1, 2), nrow = 1), widths = c(5, 1))
  colours <- grDevices::colorRampPalette(c("#FFFFFF", "#FDAE61", "#D73027"))(256L)
  graphics::par(mar = c(7, 8, 4, 1))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, ncol(scores)), ylim = c(0, nrow(scores)),
                        xaxs = "i", yaxs = "i")
  for (i in seq_len(nrow(scores))) for (j in seq_len(ncol(scores))) {
    index <- pmin(256L, 1L + floor(scores[i, j] * 255))
    graphics::rect(j - 1, nrow(scores) - i, j, nrow(scores) - i + 1,
                   col = colours[index], border = "#EEEEEE")
  }
  graphics::axis(1, at = seq_len(ncol(scores)) - 0.5, labels = colnames(scores),
                 las = 2, tick = FALSE)
  graphics::axis(2, at = nrow(scores) - seq_len(nrow(scores)) + 0.5,
                 labels = rownames(scores), las = 2, tick = FALSE)
  graphics::title(main = title, cex.main = 0.9)
  graphics::mtext("Query populations", side = 1, line = 5)
  graphics::mtext("Reference populations", side = 2, line = 6)
  graphics::par(mar = c(7, 0, 4, 4))
  graphics::plot.new()
  graphics::plot.window(xlim = c(0, 1), ylim = c(0, 1), xaxs = "i", yaxs = "i")
  for (k in seq_len(256L))
    graphics::rect(0.1, (k - 1) / 256, 0.5, k / 256, col = colours[k], border = NA)
  graphics::axis(4, at = c(0, 0.5, 1), labels = c("0", "0.5", "1"), las = 1)
  graphics::mtext("Correspondence score", side = 4, line = 2.5)
  invisible(file)
}
