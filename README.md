# Cross-dataset transcriptional correspondence

A small R workflow for asking whether independently annotated populations in one dataset are recognised by expression profiles learned from another dataset. The repository includes a completely synthetic A/B/C example; no experimental data or results are included.

## The logic

1. **Choose a direction.** In A → B, A supplies the training cells and their population labels; B supplies the cells to examine. “Reference” describes the training role. It does not mean an experimental control.
2. **Align expression features.** Retain gene identifiers present in both datasets, and put them in the same order. Each pair has its own shared feature set.
3. **Learn each reference population separately.** Fit one binary classifier for population A1 versus all other A cells, then repeat for A2, A3 and so on. A population's prevalence is included as a fixed training offset.
4. **Apply the models to query cells.** Each B cell receives a linear predictor from each A-population model. The reference-prevalence offset is set to zero at prediction.
5. **Summarise by the existing query labels.** Average the linear predictors within each B population, then apply the inverse logit. Rows are reference populations and columns are query populations.

A high score means that a query population has expression features that favour a particular reference population in its fitted model. Scores from separate reference models do not have to sum to one. They are not classification accuracies or calibrated cell-type assignment probabilities.

Reversing A → B requires training new models on B. It is not obtained by transposing the A → B matrix: the training labels, cell composition, fitted coefficients and selected regularisation can change. Reciprocal comparisons can therefore provide different views of correspondence.

## Run the synthetic example

Install R and the two dependencies:

```r
install.packages(c("Matrix", "glmnet"))
```

Open a terminal in the repository folder and run:

```sh
Rscript tests/test_core.R
Rscript example_ABC.R
```

Or, with the repository folder as the working directory in RStudio:

```r
source("tests/test_core.R")
source("example_ABC.R")
```

The example generates three matrices called A, B and C, each with three generic populations. Its counts, gene names and population labels are simulated. It then fits **B → A, B → C and A → C**. These directions demonstrate the interface; edit the `comparisons` table to request different pairs. Example PDF titles and output files explicitly identify the data as synthetic.

A new `example_output_<timestamp>` folder contains the generated example inputs, comparison table, white–orange–red PDF heatmaps and one results folder per direction. An optional first command-line argument supplies a new output folder:

```sh
Rscript example_ABC.R my_new_output_folder
```

## Supply your own data

Source `R/correspondence.R` and construct one list per dataset:

```r
source("R/correspondence.R")

# expression_A: numeric matrix, genes in rows and cells in columns.
# labels_A: one population label per cell, in the same order as its columns.
A <- list(
  expression = expression_A,
  labels = setNames(as.character(labels_A), colnames(expression_A)),
  population_order = c("population_1", "population_2", "population_3")
)
B <- list(
  expression = expression_B,
  labels = setNames(as.character(labels_B), colnames(expression_B)),
  population_order = c("population_1", "population_2", "population_3")
)

result <- fit_correspondence(A, B, reference_name = "A", query_name = "B")
write_correspondence(result, "A_to_B_results")
plot_correspondence(result, "A_to_B.pdf")
```

Input requirements:

- Supply comparable normalised expression values, with the same gene identifier system and expression transformation across datasets. The functions do not normalise counts or correct batch effects.
- Matrices may be dense or numeric sparse `Matrix` objects. Gene names and cell names must be unique and nonempty within each dataset. Cell names can overlap between datasets because each dataset is processed separately.
- Named labels must exactly match the expression matrix's cell columns. No cells are assigned new labels by this workflow.
- `population_order` must list every observed population exactly once. It determines fitting/display order; keep it fixed when reproducing a run.
- At least two genes must be shared. The reference must contain at least 12 cells, and each reference-versus-rest response must contain at least two cells in both classes. Rare populations can still cause cross-validation fits to fail; errors identify the affected comparison and population.

## Score definition and model defaults

For reference population k, let πk be its fraction of reference cells. Training uses offset `log(πk / (1 - πk))` and no free intercept. For query cell i, ηik is the fitted linear predictor with prediction offset zero. For query population q:

```text
score[k, q] = inverse_logit(mean(eta[i, k] for query cells i in q))
```

The average is taken **before** the inverse logit. Averaging per-cell probabilities gives a different quantity.

| Setting | Default |
|---|---|
| Model | Binomial elastic-net logistic regression, one population versus rest |
| Features | All pairwise-shared input genes |
| Elastic-net α | 0.99 |
| Predictor standardisation | Enabled |
| Free intercept | Disabled |
| Training offset | Reference population's log prevalence odds |
| Prediction offset | Zero |
| Cross-validation | 12 unstratified folds with approximately equal cell counts |
| λ selection | `lambda.min`, using mean misclassification error |
| λ path | 100 requested values; minimum ratio 0.01 if cells < genes, otherwise 0.0001 |
| Convergence | Threshold 1e−7; maximum iterations 100,000 |
| RNG | Mersenne-Twister / Inversion / Rejection; seed 2026 reset once per directed pair |
| Execution | Sequential |

Folds are drawn separately for each classifier and saved. Query labels are used only to group predictions for display. Query expression values and labels do not select λ; the query's available gene identifiers determine the shared feature set. Cross-validation selects regularisation within the reference dataset; it does not measure generalisation across independent subjects. Expression differences can reflect technical and biological differences, so correspondence should be interpreted with the experimental design.

## Files and results

| File | Purpose |
|---|---|
| `R/correspondence.R` | Reusable fitting, aggregation, export and plotting functions |
| `example_ABC.R` | Generate synthetic A/B/C inputs and run three comparisons |
| `tests/test_core.R` | Check gene alignment, cell-label alignment and score aggregation |
| `.github/workflows/r-smoke.yml` | Run the checks and full synthetic example on GitHub Actions |
| `scores.csv` | Population-level matrix at full exported numerical precision |
| `actual_CV_fold_ids.csv` | Reference cell identifiers and actual folds for every classifier |
| `model_settings.csv` | Fitting settings and selected λ for every classifier |
| `nonzero_coefficients.csv` | Selected nonzero gene coefficients |
| `shared_genes.csv` | Aligned feature identifiers |
| `analysis.rds` | Scores, per-cell predictors, coefficients, folds, settings and provenance |
| `sessionInfo.txt` | R and loaded package versions |

Reproduction requires the same input values and cell order, population order, code, seed and software versions. Save an environment lockfile from a successful run if distributing a fixed numerical analysis. This repository does not contain a fabricated lockfile or precomputed example scores.

## Validation status

The source and configuration have been inspected, but this initial package has not been executed in R in the preparation environment. Run the supplied checks and example before relying on it. The included GitHub Actions workflow runs both after publication; no successful CI run is claimed here.

## Method references

- Friedman J, Hastie T, Tibshirani R. Regularization Paths for Generalized Linear Models via Coordinate Descent. *Journal of Statistical Software*. 2010;33(1):1–22. https://doi.org/10.18637/jss.v033.i01
- glmnet fitting and offsets: https://glmnet.stanford.edu/reference/glmnet.html
- glmnet cross-validation: https://glmnet.stanford.edu/reference/cv.glmnet.html
