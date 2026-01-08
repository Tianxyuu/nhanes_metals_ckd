# Inspect model summaries and key metrics

models_path <- "data/models.rds"
if (!file.exists(models_path)) {
  stop("Missing data/models.rds. Run R/04_models.R first.")
}

models <- readRDS(models_path)

cat("=== fit_lm_full summary ===\n")
print(summary(models$fit_lm_full))

cat("\n=== fit_svy_full summary ===\n")
print(summary(models$fit_svy_full))

cat("\n=== fit_lm_full R-squared ===\n")
cat("R2:", summary(models$fit_lm_full)$r.squared, "\n")
cat("Adj R2:", summary(models$fit_lm_full)$adj.r.squared, "\n")

spline_path <- "data/spline_tests.rds"
if (file.exists(spline_path)) {
  cat("\n=== spline_tests ===\n")
  print(readRDS(spline_path))
}

invisible(TRUE)
