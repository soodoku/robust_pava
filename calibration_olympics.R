library(dplyr)
library(tidyr)
library(isotone)
library(mgcv)

## Methods

robust_pava <- function(target, pred, threshold = .05) {
  # Sort everything by predictions
  ord <- order(pred)
  y <- target[ord]  # target values we want to match
  w <- rep(1, length(y))  # weights
  y_iso <- y  # isotonic regression result
  
  iteration <- 0
  while(TRUE) {
    iteration <- iteration + 1
    diffs <- diff(y_iso)
    violations <- which(diffs < 0)
    violations <- violations[abs(diffs[violations]) > threshold]
    
    if(length(violations) == 0 || iteration > 1000) break
    
    i <- violations[1]
    pool_value <- (w[i] * y_iso[i] + w[i+1] * y_iso[i+1]) / (w[i] + w[i+1])
    y_iso[i] <- pool_value
    y_iso[i+1] <- pool_value
    w[i] <- w[i] + w[i+1]
    w[i+1] <- w[i]
  }
  
  result <- y_iso[order(ord)]
  return(result)
}

binned_isotonic <- function(actuals, predicted, bin_width = 200) {
  binned_data <- data.frame(actuals, predicted) %>%
    mutate(bin = floor(predicted / bin_width) * bin_width) %>%
    group_by(bin) %>%
    summarise(
      bin_predicted = mean(predicted),
      bin_actuals = mean(actuals),
      .groups = "drop"
    )
  
  iso_model <- isotone::gpava(binned_data$bin_predicted, binned_data$bin_actuals)
  fitted <- approx(
    x = binned_data$bin_predicted, 
    y = iso_model$x, 
    xout = predicted, 
    rule = 2
  )$y
  
  return(fitted)
}

gam_calibration <- function(actuals, predicted) {
  gam_model <- mgcv::gam(actuals ~ s(predicted), family = gaussian())
  calibrated <- predict(gam_model, newdata = data.frame(predicted = predicted))
  return(calibrated)
}

kde_binned_isotonic <- function(actuals, predicted, bandwidth = 0.1) {
  density_est <- density(predicted, bw = bandwidth)
  cutpoints <- density_est$x[which(diff(sign(diff(density_est$y))) == -2)]  # Peaks
  
  bins <- cut(predicted, breaks = c(-Inf, sort(cutpoints), Inf), labels = FALSE)
  
  binned_data <- data.frame(actuals, predicted, bins) %>%
    group_by(bins) %>%
    summarise(
      bin_predicted = mean(predicted),
      bin_actuals = mean(actuals),
      .groups = "drop"
    )
  
  iso_model <- isotone::gpava(binned_data$bin_predicted, binned_data$bin_actuals)
  calibrated <- approx(
    x = binned_data$bin_predicted,
    y = iso_model$x,
    xout = predicted,
    rule = 2
  )$y
  
  return(calibrated)
}

# Compute
compute_calibration_metrics <- function(actuals, predicted, calibrated, method) {
  mace <- mean(abs(predicted - actuals))
  medace <- median(abs(predicted - actuals))
  spearman_before <- cor(predicted, actuals, method = "spearman")
  spearman_after <- cor(calibrated, actuals, method = "spearman")
  
  mace_cal <- mean(abs(calibrated - actuals))
  medace_cal <- median(abs(calibrated - actuals))
  
  tibble(
    Method = method,
    MACE_Before = mace,
    MedACE_Before = medace,
    Spearman_Before = spearman_before,
    MACE_After = mace_cal,
    MedACE_After = medace_cal,
    Spearman_After = spearman_after,
    Unique_Values_Before = length(unique(predicted)),
    Unique_Values_After = length(unique(calibrated))
  )
}

# Simulate data

set.seed(123)

simulate_data <- function(n = 10000) {
  actuals <- c(rep(1, 0.05 * n), runif(0.95 * n, min = 0, max = 300000))
  noise <- rnorm(n, mean = 0, sd = 50000)
  predicted <- 1 / (1 + exp(-(actuals - 150000) / 50000)) * 300000 + noise
  predicted <- pmax(predicted, 0) # Ensure no negative predictions
  data.frame(actuals = actuals, predicted = predicted)
}

data <- simulate_data()

# Estimate
# Isotonic calibration
calibrated_isotonic <- isotone::gpava(y = data$predicted, z = data$actuals)$x

# Robust PAVA calibration
calibrated_robust <- robust_pava(data$actuals, data$predicted)

# Binned isotonic calibration
calibrated_binned <- binned_isotonic(data$actuals, data$predicted)

# GAM calibration
calibrated_gam <- gam_calibration(data$actuals, data$predicted)

# KDE Binned calibration
calibrated_kde <- kde_binned_isotonic(data$actuals, data$predicted)

# Generate metrics
metrics_isotonic <- compute_calibration_metrics(
  actuals = data$actuals,
  predicted = data$predicted,
  calibrated = calibrated_isotonic,
  method = "Isotonic"
)

metrics_robust <- compute_calibration_metrics(
  actuals = data$actuals,
  predicted = data$predicted,
  calibrated = calibrated_robust,
  method = "Robust PAVA"
)

metrics_binned <- compute_calibration_metrics(
  actuals = data$actuals,
  predicted = data$predicted,
  calibrated = calibrated_binned,
  method = "Binned Isotonic"
)

metrics_gam <- compute_calibration_metrics(
  actuals = data$actuals,
  predicted = data$predicted,
  calibrated = calibrated_gam,
  method = "GAM"
)

metrics_kde <- compute_calibration_metrics(
  actuals = data$actuals,
  predicted = data$predicted,
  calibrated = calibrated_kde,
  method = "KDE-Binned Isotonic"
)

# Combine results
final_metrics <- bind_rows(metrics_isotonic, metrics_robust, metrics_binned, metrics_gam, metrics_kde)

# Output results
print(final_metrics)
