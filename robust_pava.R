
# Traditional PAVA
traditional_pava <- function(y) {
  n <- length(y)
  y_iso <- y
  w <- rep(1, n)  # weights
  
  while(TRUE) {
    # Find violations
    violations <- which(diff(y_iso) < 0)
    if(length(violations) == 0) break
    
    # Pool adjacent violators
    i <- violations[1]
    # Weighted average for the violating pair
    pool_value <- (w[i] * y_iso[i] + w[i+1] * y_iso[i+1]) / (w[i] + w[i+1])
    y_iso[i] <- pool_value
    y_iso[i+1] <- pool_value
    # Update weights
    w[i] <- w[i] + w[i+1]
    w[i+1] <- w[i]
  }
  
  return(y_iso)
}

# Determine threshold from data
# 5th percentile of difference
determine_threshold <- function(y, percentile = 0.05) {
  sorted_y <- sort(y)
  
  diffs <- diff(sorted_y)
  
  positive_diffs <- diffs[diffs > 0]
  if(length(positive_diffs) == 0) return(0)
  
  threshold <- quantile(positive_diffs, percentile)
  return(threshold)
}

# Robust PAVA with threshold
robust_pava <- function(y, threshold = 0.1) {
  n <- length(y)
  y_iso <- y
  w <- rep(1, n)  # weights
  
  while(TRUE) {
    # Find violations
    diffs <- diff(y_iso)
    violations <- which(diffs < 0)
    # Only consider violations greater than threshold
    violations <- violations[abs(diffs[violations]) > threshold]
    
    if(length(violations) == 0) break
    
    # Pool adjacent violators
    i <- violations[1]
    pool_value <- (w[i] * y_iso[i] + w[i+1] * y_iso[i+1]) / (w[i] + w[i+1])
    y_iso[i] <- pool_value
    y_iso[i+1] <- pool_value
    w[i] <- w[i] + w[i+1]
    w[i+1] <- w[i]
  }
  
  return(y_iso)
}

# MACE
calculate_mace <- function(predicted_probs, actual_outcomes) {
  n <- length(predicted_probs)
  n_bins <- min(n %/% 50, 20)  # At least 50 samples per bin, max 20 bins
  
  # Sort by predictions
  ord <- order(predicted_probs)
  sorted_preds <- predicted_probs[ord]
  sorted_actual <- actual_outcomes[ord]
  
  # Create roughly equal-sized bins
  bin_size <- ceiling(n / n_bins)
  bins <- rep(1:n_bins, each = bin_size)[1:n]
  
  # Calculate mean prediction and actual outcome for each bin
  bin_means <- tapply(sorted_preds, bins, mean)
  bin_actuals <- tapply(sorted_actual, bins, mean)
  
  # MACE is mean absolute difference between predictions and actuals
  mace <- mean(abs(bin_means - bin_actuals))
  return(mace)
}
