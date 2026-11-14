# MODEL EVALUATION AND MONTE CARLO SUMMARIES

# Measure coefficient recovery and forecasting.
# Input: Fitted model, true matrices, simulated data, innovation SD, and selection threshold.
# Output: Subject-averaged coefficient, support, and forecast accuracy measures.
evaluate_fit <- function(fit, truth, sim, sigma, threshold = 1e-5) {
  K <- length(truth$A)
  T <- nrow(sim$observed[[1]])
  train_end <- T - 3L
  relative <- sensitivity <- specificity <- precision <- numeric(K)
  forecast <- matrix(0, K, 3L)
  for (k in seq_len(K)) {
    A <- truth$A[[k]]
    estimate <- t(fit$B[[k]])
    relative[k] <- sqrt(sum((estimate - A)^2) / sum(A^2))
    selected <- abs(estimate) > threshold
    nonzero <- A != 0
    sensitivity[k] <- mean(selected[nonzero])
    specificity[k] <- mean(!selected[!nonzero])
    precision[k] <- if (sum(selected)) mean(nonzero[selected]) else NA_real_
    pred <- sim$observed[[k]][train_end, ]
    for (h in 1:3) {
      pred <- drop(pred %*% fit$B[[k]])
      forecast[k, h] <- sqrt(mean((pred - sim$observed[[k]][train_end + h, ])^2))
    }
  }
  c(relative_frobenius = mean(relative),
    sensitivity = mean(sensitivity), specificity = mean(specificity),
    precision = if (all(is.na(precision))) NA_real_ else mean(precision, na.rm = TRUE),
    setNames(colMeans(forecast), paste0("rmsfe_h", 1:3)),
    setNames(colMeans(forecast) / sigma, paste0("rmsfe_over_sigma_h", 1:3)))
}


# Average metric across simulation replicates.
# Input: Replicate results, metric names, and scenario design.
# Output: Means, Monte Carlo SEs, counts, and 95% intervals by scenario and method.
summarize_results <- function(raw, metrics, scenarios) {
  rows <- list()
  b <- 0L
  for (s in unique(raw$scenario)) for (method in unique(raw$method)) {
    z <- raw[raw$scenario == s & raw$method == method, ]
    for (metric in metrics) {
      x <- z[[metric]]
      x <- x[is.finite(x)]
      n <- length(x)
      se <- if (n > 1L) sd(x) / sqrt(n) else NA_real_
      critical <- if (n > 1L) qt(.975, n - 1L) else NA_real_
      b <- b + 1L
      rows[[b]] <- data.frame(scenario = s, method = method, metric = metric,
                              mean = mean(x), mcse = se, n = n,
                              lower = mean(x) - critical * se,
                              upper = mean(x) + critical * se)
    }
  }
  merge(do.call(rbind, rows), scenarios, by = "scenario", sort = FALSE)
}


# Metric comparison for different scenarios.
# Input: Replicate results, metrics, and reference (baseline scenario or individual LASSO).
# Output: Mean paired differences, Monte Carlo SEs, and replicate counts.
paired_contrasts <- function(raw, metrics, scenario_comparison = FALSE) {
  out <- list()
  i <- 0L
  for (metric in metrics) {
    if (scenario_comparison) {
      base <- raw[raw$scenario == "baseline", c("replication", "method", metric)]
      names(base)[3L] <- "reference"
      joined <- merge(raw, base, by = c("replication", "method"))
      groups <- split(joined, interaction(joined$scenario, joined$method, drop = TRUE))
    } else {
      base <- raw[raw$method == "individual_lasso", c("replication", "scenario", metric)]
      names(base)[3L] <- "reference"
      joined <- merge(raw, base, by = c("replication", "scenario"))
      groups <- split(joined, interaction(joined$scenario, joined$method, drop = TRUE))
    }
    for (z in groups) {
      x <- z[[metric]] - z$reference
      x <- x[is.finite(x)]
      i <- i + 1L
      out[[i]] <- data.frame(scenario = z$scenario[1L], method = z$method[1L],
                             metric = metric, difference = mean(x),
                             mcse = if (length(x) > 1L) sd(x) / sqrt(length(x)) else NA_real_,
                             n = length(x))
    }
  }
  do.call(rbind, out)
}
