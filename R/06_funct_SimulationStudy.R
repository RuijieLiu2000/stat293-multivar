# STUDY DESIGN

# Define the baseline and 15 one-factor alternatives.
# Input: None.
# Output: Table of settings for the 16 simulation conditions.
make_scenarios <- function() {
  base <- data.frame(scenario = "baseline", factor = "baseline", value = 0,
                     d = 20L, K = 10L, T = 100L, density = .05,
                     common_fraction = .5, sigma = 1, cross_rho = 0,
                     serial_rho = 0, measurement_sd = 0, target_radius = .8)
  variations <- list(
    innovation_scale = list(column = "sigma", values = c(.5, 2)),
    dimension = list(column = "d", values = c(10L, 30L, 40L)),
    density = list(column = "density", values = c(.025, .10, .15)),
    contemporaneous_correlation = list(column = "cross_rho", values = c(.3, .6)),
    serial_correlation = list(column = "serial_rho", values = c(.3, .6)),
    measurement_noise = list(column = "measurement_sd", values = c(.5, 1)),
    series_length = list(column = "T", values = 50L)
  )
  rows <- list(base)
  for (factor in names(variations)) {
    variation <- variations[[factor]]
    for (value in variation$values) {
      row <- base
      row$scenario <- paste0(factor, "_", value)
      row$factor <- factor
      row$value <- value
      row[[variation$column]] <- value
      rows[[length(rows) + 1L]] <- row
    }
  }
  do.call(rbind, rows)
}


# Computational settings: replication counts, seed, origins etc.
# Input: Study profile and replication count.
# Output: List of simulation and tuning settings.
study_config <- function(profile = "preliminary", reps = 5L) {
  fine_grid <- profile == "thorough"
  list(
    profile = profile,
    reps = reps,
    seed = 20261111L,
    fractions = exp(seq(0, log(.001), length.out = if (fine_grid) 20L else 10L)),
    cv_origins = if (fine_grid) Inf else 5L,
    ratio_multipliers = if (fine_grid) c(.5, 1.1, 2, 4, 8, 16) else c(.5, 2, 8),
    burn = 500L,
    support_threshold = 1e-5
  )
}


# Choose expanding-window validation origins used for choosing penalty
# Input: Series length and number of validation origins (Inf for all).
# Output: Time indices for one-step validation forecasts.
rolling_origins <- function(total_times, count) {
  eligible <- seq.int(max(5L, floor(total_times / 3)), total_times - 4L)
  if (is.infinite(count)) return(eligible)
  unique(round(seq(min(eligible), max(eligible), length.out = count)))
}


# Generate one dataset with random seeds across scenarios.
# Input: One scenario row, replication number, and study settings.
# Output: True coefficient matrices and simulated subject series.
generate_dataset <- function(scenario, replication, config) {
  seed <- config$seed + replication
  set.seed(seed)
  truth <- make_truth(
    scenario$d, scenario$K, scenario$density,
    scenario$common_fraction, scenario$target_radius
  )
  set.seed(seed + 100000L)
  sim <- simulate_data(
    truth, T = scenario$T, sigma = scenario$sigma,
    cross_rho = scenario$cross_rho, serial_rho = scenario$serial_rho,
    measurement_sd = scenario$measurement_sd, burn = config$burn
  )
  list(truth = truth, sim = sim)
}


# Fit the four estimators on same training dataset.
# Input: Training series, validation origins, and study settings.
# Output: Named list of fits for the four methods.
fit_methods <- function(training, origins, config) {
  ratios <- config$ratio_multipliers / length(training)
  list(
    multiVAR_FISTA = tune_multivar(
      training, origins, config$fractions, ratios, custom_final = TRUE
    ),
    adaptive_multiVAR = tune_multivar(
      training, origins, config$fractions, ratios, adaptive = TRUE
    ),
    individual_lasso = tune_lasso(training, origins, config$fractions),
    pooled_lasso = tune_lasso(training, origins, config$fractions, pooled = TRUE)
  )
}


# Simulate, fit, and evaluate one replicate of a condition.
# Input: One scenario row, replication number, and study settings.
# Output: One row of accuracy measures per method.
run_replication <- function(scenario, replication, config) {
  dataset <- generate_dataset(scenario, replication, config)
  train_end <- scenario$T - 3L
  # Input: One subject's observed series x. Output: Training observations through train_end.
  training <- lapply(dataset$sim$observed, function(x) {
    x[seq_len(train_end), , drop = FALSE]
  })
  origins <- rolling_origins(scenario$T, config$cv_origins)
  fits <- fit_methods(training, origins, config)
  # Input: Fitted method name. Output: One row of replication details and accuracy measures.
  rows <- lapply(names(fits), function(method) {
    metrics <- evaluate_fit(
      fits[[method]], dataset$truth, dataset$sim,
      scenario$sigma, config$support_threshold
    )
    cbind(
      data.frame(scenario = scenario$scenario, replication = replication,
                 seed = config$seed + replication, method = method),
      as.data.frame(as.list(metrics))
    )
  })
  do.call(rbind, rows)
}


# Run the full simulation study and export figures and tables.
# Input: Study profile, replication count, optional scenario names, and output folder.
# Output: Replicate results, summaries, scenarios, and settings. Saves figures and tables.
run_study <- function(profile = "preliminary", reps = 5L, selected = NULL,
                      output_dir = "results") {
  config <- study_config(profile, reps)
  scenarios <- make_scenarios()
  if (!is.null(selected)) {
    scenarios <- scenarios[scenarios$scenario %in% selected, ]
  }
  rows <- list()
  for (scenario_index in seq_len(nrow(scenarios))) {
    scenario <- scenarios[scenario_index, ]
    for (replication in seq_len(reps)) {
      message(sprintf("%s: replication %d/%d", scenario$scenario, replication, reps))
      rows[[length(rows) + 1L]] <- run_replication(scenario, replication, config)
    }
  }
  raw <- do.call(rbind, rows)
  metrics <- setdiff(names(raw), c("scenario", "replication", "seed", "method"))
  summary <- summarize_results(raw, metrics, scenarios)
  write_key_results(raw, summary, scenarios, output_dir)
  invisible(list(raw = raw, summary = summary, scenarios = scenarios, config = config))
}
