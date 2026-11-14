# ROLLING-ORIGIN TUNING

# Calculate validation MSE.
# Input: Subject coefficients B, observed series, and forecast origin.
# Output: One-step mean squared forecast error for each subject.
prediction_loss <- function(B, series, origin) {
  # Input: Subject index k. Output: Mean squared forecast error across subject k's variables.
  vapply(seq_along(B), function(k) {
    err <- series[[k]][origin + 1L, ] - drop(series[[k]][origin, ] %*% B[[k]])
    mean(err^2)
  }, numeric(1))
}


# Penalty selection and final refit. (For Multi-VAR and Adaptive Multi-VAR)
# Optional pilot_fits supplies window-specific and final standard multi-VAR pilots.
# Input: Series, origins, penalty fractions and ratios, adaptive and solver options, optional pilots.
# Output: Fitted coefficients, selected penalty settings, and validation MSE.
tune_multivar <- function(series, origins, fractions, ratios, adaptive = FALSE,
                          custom_final = FALSE, pilot_fits = NULL) {
  d <- ncol(series[[1]])
  K <- length(series)
  score <- matrix(0, length(ratios), length(fractions))
  for (origin in origins) {
    dat <- lagged(series, origin)
    aug <- augmented_design(dat)
    aw <- if (adaptive && is.null(pilot_fits)) adaptive_weights(dat) else NULL
    if (adaptive && !is.null(pilot_fits)) {
      aw <- weights_from_pilot(pilot_fits$validation[[as.character(origin)]]$B)
    }
    for (r in seq_along(ratios)) {
      W <- penalty_weights(d, K, ratios[r], aw)
      lam <- lambda_max(aug$Z, aug$Y, W) * fractions
      path <- glmnet_path(aug$Z, aug$Y, W, lam)
      for (ell in seq_along(fractions))
        score[r, ell] <- score[r, ell] + mean(prediction_loss(
          theta_to_B(path[, , ell], d, K), series, origin)) / length(origins)
    }
  }
  best <- arrayInd(which.min(score), dim(score))
  r <- best[1L]
  ell <- best[2L]
  dat <- lagged(series, nrow(series[[1]]))
  aug <- augmented_design(dat)
  aw <- if (adaptive && is.null(pilot_fits)) adaptive_weights(dat) else NULL
  if (adaptive && !is.null(pilot_fits)) {
    aw <- weights_from_pilot(pilot_fits$final$B)
  }
  W <- penalty_weights(d, K, ratios[r], aw)
  lam <- lambda_max(aug$Z, aug$Y, W) * fractions[ell]
  if (custom_final) {
    theta <- fit_fista(dat, lam, W)$theta
  } else {
    theta <- glmnet_path(aug$Z, aug$Y, W, lam, thresh = 1e-11)[, , 1L]
  }
  list(B = theta_to_B(theta, d, K), theta = theta,
       lambda = lam, ratio = ratios[r], fraction = fractions[ell],
       cv_mse = score[r, ell])
}


# Fit standard multi-VAR pilots within the adaptive validation windows.
# Input: Training series, forecast origins, penalty grids, and the final pilot fit.
# Output: Pilot fits by validation origin and the supplied final pilot.
fit_multivar_pilots <- function(series, origins, fractions, ratios, final_fit) {
  # Input: Outer forecast origin. Output: Pilot fit using data through that origin.
  validation <- lapply(origins, function(origin) {
    # Input: One subject's training series x. Output: Observations through the outer forecast origin.
    prefix <- lapply(series, function(x) x[seq_len(origin), , drop = FALSE])
    inner_origins <- unique(round(seq(
      max(5L, floor(origin / 3)), origin - 1L, length.out = length(origins)
    )))
    tune_multivar(prefix, inner_origins, fractions, ratios, custom_final = TRUE)
  })
  list(validation = setNames(validation, as.character(origins)), final = final_fit)
}


# Penalty selection and final refit. (For Individual LASSO and Pooled LASSO)
# Input: Training series, forecast origins, penalty fractions, and pooling option.
# Output: Subject coefficient matrices, mean selected fraction, and mean validation MSE.
tune_lasso <- function(series, origins, fractions, pooled = FALSE) {
  K <- length(series)
  d <- ncol(series[[1]])
  groups <- if (pooled) list(seq_len(K)) else as.list(seq_len(K))
  B <- vector("list", K)
  chosen <- numeric(length(groups))
  cv <- numeric(length(groups))
  for (g in seq_along(groups)) {
    ids <- groups[[g]]
    scores <- numeric(length(fractions))
    for (origin in origins) {
      dat <- lagged(series[ids], origin)
      X <- do.call(rbind, dat$X)
      Y <- do.call(rbind, dat$Y)
      W <- matrix(1, d, d)
      path <- glmnet_path(X, Y, W, lambda_max(X, Y, W) * fractions)
      for (ell in seq_along(fractions))
        scores[ell] <- scores[ell] + mean(prediction_loss(
          rep(list(path[, , ell]), length(ids)), series[ids], origin)) / length(origins)
    }
    ell <- which.min(scores)
    chosen[g] <- fractions[ell]
    cv[g] <- scores[ell]
    dat <- lagged(series[ids], nrow(series[[1]]))
    X <- do.call(rbind, dat$X)
    Y <- do.call(rbind, dat$Y)
    W <- matrix(1, d, d)
    coefs <- glmnet_path(X, Y, W, lambda_max(X, Y, W) * fractions[ell], thresh = 1e-11)[, , 1L]
    for (k in ids) B[[k]] <- coefs
  }
  list(B = B, fraction = mean(chosen), cv_mse = mean(cv))
}
