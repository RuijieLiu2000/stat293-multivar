# MODEL FITTING

# Construct adaptive weights from preliminary subject coefficient matrices.
# Input: Subject pilot coefficient matrices B and denominator floor epsilon.
# Output: Stacked penalty weights for shared coefficients and deviations.
weights_from_pilot <- function(B, epsilon = 1e-3) {
  d <- nrow(B[[1]])
  med <- matrix(apply(simplify2array(B), c(1, 2), median), d, d)
  do.call(rbind, c(list(1 / pmax(abs(med), epsilon)),
                   # Input: One subject's pilot matrix b. Output: Weights for deviations from the pilot median.
                   lapply(B, function(b) 1 / pmax(abs(b - med), epsilon))))
}

# Preliminary-estimation and construct weights.
# Input: Lagged subject data and the denominator floor epsilon.
# Output: Adaptive penalty weights estimated from the pilot fits.
adaptive_weights <- function(dat, epsilon = 1e-3) {
  d <- ncol(dat$X[[1]])
  K <- length(dat$X)
  # Input: One subject's predictor matrix x. Output: Whether x has full column rank and more rows than columns.
  ols_ok <- all(vapply(dat$X, function(x) nrow(x) > d && qr(x)$rank == d, logical(1)))
  if (ols_ok) {
    # Input: Subject index k. Output: OLS coefficient matrix for subject k.
    B <- lapply(seq_len(K), function(k) qr.solve(dat$X[[k]], dat$Y[[k]]))
  } else {
    aug <- augmented_design(dat)
    W <- penalty_weights(d, K, 2 / K)
    lam <- .1 * lambda_max(aug$Z, aug$Y, W)
    theta <- glmnet_path(aug$Z, aug$Y, W, lam)[, , 1L]
    B <- theta_to_B(theta, d, K)
  }
  weights_from_pilot(B, epsilon)
}

# Fit weighted LASSO.
# Input: Design Z, responses Y, penalty weights W, penalty values, and solver tolerance.
# Output: Coefficient array indexed by predictor, response, and penalty.
glmnet_path <- function(Z, Y, W, lambdas, thresh = 1e-7) {
  q <- ncol(Z)
  d <- ncol(Y)
  L <- length(lambdas)
  out <- array(0, c(q, d, L))
  for (j in seq_len(d)) {
    fit <- glmnet::glmnet(Z, Y[, j], family = "gaussian", alpha = 1,
                         lambda = lambdas * mean(W[, j]),
                         penalty.factor = W[, j], intercept = FALSE,
                         standardize = FALSE, thresh = thresh, maxit = 1000000L)
    if (fit$jerr != 0L || length(fit$lambda) != L)
      stop("glmnet fails")
    out[, j, ] <- as.matrix(fit$beta)
  }
  out
}


# Calculate the objective.
# Input: Lagged data, stacked coefficients theta, weights W, and penalty lambda.
# Output: Scaled squared-error loss plus the weighted LASSO penalty.
objective <- function(dat, theta, W, lambda) {
  d <- ncol(dat$X[[1]])
  K <- length(dat$X)
  B <- theta_to_B(theta, d, K)
  n <- sum(vapply(dat$X, nrow, integer(1)))
  # Input: Subject index k. Output: Sum of squared residuals for subject k.
  sum(vapply(seq_len(K), function(k)
    sum((dat$Y[[k]] - dat$X[[k]] %*% B[[k]])^2), numeric(1))) / (2 * n) +
    lambda * sum(W * abs(theta))
}


# Fit multi-VAR using accelerated proximal gradient descent.
# Input: Lagged data, penalty lambda, weights W, tolerance, and iteration limit.
# Output: Stacked coefficients theta and fitted matrices B for each subject.
fit_fista <- function(dat, lambda, W, tol = 1e-6, maxit = 12000L) {
  d <- ncol(dat$X[[1]])
  K <- length(dat$X)
  n <- sum(vapply(dat$X, nrow, integer(1)))
  # Input: One subject's predictor matrix x. Output: Predictor cross-product scaled by the total row count.
  G <- lapply(dat$X, function(x) crossprod(x) / n)
  # Input: Subject index k. Output: Predictor-response cross-product scaled by the total row count.
  H <- lapply(seq_len(K), function(k) crossprod(dat$X[[k]], dat$Y[[k]]) / n)
  # Input: One subject's scaled predictor cross-product g. Output: Largest eigenvalue of g.
  L <- (K + 1) * max(vapply(G, function(g)
    max(eigen(g, symmetric = TRUE, only.values = TRUE)$values), numeric(1)))
  # Input: Stacked coefficient matrix theta.
  # Output: Gradient of the squared-error loss for all coefficient blocks.
  grad <- function(theta) {
    # Input: Subject index k. Output: Subject k's contribution to the gradient.
    D <- lapply(seq_len(K), function(k)
      G[[k]] %*% (theta[seq_len(d), , drop = FALSE] +
                     theta[k * d + seq_len(d), , drop = FALSE]) - H[[k]])
    do.call(rbind, c(list(Reduce(`+`, D)), D))
  }
  # Input: Coefficient values x and thresholds t.
  # Output: Soft-thresholded values.
  soft <- function(x, t) sign(x) * pmax(abs(x) - t, 0)
  theta <- extrap <- matrix(0, (K + 1L) * d, d)
  a <- 1
  obj <- objective(dat, theta, W, lambda)
  grad_scale <- max(abs(grad(theta)), .Machine$double.eps)
  converged <- FALSE
  for (iteration in seq_len(maxit)) {
    candidate <- soft(extrap - grad(extrap) / L, lambda * W / L)
    new_obj <- objective(dat, candidate, W, lambda)
    if (new_obj > obj + 1e-12 * max(1, abs(obj))) {
      extrap <- theta
      a <- 1
      candidate <- soft(extrap - grad(extrap) / L, lambda * W / L)
      new_obj <- objective(dat, candidate, W, lambda)
    }
    a_new <- (1 + sqrt(1 + 4 * a^2)) / 2
    extrap <- candidate + ((a - 1) / a_new) * (candidate - theta)
    theta <- candidate
    a <- a_new
    obj <- new_obj
    if (iteration %% 10L == 0L) {
      g <- grad(theta)
      violation <- ifelse(abs(theta) > 1e-10,
                          abs(g + lambda * W * sign(theta)),
                          pmax(abs(g) - lambda * W, 0))
      kkt_relative <- max(violation) / grad_scale
      if (kkt_relative < tol) {
        converged <- TRUE
        break
      }
    }
  }
  if (!converged) stop("FISTA fails")
  list(theta = theta, B = theta_to_B(theta, d, K))
}



