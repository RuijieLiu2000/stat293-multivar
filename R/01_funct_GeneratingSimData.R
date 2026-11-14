# DATA GENERATION

# Compute the largest eigenvalue for stable transition.
# Input: Square transition matrix A.
# Output: Largest absolute eigenvalue of A.
spectral_radius <- function(A) max(Mod(eigen(A, only.values = TRUE)$values))


# Generate shared and individual VAR(1) coefficients.
# Input: Variable and subject counts, nonzero and shared fractions, target spectral radius.
# Output: Shared, subject-specific, and total transition matrices.
make_truth <- function(d, K, density = .05, common_fraction = .5,
                       target_radius = .8) {
  s <- max(2L, as.integer(round(density * d^2)))
  sc <- max(1L, min(s - 1L, as.integer(round(common_fraction * s))))
  su <- s - sc
  for (attempt in seq_len(1000L)) {
    positions <- sample.int(d^2, sc + K * su)
    common <- matrix(0, d, d)
    common[positions[seq_len(sc)]] <- runif(sc, .1, .9)
    # Input: Subject index k. Output: Unscaled deviation matrix for subject k.
    unique <- lapply(seq_len(K), function(k) {
      U <- matrix(0, d, d)
      U[positions[sc + (k - 1L) * su + seq_len(su)]] <- runif(su, .1, .9)
      U
    })
    # Input: One subject's deviation matrix U. Output: Shared matrix plus U.
    A <- lapply(unique, function(U) common + U)
    r <- max(vapply(A, spectral_radius, numeric(1)))
    if (r > .05) break
  }
  if (r <= .05) stop("Not Stable.")
  scale_factor <- target_radius / r
  common <- common * scale_factor
  # Input: One subject's deviation matrix U. Output: U scaled by the common stability factor.
  unique <- lapply(unique, function(U) U * scale_factor)
  # Input: One subject's scaled deviation matrix U. Output: Scaled total transition matrix.
  A <- lapply(unique, function(U) common + U)
  list(A = A, common = common, unique = unique)
}


# Simulate latent and observed series.
# Input: True matrices, series length, noise scales, error correlations, and burn-in length.
# Output: Observed series, latent series, and innovations for each subject.
simulate_data <- function(truth, T = 100L, sigma = 1, cross_rho = 0,
                          serial_rho = 0, measurement_sd = 0, burn = 500L) {
  d <- nrow(truth$common)
  K <- length(truth$A)
  # Input: Variable indices i and j. Output: Same-time error correlation for that pair.
  R <- outer(seq_len(d), seq_len(d), function(i, j) cross_rho^abs(i - j))
  C <- chol(R)
  latent <- observed <- innovations <- vector("list", K)
  for (k in seq_len(K)) {
    z <- matrix(rnorm((burn + T) * d), burn + T, d) %*% C
    e <- sigma * drop(matrix(rnorm(d), 1L) %*% C)
    x <- numeric(d)
    xx <- ee <- matrix(0, burn + T, d)
    for (tt in seq_len(burn + T)) {
      e <- serial_rho * e + sigma * sqrt(1 - serial_rho^2) * z[tt, ]
      x <- drop(truth$A[[k]] %*% x) + e
      xx[tt, ] <- x
      ee[tt, ] <- e
    }
    latent[[k]] <- xx[burn + seq_len(T), , drop = FALSE]
    innovations[[k]] <- ee[burn + seq_len(T), , drop = FALSE]
    observed[[k]] <- latent[[k]] + measurement_sd * matrix(rnorm(T * d), T, d)
    colnames(observed[[k]]) <- colnames(latent[[k]]) <- paste0("V", seq_len(d))
  }
  list(observed = observed, latent = latent, innovations = innovations)
}
