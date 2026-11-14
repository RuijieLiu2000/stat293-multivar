# REGRESSION DESIGN

# Form one-lag regression pairs.
# Input: Subject series and the last time point to include.
# Output: Lagged predictors X and next-step responses Y for each subject.
lagged <- function(series, end) {
  # Input: One subject's series x. Output: Predictor rows through time end - 1.
  list(X = lapply(series, function(x) x[seq_len(end - 1L), , drop = FALSE]),
       # Input: One subject's series x. Output: Response rows from time 2 through end.
       Y = lapply(series, function(x) x[2:end, , drop = FALSE]))
}


# Build the block regression design.
# Input: Lagged predictors and responses for each subject.
# Output: Block design Z, stacked responses Y, and counts d and K.
augmented_design <- function(dat) {
  pooled <- do.call(rbind, dat$X)
  Z <- cbind(Matrix::Matrix(pooled, sparse = TRUE),
             # Input: One subject's predictor matrix x. Output: Sparse copy of x.
             Matrix::bdiag(lapply(dat$X, function(x) Matrix::Matrix(x, sparse = TRUE))))
  list(Z = Z, Y = do.call(rbind, dat$Y), d = ncol(pooled), K = length(dat$X))
}


# Build weight matrix.
# Input: Variable and subject counts, deviation penalty ratio, and optional adaptive weights.
# Output: Matrix of shared and subject-specific penalty weights.
penalty_weights <- function(d, K, ratio, adaptive = NULL) {
  W <- if (is.null(adaptive)) matrix(1, (K + 1L) * d, d) else adaptive
  W[-seq_len(d), ] <- W[-seq_len(d), ] * ratio
  W
}


# Find the penalty bound(every coefficient is zero)
# Input: Design matrix Z, responses Y, and positive penalty weights W.
# Output: Smallest penalty for which all coefficients can be zero.
lambda_max <- function(Z, Y, W) max(abs(as.matrix(crossprod(Z, Y))) / nrow(Z) / W)


# Combine shared and individual blocks.
# Input: Stacked coefficients theta, variable count d, and subject count K.
# Output: One fitted coefficient matrix per subject.
theta_to_B <- function(theta, d, K) {
  # Input: Subject index k. Output: Shared coefficients plus subject k's deviations.
  lapply(seq_len(K), function(k)
    theta[seq_len(d), , drop = FALSE] + theta[k * d + seq_len(d), , drop = FALSE])
}
