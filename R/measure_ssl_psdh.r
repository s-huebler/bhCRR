#' measure_ssl_psdh
#'
#' @param y_true
#' @param risk_score
#' @param evaluation_time
#'
#' @returns
#' @export
#'
#' @examples
measure_ssl_psdh <- function(y_true, risk_score, evaluation_time) {

  # Extract times and status from y_true array with (TTE, Status)
  # 1=Event, 2=Competing, 0=Censored
  T_i <- y_true[, 1]
  D_i <- y_true[, 2]
  n <- length(T_i)

  # --- Step 1: Estimate Censoring Distribution G(t) ---
  # Switch status: 0 becomes event (censoring), 1&2 become censored
  cens_status <- ifelse(D_i == 0, 1, 0)

  # Compute Kaplan-Meier for Censoring (G)
  # Need G(t) and G(t-) for every observed time T_i
  km_cens <- survfit(Surv(T_i, cens_status) ~ 1)

  # Helper to get G(t) from KM object
  get_G <- function(t, km_obj) {
    # Find the index of the largest time in km_obj$time <= t
    idx <- findInterval(t, km_obj$time)
    # If t is smaller than earliest observed time, prob is 1
    probs <- rep(1, length(t))
    valid <- idx > 0
    if (any(valid)) {
      probs[valid] <- km_obj$surv[idx[valid]]
    }
    return(probs)
  }

  # Pre-calculate G(T_i) and G(T_i - epsilon) for weights
  # G_hat_Ti  <- G(T_i)
  # G_hat_Tim <- G(T_i-)
  G_Ti  <- get_G(T_i, km_cens)
  # For G(T_i-), we effectively look up T_i - epsilon
  # In findInterval, this is equivalent to looking for strict inequality
  # or simply using the previous index in step function.
  # Simpler approach: findInterval with left.open=TRUE logic manually
  idx_minus <- findInterval(T_i - 1e-08, km_cens$time)
  G_Ti_minus <- rep(1, n)
  valid_m <- idx_minus > 0
  G_Ti_minus[valid_m] <- km_cens$surv[idx_minus[valid_m]]


  # --- Step 2: Calculate Numerator and Denominator (Eq 3.4) ---
  numerator <- 0
  denominator <- 0

  # Iterate over all pairs (i, j)
  for (i in 1:n) {

    # We only care about subject i if they are a "Case" at time t
    # Def: N_i^1(t) = I(T_i <= t, D_i = 1)
    if (!(T_i[i] <= evaluation_time && D_i[i] == 1)) {
      next
    }

    for (j in 1:n) {
      if (i == j) next

      # Determine if pair (i,j) is evaluable
      # A_ij = I(T_i < T_j)
      # B_ij = I(T_i >= T_j and D_j = 2) (Competing event happens before or at T_i)


      is_A <- T_i[i] < T_i[j]
      is_B <- (T_i[i] >= T_i[j]) && (D_i[j] == 2)

      if (!is_A && !is_B) {
        # Pair is not comparable/evaluable
        next
      }

      # Calculate Weight W_ij
      # W_ij,1 = G(T_i-) * G(T_i)  [If is_A]
      # W_ij,2 = G(T_i-) * G(T_j-) [If is_B]

      weight <- 0
      if (is_A) {
        # Check divisibility to avoid NaNs
        w_val <- G_Ti_minus[i] * G_Ti[i]
        if (w_val > 0) weight <- 1 / w_val
      } else if (is_B) {
        w_val <- G_Ti_minus[i] * G_Ti_minus[j]
        if (w_val > 0) weight <- 1 / w_val
      }

      if (weight == 0) next

      # Check Concordance
      # Q_ij(t) = I(M_i > M_j)
      # Higher risk score for 'i' (the one with the event) is concordant
      concordant <- 0
      if (risk_score[i] > risk_score[j]) {
        concordant <- 1
      } else if (risk_score[i] == risk_score[j]) {
        concordant <- 0.5 # Handle ties
      }

      # Accumulate
      numerator   <- numerator + (weight * concordant)
      denominator <- denominator + weight
    }
  }

  if (denominator == 0) return(NA)

  return(numerator / denominator)
}
