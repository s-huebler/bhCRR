#' IPCW concordance index for competing risks
#'
#' Computes the inverse-probability-of-censoring-weighted (IPCW) C-index
#' for a cause-1 competing-risks model at a specified evaluation time.  The
#' censoring distribution \eqn{G(t)} is estimated by Kaplan-Meier; pairs
#' \eqn{(i, j)} are deemed comparable if subject \eqn{i} had a cause-1
#' event before the evaluation time and either (A) subject \eqn{j}'s
#' event/censoring time is later, or (B) subject \eqn{j} had a competing
#' event before subject \eqn{i}.
#'
#' @param y_true Numeric matrix of dimensions \eqn{n \times 2}.  Column 1
#'   contains observed event/censoring times; column 2 contains status codes
#'   (\code{0} = censored, \code{1} = cause-1 event, \code{2} = competing
#'   event).
#' @param risk_score Numeric vector of length \eqn{n}.  Predicted absolute
#'   risk (cumulative incidence) at \code{evaluation_time} for each subject,
#'   as returned by \code{\link{predict_from_ssl_psdh}}.  \code{NA} values
#'   are skipped when evaluating concordance for a pair.
#' @param evaluation_time Numeric scalar.  Time horizon at which the C-index
#'   is evaluated; only subjects with a cause-1 event at or before this time
#'   contribute as index cases.
#'
#' @returns Numeric scalar in \eqn{[0, 1]}.  The IPCW-weighted C-index;
#'   \code{NA} if no comparable pairs exist (e.g. all predictions are
#'   missing or the denominator is zero).
#'
#' @importFrom survival survfit Surv
#'
#' @seealso \code{\link{predict_from_ssl_psdh}}, \code{\link{cv_ssl_psdh}}
#'
#' @export
#'
#' @examples
#' \dontrun{
#' lp  <- predict_from_ssl_psdh(fit, newx = x, prediction_time = 30)
#' measure_ssl_psdh(y_true = y, risk_score = lp, evaluation_time = 30)
#' }
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
  km_cens <- survival::survfit(survival::Surv(T_i, cens_status) ~ 1)

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
  G_Ti <- get_G(T_i, km_cens)
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
      if (i == j) {
        next
      }

      if (is.na(risk_score[i]) | is.na(risk_score[j])) {
        next
      }

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

      if (weight == 0) {
        next
      }

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
      numerator <- numerator + (weight * concordant)
      denominator <- denominator + weight
    }
  }

  if (denominator == 0) {
    return(NA)
  }

  return(numerator / denominator)
}
