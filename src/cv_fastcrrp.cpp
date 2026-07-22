// RcppArmadillo backend for cross-validated tuning of fastCrrp models.
//
// This source file provides the heavy numerical kernels used by the R wrapper
// `cv_fastCrrp_cpp()`. The model fitting itself is still performed in R by
// `fastcmprsk::fastCrrp()` (an external, already-compiled routine); what is
// ported to C++ here is the cross-validation metric computation, which is the
// genuine O(n^2) bottleneck when the number of lambda values and/or the fold
// sizes are large.
//
// Two concordance measures are provided, mirroring the two `tuning` options of
// the pure-R `cv_fastCrrp()`:
//
//   * "normal"  -> harrell_cindex():  a standard right-censored C-index that
//                  matches survival::concordance(Surv(time, status==1) ~ rs)
//                  for distinct event times (see notes on ties below).
//   * "wolbers" -> ipcw_cindex():     an exact translation of wolbers_c(),
//                  the IPCW competing-risks C-index, with the Kaplan-Meier
//                  censoring distribution G(t) estimated once per fold.
//
// Author: generated for Sophia Huebler's bhCRR package (Project 2).

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>

using namespace Rcpp;

// ---------------------------------------------------------------------------
// Kaplan-Meier estimate of the censoring distribution G(t).
//
// Faithful to wolbers_c(): the "event" for the censoring KM is the
// original censoring indicator (status == 0). We return, for every subject i,
//   G_Ti[i]       = G(T_i)        (right-continuous step value at T_i)
//   G_Ti_minus[i] = G(T_i - eps)  (value just before T_i)
// matching the findInterval()-based lookups in the R reference.
// ---------------------------------------------------------------------------
static void km_censoring(const arma::vec& time,
                         const arma::vec& status,
                         arma::vec& G_Ti,
                         arma::vec& G_Ti_minus) {
  const arma::uword n = time.n_elem;

  // Unique observed times, sorted ascending. survfit() lists a row for every
  // unique observed time, so indexing into `ut` reproduces findInterval().
  arma::vec ut = arma::unique(time);      // sorted ascending
  const arma::uword m = ut.n_elem;

  arma::vec surv(m);
  double S = 1.0;
  for (arma::uword k = 0; k < m; ++k) {
    const double tk = ut[k];
    arma::uword nrisk = 0;   // # at risk: T_i >= tk
    arma::uword d     = 0;   // # censoring "events" at tk: T_i == tk & status == 0
    for (arma::uword i = 0; i < n; ++i) {
      if (time[i] >= tk) ++nrisk;
      if (time[i] == tk && status[i] == 0.0) ++d;
    }
    if (nrisk > 0) S *= (1.0 - static_cast<double>(d) / static_cast<double>(nrisk));
    surv[k] = S;
  }

  const double eps = 1e-08;  // matches the R reference's T_i - 1e-08
  for (arma::uword i = 0; i < n; ++i) {
    const double t  = time[i];
    const double tm = t - eps;

    // findInterval(t, ut): number of ut entries <= t.
    arma::uword idx = 0;
    for (arma::uword k = 0; k < m; ++k) {
      if (ut[k] <= t) idx = k + 1; else break;
    }
    G_Ti[i] = (idx > 0) ? surv[idx - 1] : 1.0;

    arma::uword idxm = 0;
    for (arma::uword k = 0; k < m; ++k) {
      if (ut[k] <= tm) idxm = k + 1; else break;
    }
    G_Ti_minus[i] = (idxm > 0) ? surv[idxm - 1] : 1.0;
  }
}

//' Kaplan-Meier censoring distribution (C++ helper, for verification)
//'
//' Exposes the internal censoring-distribution estimate used by the IPCW
//' C-index so it can be checked against the R reference in
//' \code{wolbers_c()}.
//'
//' @param time Numeric vector of observed event/censoring times.
//' @param status Numeric vector of status codes (0 = censored, 1 = cause-1
//'   event, 2 = competing event).
//' @return A list with elements \code{G_Ti} = \eqn{G(T_i)} and
//'   \code{G_Ti_minus} = \eqn{G(T_i^-)}.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::List cpp_km_censoring(const Rcpp::NumericVector& time_,
                            const Rcpp::NumericVector& status_) {
  const arma::vec time   = Rcpp::as<arma::vec>(time_);
  const arma::vec status = Rcpp::as<arma::vec>(status_);
  const arma::uword n = time.n_elem;
  arma::vec G_Ti(n), G_Ti_minus(n);
  km_censoring(time, status, G_Ti, G_Ti_minus);
  // Build plain (1-D) numeric vectors so the result has no dim attribute.
  Rcpp::NumericVector g(G_Ti.begin(), G_Ti.end());
  Rcpp::NumericVector gm(G_Ti_minus.begin(), G_Ti_minus.end());
  return Rcpp::List::create(Rcpp::Named("G_Ti") = g,
                            Rcpp::Named("G_Ti_minus") = gm);
}

// ---------------------------------------------------------------------------
// Harrell-style right-censored C-index for a single risk-score vector.
//
// Orientation: `reverse = false` matches the survival::concordance() default,
// in which a pair is concordant when the subject with the LARGER risk score
// also has the LARGER time (i.e. the score behaves "protectively"). Set
// `reverse = true` to flip this (larger score <-> shorter time concordant),
// which is the usual orientation for a hazard/risk score.
//
// Only cause-1 events count as events: the event indicator is (status == 1);
// status 0 and 2 are both treated as censored, mirroring the R call
// survival::concordance(Surv(time, status == 1) ~ rs).
//
// Tie handling: tied *times* are treated as non-comparable (skipped). For
// distinct event times this reproduces survival::concordance() exactly; with
// tied event times survival applies a small partial-credit correction, so
// confirm equivalence on your data via the test in the handoff notes.
// ---------------------------------------------------------------------------
static double harrell_cindex(const arma::vec& time,
                             const arma::vec& status,
                             const arma::vec& risk,
                             bool reverse) {
  const arma::uword n = time.n_elem;
  double conc = 0.0, disc = 0.0, tied = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    for (arma::uword j = i + 1; j < n; ++j) {
      if (std::isnan(risk[i]) || std::isnan(risk[j])) continue;
      if (time[i] == time[j]) continue;  // tied times: not comparable

      // Identify earlier/later subject by time.
      arma::uword a, b;  // a = earlier time, b = later time
      if (time[i] < time[j]) { a = i; b = j; } else { a = j; b = i; }

      // Comparable only if the earlier subject had a cause-1 event.
      if (status[a] != 1.0) continue;

      const double ra = risk[a];
      const double rb = risk[b];

      if (ra == rb) {
        tied += 1.0;
      } else {
        // reverse = false: concordant when later time has the larger score.
        bool concordant = reverse ? (ra > rb) : (ra < rb);
        if (concordant) conc += 1.0; else disc += 1.0;
      }
    }
  }

  const double denom = conc + disc + tied;
  if (denom == 0.0) return NA_REAL;
  return (conc + 0.5 * tied) / denom;
}

// ---------------------------------------------------------------------------
// IPCW competing-risks C-index for a single risk-score vector.
//
// Exact translation of wolbers_c(). G_Ti and G_Ti_minus are precomputed
// once per fold (they do not depend on the risk score) and passed in.
// ---------------------------------------------------------------------------
static double ipcw_cindex(const arma::vec& time,
                          const arma::vec& status,
                          const arma::vec& risk,
                          double evaluation_time,
                          const arma::vec& G_Ti,
                          const arma::vec& G_Ti_minus) {
  const arma::uword n = time.n_elem;
  double numerator = 0.0;
  double denominator = 0.0;

  for (arma::uword i = 0; i < n; ++i) {
    // Subject i must be a cause-1 "case" at or before the evaluation time.
    if (!(time[i] <= evaluation_time && status[i] == 1.0)) continue;

    for (arma::uword j = 0; j < n; ++j) {
      if (i == j) continue;
      if (std::isnan(risk[i]) || std::isnan(risk[j])) continue;

      const bool is_A = time[i] < time[j];
      const bool is_B = (time[i] >= time[j]) && (status[j] == 2.0);
      if (!is_A && !is_B) continue;

      double weight = 0.0;
      if (is_A) {
        const double w_val = G_Ti_minus[i] * G_Ti[i];
        if (w_val > 0.0) weight = 1.0 / w_val;
      } else { // is_B
        const double w_val = G_Ti_minus[i] * G_Ti_minus[j];
        if (w_val > 0.0) weight = 1.0 / w_val;
      }
      if (weight == 0.0) continue;

      double concordant = 0.0;
      if (risk[i] > risk[j]) {
        concordant = 1.0;
      } else if (risk[i] == risk[j]) {
        concordant = 0.5;
      }

      numerator   += weight * concordant;
      denominator += weight;
    }
  }

  if (denominator == 0.0) return NA_REAL;
  return numerator / denominator;
}

//' Predicted risk scores X %*% beta (RcppArmadillo)
//'
//' @param x_test Numeric matrix (n_test x p) of test-fold predictors.
//' @param coef Numeric matrix (p x n_lambda) of coefficients across the
//'   lambda path, e.g. \code{fit$coef} from \code{fastcmprsk::fastCrrp()}.
//' @return Numeric matrix (n_test x n_lambda) of linear-predictor risk scores.
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericMatrix cpp_predict_risk(const Rcpp::NumericMatrix& x_test_,
                                     const Rcpp::NumericMatrix& coef_) {
  const arma::mat x_test = Rcpp::as<arma::mat>(x_test_);
  const arma::mat coef   = Rcpp::as<arma::mat>(coef_);
  arma::mat res = x_test * coef;
  return Rcpp::wrap(res);
}

//' Cross-validation C-index across a lambda path (RcppArmadillo backend)
//'
//' Computes, for each column of \code{coef} (one lambda), the test-fold
//' concordance of the risk scores \code{x_test \%*\% coef[, j]}. This is the
//' compiled kernel behind \code{cv_fastCrrp_cpp()}.
//'
//' @param x_test Numeric matrix (n_test x p) of held-out predictors.
//' @param coef Numeric matrix (p x n_lambda) of fitted coefficients along the
//'   lambda path.
//' @param time_test Numeric vector of held-out event/censoring times.
//' @param status_test Numeric vector of held-out status codes (0 = censored,
//'   1 = cause-1 event, 2 = competing event).
//' @param tuning Either \code{"normal"} (Harrell C-index, matching
//'   \code{survival::concordance}) or \code{"wolbers"} (IPCW competing-risks
//'   C-index, matching \code{wolbers_c}).
//' @param evaluation_time Numeric scalar; the time horizon used only by the
//'   \code{"wolbers"} tuning. Ignored when \code{tuning = "normal"}.
//' @param reverse Logical; orientation of the Harrell C-index for
//'   \code{tuning = "normal"}. \code{FALSE} (default) matches the
//'   \code{survival::concordance} default. Ignored when
//'   \code{tuning = "wolbers"}.
//' @return Numeric vector of length \code{ncol(coef)} of C-index values
//'   (\code{NA} where no comparable pairs exist).
//' @keywords internal
// [[Rcpp::export]]
Rcpp::NumericVector cpp_cv_cindex(const Rcpp::NumericMatrix& x_test_,
                                  const Rcpp::NumericMatrix& coef_,
                                  const Rcpp::NumericVector& time_test_,
                                  const Rcpp::NumericVector& status_test_,
                                  std::string tuning,
                                  double evaluation_time,
                                  bool reverse = false) {
  const arma::mat x_test      = Rcpp::as<arma::mat>(x_test_);
  const arma::mat coef        = Rcpp::as<arma::mat>(coef_);
  const arma::vec time_test   = Rcpp::as<arma::vec>(time_test_);
  const arma::vec status_test = Rcpp::as<arma::vec>(status_test_);
  const arma::uword n_lambda = coef.n_cols;
  arma::vec out(n_lambda);
  out.fill(NA_REAL);

  const bool is_normal  = (tuning == "normal");
  const bool is_wolbers = (tuning == "wolbers");
  if (!is_normal && !is_wolbers) {
    Rcpp::stop("`tuning` must be 'normal' or 'wolbers'.");
  }

  // For the IPCW path, estimate G(t) once per fold (independent of lambda).
  arma::vec G_Ti, G_Ti_minus;
  if (is_wolbers) {
    G_Ti.set_size(time_test.n_elem);
    G_Ti_minus.set_size(time_test.n_elem);
    km_censoring(time_test, status_test, G_Ti, G_Ti_minus);
  }

  // Predictions for all lambda at once, then evaluate column by column.
  arma::mat pred_risk = x_test * coef;  // (n_test x n_lambda)

  for (arma::uword j = 0; j < n_lambda; ++j) {
    arma::vec rs = pred_risk.col(j);
    if (is_normal) {
      out[j] = harrell_cindex(time_test, status_test, rs, reverse);
    } else {
      out[j] = ipcw_cindex(time_test, status_test, rs,
                           evaluation_time, G_Ti, G_Ti_minus);
    }
  }

  // Return a plain (1-D) numeric vector, one C-index per lambda.
  return Rcpp::NumericVector(out.begin(), out.end());
}
