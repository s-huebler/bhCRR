# Algorithm–Code Synthesis: bhCRR

**Date:** 2026-04-07  
**Files compared:** `manuscript/main.tex` (Algorithms 1–3 and their surrounding math) vs. all `.r` files in `R/`

---

## Overview

The manuscript contains three formal algorithms (EM with CCD, Pre-validation, and the IPCW C-index), plus supporting mathematical derivations for the E-step quantities. The table below summarizes alignment status before the detailed notes.

| Component | Manuscript Section | Code File(s) | Alignment |
|---|---|---|---|
| EM loop structure (E/M order) | §2.3, Algorithm 1 | `fit_ssl_psdh.r` | ✅ Matches |
| Inclusion probability (E-step) | Eq for p_j^(k) | `expected_inclusion_probs.r` | ⚠️ **Bug: slab mean = 1, not 0** |
| Penalty weights (E-step) | λ_j^(k) formula | `expected_penalty_weights.r` | ✅ Matches exactly |
| Mixture proportion update (M-step) | θ^(k+1) = mean(p_j) | `fit_ssl_psdh.r`, `update_mixture_prob.r` | ✅ Matches |
| Beta update (M-step) | fastCrrp CCD | `update_betas.r` | ✅ Matches |
| Convergence criterion | Algorithm 1 | `fit_ssl_psdh.r` | ⚠️ **Denominator sign discrepancy** |
| Initialization of penalty weights | §2.3 (init) | `fit_ssl_psdh.r` | ⚠️ **Dimension mismatch on init** |
| Pre-validation loop | Algorithm 2 | `cv_ssl_psdh.r` | ✅ Logic matches; note on output below |
| Pre-validation output | Algorithm 2: η = Xβ | `predict_from_ssl_psdh.r` | ⚠️ **Returns CIF, not linear predictor** |
| C-index: case and pair definition | Algorithm 3 | `measure_ssl_psdh.r` | ✅ Matches |
| C-index: IPCW weights | Algorithm 3 | `measure_ssl_psdh.r` | ✅ Matches |
| C-index: concordance + tie handling | Algorithm 3 | `measure_ssl_psdh.r` | ✅ Matches |
| C-index: evaluation time | §2.3.4 | `cv_ssl_psdh.r` | ⚠️ **Note: includes competing event times** |

---

## Detailed Findings

### 1. EM Algorithm (Algorithm 1) — `fit_ssl_psdh.r`

**Structure:** The E-step → M-step ordering in the code matches the manuscript exactly. Within the E-step, inclusion probabilities are computed first (`expected_inclusion_probs`), then penalty weights (`expected_penalty_weights`). In the M-step, the mixture proportion (`mean(inclusion_probs)`) is updated first, then betas via `update_betas`. This mirrors Algorithm 1 faithfully.

---

### 2. ⚠️ Bug: Slab distribution mean — `expected_inclusion_probs.r`

**Manuscript** (§2.2, hierarchy specification):
```
β_j | γ_j, s₀, s₁ ~ (1 - γ_j) DE(β_j | 0, s₀) + γ_j DE(β_j | 0, s₁)
```
Both spike and slab are **centered at 0**.

**Code** (`expected_inclusion_probs.r`, line 16):
```r
dens_Slab <- dlaplace(betas, mu = 1, b = s1)   # <-- mu = 1
dens_Spike <- dlaplace(betas, mu = 0, b = s0)  # correct
```

The slab density is evaluated with `mu = 1` rather than `mu = 0`. This means the inclusion probability calculation biases coefficients near 1 toward the slab (high inclusion), and coefficients near 0 would be pulled toward the spike more than the model intends. The closed-form inclusion probability derivation in §2.3.2 of the manuscript is derived assuming both densities are centered at 0. Setting `mu = 1` for the slab is inconsistent with all the prior specification math.

**Fix:** Change `mu = 1` to `mu = 0` in `expected_inclusion_probs.r`.

---

### 3. ⚠️ Convergence Criterion Denominator — `fit_ssl_psdh.r`

**Manuscript** (Algorithm 1, stopping rule):
```
|d^(t) - d^(t-1)| / (0.1 - |d^(t)|) < ε
```

**Code** (`fit_ssl_psdh.r`, line 80):
```r
if(abs(logLik - devold)/(0.1 + abs(logLik)) < epsilon & iter > 5)
```

The denominator uses `+` where the manuscript has `-`. The manuscript itself flags this in red italics: *"This stopping criteria seems incorrectly implemented."*

For context: in the Tang (2017) formulation this is adapted from, the stopping rule is typically `|d^(t) - d^(t-1)| / (0.1 + |d^(t)|) < ε`, which is what the code actually implements. Since log-likelihoods are negative in survival models, using `0.1 - |d|` in the denominator would yield a negative denominator (which flips the inequality), making it unusable. The **code is likely correct** and the manuscript has the wrong sign. The manuscript needs to be corrected to `0.1 + |d^(t)|`.

---

### 4. ⚠️ Penalty Weight Initialization Dimension — `fit_ssl_psdh.r`

**Manuscript:** Initialization is described abstractly (initial β₀, θ₀). The natural reading is that θ is a scalar and penalty weights λ_j are length-p (one per feature).

**Code** (`fit_ssl_psdh.r`, lines 25–28):
```r
current_mixture_prob <- rep(initial_sparsity, nrow(x))  # length n (observations!)
init_mixture_scale   <- (1-current_mixture_prob)*ss0 + current_mixture_prob*ss1
current_penalty_weights <- 1/init_mixture_scale          # also length n
```

`nrow(x)` is the number of **observations** (n), but penalty weights should be length **p** (number of features, `ncol(x)`). The initial penalty weights vector is length-n instead of length-p. This is passed to `fastCrrp` as `penalty.factor`, which expects a length-p vector. R's recycling will silently handle mismatched lengths in some cases, but the behavior is undefined and likely wrong. After the first E-step, `expected_penalty_weights` correctly produces length-p weights (from length-p betas), so this only affects initialization.

**Fix:** Change `nrow(x)` to `ncol(x)` in the initialization block:
```r
current_mixture_prob    <- rep(initial_sparsity, ncol(x))  # was nrow(x)
init_mixture_scale      <- (1 - current_mixture_prob)*ss0 + current_mixture_prob*ss1
current_penalty_weights <- 1/init_mixture_scale
```

---

### 5. ✅ Penalty Weight Formula — `expected_penalty_weights.r`

**Manuscript** (§2.3.2):
```
λ_j^(k) = (1 - p_j^(k)) / s₀ + p_j^(k) / s₁
```

**Code:**
```r
expected_penalty_weights <- function(s1, s0, p){
  (1 - p) / s0 + p / s1
}
```
Exact match.

---

### 6. ✅ Mixture Proportion Update — `fit_ssl_psdh.r`

**Manuscript** (M-step):
```
θ^(k+1) = (1/J) Σ p_j^(k)
```

**Code:**
```r
current_mixture_prob <- mean(current_inclusion_probs)
```
Exact match (`mean` = sum/J).

---

### 7. ✅ Beta Update — `update_betas.r`

**Manuscript:** Maximize Q₁ via CCD for a LASSO-penalized PSH model with per-coefficient weights λ_j^(k), using the fastcmprsk implementation.

**Code:**
```r
fastcmprsk::fastCrrp(..., penalty.factor = penalty_weights, penalty = "LASSO", lambda = lambda)
```
Matches exactly. The global `lambda` is computed in `fit_ssl_psdh.r` as `sum(Pf) / (nrow(x) * ncol(x))`, which scales the per-feature penalty factors into a valid global multiplier.

---

### 8. ⚠️ Pre-validation Output: CIF vs. Linear Predictor — `predict_from_ssl_psdh.r`

**Manuscript** (Algorithm 2, Pre-validation):
```
η̂ᵢ = Xᵢ β̂^(-k)   (linear predictor)
```

**Code** (`predict_from_ssl_psdh.r`):
```r
lp <- as.vector(newx %*% beta)
# ... computes baseline cumulative hazard ...
predicted_risk <- 1 - exp(-cum_base_haz * exp(lp))   # CIF, not linear predictor
```

The function returns the absolute risk (CIF) evaluated at `prediction_time`, not the linear predictor. For the purpose of **ranking** (which is all that the c-index requires), this is monotonically equivalent to the linear predictor when the baseline hazard is fixed — so the c-index values will be correct. However:

- The manuscript describes this as computing a prognostic *index*, which is the linear predictor.
- If CVPL is ever implemented using these predictions, the CIF cannot be substituted for the linear predictor.
- The manuscript should clarify that the cross-validated score used is the risk at the evaluation time, or the code should optionally expose the linear predictor.

---

### 9. ✅ C-index Algorithm (Algorithm 3) — `measure_ssl_psdh.r`

The IPCW C-index implementation is a close match to the manuscript's Algorithm 3.

**Censoring distribution:** KM fit with inverted status (correct).

**Case definition:** Subject i is a case if `T_i ≤ t` and `D_i = 1` (event of interest). Code:
```r
if (!(T_i[i] <= evaluation_time && D_i[i] == 1)) next
```
Matches.

**Pair evaluability:** 
- `A_ij = I(T_i < T_j)` — code: `is_A <- T_i[i] < T_i[j]` ✅
- `B_ij = I(T_i ≥ T_j and D_j = 2)` — code: `is_B <- (T_i[i] >= T_i[j]) && (D_i[j] == 2)` ✅

**IPCW weights:**
- If A: `W_ij = [G(T_i⁻) × G(T_i)]⁻¹` — code uses `G_Ti_minus[i] * G_Ti[i]` ✅
- If B: `W_ij = [G(T_i⁻) × G(T_j⁻)]⁻¹` — code uses `G_Ti_minus[i] * G_Ti_minus[j]` ✅

**Concordance and tie handling:** Correct (ties get 0.5, higher risk for case = concordant). ✅

---

### 10. ⚠️ C-index Evaluation Time Composition — `cv_ssl_psdh.r`

**Manuscript** (§2.3.4, note in red):
> "the marginal distribution of event times is comprised of times corresponding to both the event of interest as well as the competing event, but not the censoring times. *Need to change in code to only evaluate event times.*"

**Code** (`cv_ssl_psdh.r`, line 21):
```r
eval_time <- as.numeric(quantile(y[,1], eval_quantile))
```

`y[,1]` is all observation times, including censored observations. The manuscript already flags that this should use only non-censored event times (or even only event-of-interest times). This is a known TODO in both the manuscript and code.

---

## Summary of Action Items

| Priority | Issue | File | Fix |
|---|---|---|---|
| 🔴 High | Slab density mean = 1 (should be 0) | `expected_inclusion_probs.r` | Change `mu = 1` → `mu = 0` |
| 🔴 High | Init penalty weights length n (should be p) | `fit_ssl_psdh.r` | Change `nrow(x)` → `ncol(x)` |
| 🟡 Medium | Convergence denominator sign in manuscript | `manuscript/main.tex` | Change `0.1 - |d|` → `0.1 + |d|` |
| 🟡 Medium | Eval time uses all obs times, not event-only | `cv_ssl_psdh.r` | Filter `y[y[,2] != 0, 1]` for eval time |
| 🟢 Low | Pre-validation returns CIF not linear predictor | `predict_from_ssl_psdh.r` | Clarify in manuscript; expose lp option in function |
