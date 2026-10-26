# bhcrr_autotune() — Bridge Design Note

**Status:** Temporary bridge. `bhcrr_autotune()` is rebuilt on `bhcrr_cv()`
to unblock the `sim/` rebuild. It will be replaced once the simulation layer is
stable. See §5 (WHEN THIS GETS REPLACED).

**Primary caller:** `sim/R/fit.r` (and only that file).

---

## 1. STEPS

A numbered walkthrough of what `bhcrr_autotune()` does end-to-end. All line
references are to `R/autotune.r`.

**1. Input validation** (lines 379–393)

`bhcrr_autotune()` opens with two `inherits()` guards: `mod` must be class
`"ssl_psdh"` (output of `fit_ssl_psdh()`), and `control` must be class
`"bhcrr_cv_control"` (output of `bhcrr_cv_control()`). Failing either stops
immediately with a descriptive message rather than propagating a cryptic
downstream error. The design matrix and outcome are then pulled as `x <- mod$x`
and `y <- mod$y` and checked for existence and conformability.

**2. Pre-flight sizing** (lines 397–420)

`bhcrr_tune_preflight()` (lines 65–179) is called with `mod`, `beta_min`,
`beta_floor`, `fix_s1`, and any caller-supplied `s0_seq`/`s1_seq`. It:
- estimates `theta` (proportion of non-zero coefficients) from
  `mod$coefficients$Estimate`,
- computes the feasible `s1` band via `s1_range(n, p)`,
- solves for the `s0` values that place the zero gap at `beta_min` and at
  `beta_floor` (via `.safe_solve_s0()`),
- assembles recommended `s0_seq` / `s1_seq` sequences that bracket the
  clinical target in zero-gap space.

If the caller supplied sequences they are validated by `bhcrr_check_grid()`;
otherwise the recommendation is used. The feasibility caps for widening are set
here: `s0_ceiling` (tightest `s0_upper_range` cap × 0.98, further clamped to
`fix_s1 × 0.98` when `s1` is coerced), `s1_floor`, `s1_ceiling`.

**3. Grid resolution** (lines 402–420)

After preflight, `s0_seq` and `s1_seq` are finalised: `NULL` inputs are replaced
by the recommendation; a coerced `fix_s1` overwrites `s1_seq` unconditionally.
`s1_fixed` is set to `TRUE` when `fix_s1` is non-NULL; this flag suppresses the
`s1` edge-test in every subsequent `bhcrr_tune_validate()` call and disables
`s1` widening in the loop.

**4. Fold generation** (line 423)

Folds are generated **once**, before the retry loop, by calling
`bhcrr_make_folds(y, control)` when `folds` is `NULL`. Because the folds object
is created here and passed by reference into every `bhcrr_cv()` call, all
attempts in the retry loop use exactly the same train/test split assignments.
A caller who needs exact reproducibility across `bhcrr_autotune()` calls should
supply an explicit `folds` argument.

**5. Attempt loop** (lines 431–491)

The loop is a `repeat {}` with a maximum of `max_widen + 1` iterations.

  **5a. CV call** (lines 434–437): `bhcrr_cv(x, y, s0_seq, s1_seq, control,
  folds, fold_inits)` is called. On the first attempt `fold_inits` is `NULL`;
  `bhcrr_cv()` then runs its initialisation method once per (rep, fold). The
  resulting `$fold_inits` are captured immediately after (step 5b).

  **5b. fold_inits cache** (lines 439–443): after the first CV call, the
  `$fold_inits` list is saved and reused on every subsequent retry. This is
  fold-honest because those inits were each computed on that fold's own
  training rows, and the folds are the same across attempts. Reusing them
  avoids re-running initialisation (LASSO or LASSO_bic), which is typically
  the most expensive part of CV.

  **5c. Validation** (lines 444–449): `bhcrr_tune_validate()` (lines 238–299)
  inspects `cv$tuning`. It checks (i) whether there is a unique maximiser,
  (ii) whether that maximiser is interior to the grid (not on any edge), and
  (iii) whether the zero gap at the optimum lies within `gap_band × beta_min`.
  Results are saved into `history[[attempt]]` alongside `s0_seq`, `s1_seq`,
  `tuning`, and the `cv` object itself.

  **5d. Boundary test and break** (line 450): if the optimum is interior
  (`!val$on_boundary`) or the attempt budget is exhausted (`attempt > max_widen`),
  the loop exits immediately.

  **5e. Widening** (lines 452–468): for each boundary dimension reported by
  `val$boundary_dims`, `.widen_seq()` (lines 304–312) extends the sequence by 2
  points at the sequence's own spacing, in the direction of the boundary,
  clamped by the feasibility caps. `s1` widening is skipped when `s1_fixed` is
  `TRUE` (the guard at line 463 is preserved so the intent survives future
  changes to the boundary test, even though the validate function never reports
  an `s1` boundary when `s1` is coerced).

  **5f. Stall bail-out** (lines 475–491): after widening, `.seq_unchanged()`
  (lines 318–322) checks whether either sequence actually grew. If neither did
  (every candidate point was clamped by the floor or ceiling), `widen_stalled`
  is set to `TRUE`, a flag message is appended, and the loop breaks. The attempt
  just completed stands as the final result; no identical re-fit is attempted.

  If widening did add new points, a flag message is appended (line 492) and the
  loop continues to the next attempt.

**6. Final attempt selection** (lines 493–518)

The last element of `history` is the final attempt. Post-loop flags are
appended for a boundary maximum that survived the full widen budget (line 495)
or one that could not be widened at all (line 497). If the final validation did
not pass all checks, a `warning()` is emitted (line 499). The structured
`"bhcrr_autotune"` list is assembled from the final attempt's fields plus
`$preflight` and a top-level `$cv` pointing to the final `bhcrr_cv` object.

---

## 2. DECISIONS

Each decision is listed with the alternative considered and the reason it was
rejected.

**Control object instead of `nfolds`/`ncv`/`foldid` formals**

The original `bhcrr_autotune()` accepted `nfolds`, `ncv`, and `foldid` as
explicit arguments and threaded them into `tune_ssl_psdh()`. Now that all CV
settings live on a `bhcrr_cv_control` object, duplicating those three fields as
formals would mean two authoritative sources for the same setting (and a silent
precedence question if both are supplied). Collapsing them onto `control` keeps
the interface consistent with `bhcrr_cv()`.

**No `...`**

The original forwarded `...` into `tune_ssl_psdh()`. Dropping it means a call
site that still passes old named arguments (e.g. `nfolds = 5`) fails
immediately with "unused argument" rather than silently doing nothing.
Explicitness is worth the one-time migration cost.

**Data pulled from `mod$x` / `mod$y` rather than new `x`/`y` arguments**

Adding `x` and `y` formals would require the caller to pass the same data twice
(once to fit the initial model, once to autotune). `fit_ssl_psdh()` already
stores both on the returned object, so extracting them from `mod` is the only
source of truth and eliminates the mismatch risk. An error is raised if either
is missing, because an `ssl_psdh` object without `$x`/`$y` is not a supported
state.

**Folds generated once before the retry loop**

Generating new folds inside each attempt would mean comparing CV scores from
different train/test splits — apples to oranges. The optimum could land on a
boundary purely because of fold noise, trigger a widen, and the wider grid could
then appear better or worse purely because it was evaluated on a different split.
Generating folds once makes the score differences across attempts attributable
solely to grid placement.

**`fold_inits` reused across attempts**

Initialisation (LASSO_cv or LASSO_bic) is typically the most expensive part of
CV. Because the folds are identical across attempts, the init vector computed
for fold *k* in attempt 1 was computed on exactly the same training rows as fold
*k* in attempt 2. Reusing it does not bias the CV score comparison; it simply
avoids paying the initialisation cost again for each retry.

**No `tune_ssl_psdh` shim**

A deprecated wrapper returning `bhcrr_cv()$tuning` would have silently broken
any caller that relied on the old field structure (e.g. `$best_s0`, `$best_s1`
vs. the new `$tuning`-based contract). A missing symbol error is more
informative than a silent structural mismatch. The sim rebuild is the only
intended caller, so a clean break is correct.

**`bhcrr_tune_preflight` and `bhcrr_tune_validate` copied unchanged**

Both functions read only `mod$coefficients$Estimate`, `mod$x`, and
`cv$tuning$s0`/`s1`/`score_mean` respectively. Those interfaces are stable
across the `tune_ssl_psdh` → `bhcrr_cv` transition. Copying them unchanged
preserves the exact grid-heuristic logic and validation logic; any change would
be a functional departure from the bridge's stated goal of faithful porting.

**The archive copy left on disk**

`R/Archive/tuning_diagnostics.r` is excluded from the package build via
`.Rbuildignore`. Deleting it would remove the reference implementation. Keeping
it costs nothing and makes the diff between the archive and the bridge
inspectable by `diff` or `git diff`.

**`$cv` added to history and top-level return**

The original history elements did not carry the `bhcrr_cv` object, which made
post-hoc inspection of per-attempt CV details impossible. Adding `$cv` to each
history element and to the top-level return adds no compatibility risk and makes
the bridge self-documenting. It does not change the five stable-surface fields.

---

## 3. I/O CONTRACT

### Arguments

| Argument | Class / type | Shape | Permitted values |
|----------|-------------|-------|-----------------|
| `mod` | `ssl_psdh` (list) | — | Output of `fit_ssl_psdh()`; must carry `$x` (numeric matrix n×p), `$y` (numeric matrix n×2), `$coefficients$Estimate` (numeric vector length p) |
| `beta_min` | numeric | scalar | finite positive; default `0.1` |
| `beta_floor` | numeric | scalar | finite positive, < `beta_min`; default `0.01` |
| `s0_seq` | numeric or `NULL` | vector ≥ 1 | positive values, should satisfy `s0 < s1` for each pair; `NULL` uses preflight recommendation |
| `s1_seq` | numeric or `NULL` | vector ≥ 1 | positive values; `NULL` uses preflight recommendation; overridden by `fix_s1` |
| `fix_s1` | numeric or `NULL` | scalar | single finite positive; `NULL` searches both dimensions |
| `control` | `bhcrr_cv_control` | list | Output of `bhcrr_cv_control()`; governs `nfolds`, `ncv`, init method, parallelism, etc. |
| `folds` | list or `NULL` | output of `bhcrr_make_folds()` | `NULL` generates folds internally once before the retry loop |
| `max_widen` | integer | scalar | non-negative; `0L` disables widening (exactly 1 attempt); default `1L` |
| `gap_band` | numeric | length-2 vector | `c(lo, hi)` with `0 < lo < hi`; multiplicative band around `beta_min` for the zero-gap check; default `c(0.5, 2)` |

### Return fields

The return value is a list of class `"bhcrr_autotune"`.

| Field | Class / type | Shape | Notes |
|-------|-------------|-------|-------|
| `$preflight` | `bhcrr_preflight` | list | from `bhcrr_tune_preflight()`; contains `$n`, `$p`, `$theta`, `$recommended`, `$recommended_check`, etc. |
| `$cv` | `bhcrr_cv` | list | the final `bhcrr_cv()` result (last attempt) |
| `$tuning` | data.frame | n_pairs × 7 | same as `$cv$tuning`; columns `s0`, `s1`, `score_mean`, `score_sd`, `score_fold_sd`, `n_failed_folds`, `n_not_converged` |
| `$validation` | list | — | from `bhcrr_tune_validate()`; contains `$checks`, `$all_pass`, `$best`, `$on_boundary`, `$boundary_dims`, `$gap_at_opt`, `$s1_fixed` |
| `$best` | data.frame or `NULL` | 1 row | one row of `$tuning` at the maximiser; columns identical to `$tuning`; `NULL` if every fit failed |
| `$fix_s1` | numeric or `NULL` | scalar | the coerced `s1` value after `.check_fix_s1()`; `NULL` when `s1` was searched |
| `$s0_seq` | numeric | vector | final `s0` grid used in the last attempt (may be wider than the starting grid) |
| `$s1_seq` | numeric | vector | final `s1` grid used in the last attempt |
| `$history` | list | length = `n_attempts` | each element is a list with `$s0_seq` (numeric vector), `$s1_seq` (numeric vector), `$tuning` (data.frame), `$validation` (list), `$cv` (`bhcrr_cv`); history[[1]] is the first attempt, history[[n_attempts]] is the final attempt |
| `$n_attempts` | integer | scalar | number of CV runs performed (1 when `max_widen = 0`, or when validation passes on the first try, or when widen stalls) |
| `$widen_stalled` | logical | scalar | `TRUE` when the loop ended because widening added no new grid points |
| `$flags` | character | vector, length ≥ 0 | human-readable messages about boundary hits, stalls, and any remaining boundary after the budget was exhausted; empty when everything is clean |
| `$all_pass` | logical | scalar | `TRUE` iff the final `bhcrr_tune_validate()` call passed all checks (no critical failures, no warnings) |

---

## 4. STABLE SURFACE

The `sim/` layer (specifically `sim/R/fit.r`) may depend on **only** these five
fields:

| Field | Type | Meaning |
|-------|------|---------|
| `$best$s0` | numeric scalar | spike scale at the CV maximiser |
| `$best$s1` | numeric scalar | slab scale at the CV maximiser |
| `$all_pass` | logical scalar | `TRUE` iff all validation checks passed |
| `$flags` | character vector | human-readable notes; length 0 when clean |
| `$n_attempts` | integer scalar | how many CV attempts were made |

Everything else — `$preflight`, `$cv`, `$tuning`, `$validation`, `$history`,
`$widen_stalled`, `$fix_s1`, `$s0_seq`, `$s1_seq` — is internal to this bridge
implementation and may change in a future redesign without notice.

A future redesign must preserve these five fields and nothing more.

---

## 5. WHEN THIS GETS REPLACED

### What a real redesign would change

- **Score aggregation.** The current bridge inherits `bhcrr_cv()`'s pooled
  Wolbers C-index as the sole criterion. A redesign might add per-fold
  score summaries, 1-SE rules, or other selectors, and expose them cleanly
  through `control` rather than burying them in post-processing.

- **Grid search strategy.** `.widen_seq()` extends the grid arithmetically at
  the sequence's own spacing (see §6). A redesign would replace this with a
  principled adaptive grid: Bayesian optimisation, a finer initial recommendation,
  or a bracketing strategy that guarantees coverage of the zero-gap target.

- **Separation of preflight and tuning.** `bhcrr_tune_preflight()` and
  `bhcrr_tune_validate()` were designed when the old `tune_ssl_psdh()` call was
  the bottleneck. A redesign can merge the sizing heuristics directly into the
  CV control object or split them into a dedicated `bhcrr_grid()` helper.

- **History format.** Storing a full `bhcrr_cv` object per attempt makes
  `$history` large in memory (especially with `keep_coefs = TRUE`). A redesign
  would store only what is needed for provenance.

### What a redesign must not break

- The five stable-surface fields listed in §4.
- The `"bhcrr_autotune"` class (for `print()` dispatch and `is/inherits` guards
  in `sim/R/fit.r`).
- The argument names `mod`, `control`, `folds`, `max_widen` (callers in
  `sim/R/fit.r` use these by name).

### Intended only caller

`sim/R/fit.r` is the only file that should call `bhcrr_autotune()`. No code in
`analyses/`, `chpc/R/`, or vignettes should reference it (those are being
replaced by the sim rebuild). Before redesigning, audit `sim/R/fit.r` to confirm
which stable-surface fields it actually reads.

---

## 6. CARRIED-OVER LIMITATIONS

The following limitations exist in the bridge and would need to be addressed in
a real redesign. They are documented here so the redesigner does not need to
rediscover them.

**Arithmetic widening cannot reach below the floor.**
`.widen_seq()` adds points by stepping outward at the sequence's own spacing.
When the sequence has only one point, the step is `max(s) × 0.1`. When the
sequence already starts near the floor (`1e-4`), the new low points fall at or
below the floor and are clamped away, producing `widen_stalled = TRUE`. A
low-edge optimum near the floor will always stall and the returned `$best$s0`
will be the smallest feasible value the grid could reach, not the true optimum.

**Score metric: higher is better.**
The Wolbers C-index is used as the CV score. It is bounded in (0, 1) and higher
values indicate better discriminative ability. The maximiser of `score_mean` is
selected as `$best`. There is no 1-SE rule or complexity penalty.

**Evaluation horizon resolution.**
`bhcrr_cv()` resolves the evaluation horizon `τ` from `control$eval_time` when
supplied, or falls back to a quantile of the full-data cause-1 event times via
`control$eval_quantile` (default 0.5). The fallback uses all of `y`, which is a
marginal, model-free quantity, but it still consumes information from outside any
given fold. For results going into the manuscript, **always set
`control$eval_time` explicitly** so `τ` is a fixed, pre-specified value and not
a data-derived one.

**`theta` is estimated once from the initial model.**
The `theta` passed to `bhcrr_tune_validate()` and the zero-gap heuristics is
estimated from `mod$coefficients$Estimate` — the initial, untuned fit. If the
tuned model would have a very different sparsity, the gap bounds may be
miscalibrated. This is inherent to using a pre-tuning fit for grid sizing.

**`fix_s1` disables s1 widening unconditionally.**
When `fix_s1` is supplied, `s1_seq` is a single point and the validate function
never reports an `s1` boundary. The guard at line 463 ensures `s1` is not
widened even if this logic changes. This is by design: a coerced `s1` is a
deliberate act, not a grid-search outcome.
