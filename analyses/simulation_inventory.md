# Simulation Vignette Inventory
**Source:** `vignettes/Simulations.Rmd`  
**Date inventoried:** 2026-04-17

---

## Packages loaded (setup chunk, lines 25–31)

| Package | Role |
|---------|------|
| `dplyr` | Data manipulation (tidyverse) |
| `ggplot2` | Plotting (tidyverse) |
| `simsurv` | Parametric survival time simulation |
| `patchwork` | Composing ggplot panels |
| `survival` | Cox models, `Surv()`, `survfit()` (loaded later, line 930) |
| `fastcmprsk` | Fine-Gray / PSDH competing risks (loaded line 1621, with explicit `lib.loc`) |
| `BhGLM` | Bayesian spike-and-slab GLM (`bmlasso`, `cv.bh`) (loaded line 1622, with explicit `lib.loc`) |

The vignette also `source()`s the local bhCRR R files (lines 1626–1636):  
`fit_ssl_psdh`, `cv_ssl_psdh`, `update_betas`, `expected_inclusion_probs`,  
`expected_penalty_weights`, `update_mixture_prob`, `helpers`, `tune_ssl_psdh`,  
`measure_ssl_psdh`, `predict_from_ssl_psdh`, `generate_foldid`.

---

## Functions defined in the vignette

### Utility / infrastructure

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `` `%||%` `` | 73 | `(a, b)` | base only | No |
| `sim_merge_lists` | 75 | `(default, user)` | base only | No |

---

### Simulation specification

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `sim_spec` | 100 | `(n, p, x, truth, risks, censor, seed)` | base, calls `sim_merge_lists` | No |

Returns a `sim_spec` S3 object (list with class attribute).

---

### Simulation generators

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `sim_generate` | 186 | `(spec, custom_design_matrix, custom_predictor_set, custom_betas, seed)` | base; calls `sim_generate_x`, `sim_generate_truth`, `sim_compute_effects`, `sim_weibull_ph`, `sim_generate_censor_time` | No |
| `sim_generate_x` | 370 | `(xspec, n, p)` | `stats` (`rnorm`, `qnorm`, `rpois`, `runif`, `sd`); base (`switch`, `matrix`, `apply`, `sweep`, `colMeans`) | No |
| `sim_generate_truth` | 477 | `(trspec, p, seed)` | **`dplyr::intersect`** (line 498); base (`sample`, `ceiling`) | **Yes — `dplyr::intersect`** |
| `sim_compute_effects` | 518 | `(X, truth)` | base (`scale`, `rnorm`, `sd`, `mean`, `optimize`, `match`, `numeric`) | No |

`sim_compute_effects` defines one nested anonymous function:

- **`objective`** (line 577) — nested inside `sim_compute_effects`; base only.

---

### Outcome / censoring generators

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `sim_weibull_ph` | 742 | `(shape, scale, eta_val)` | **`simsurv::simsurv`** | No |
| `sim_generate_censor_time` | 760 | `(cspec, n, maxtime)` | base (`rexp`, `rep`, `switch`) | No |

---

### Data loading helper

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `permute_rows` | 797 | `(df, seed)` | base (`as.matrix`, `t`, `apply`, `sample.int`, `colnames`, `rownames`, `as.data.frame`) | No |

---

### Diagnostics

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `outcome_simulation_diagnostic_plots` | 931 | `(sim, shapes, scales)` | **`ggplot2`**, **`survival`** (`coxph`, `survfit`, `Surv`, `residuals`), **`patchwork`**, **`dplyr`**, **`tidyr`** | **Yes — heavy** (see detail below) |

`outcome_simulation_diagnostic_plots` defines one nested function:

- **`cox_snell_func`** (line 1011) — nested inside `outcome_simulation_diagnostic_plots`.  
  Deps: **`survival`** (`coxph`, `survfit`, `Surv`), **`ggplot2`**, **`dplyr`** (`%>%`, `select`).  
  Tidyverse: `%>%` pipe (magrittr/dplyr), bare `select()`.

---

### Modeling wrappers

| Function | Line | Signature | Package deps | Tidyverse? |
|----------|------|-----------|--------------|------------|
| `cs_oracle_func` | 1672 | `(df)` | **`survival`** (`coxph`, `Surv`, `summary`), **`dplyr`** (`filter`, `mutate`, `select`) | **Yes — `dplyr::mutate`, `dplyr::select`, bare `filter`** |
| `psdh_oracle_func` | 1710 | `(x, y_time, y_status, failcode_num)` | **`fastcmprsk`** (`fastCrr`, `Crisk`), **`dplyr`** (`filter`, `mutate`, `select`) | **Yes — `dplyr::mutate`, `dplyr::select`, bare `filter`** |
| `cs_ssl_func` | 1749 | `(x, y)` | **`BhGLM`** (`bmlasso`, `cv.bh`), base (`rbind`, `data.frame`, `sample`) | No |
| `psdh_ssl_func` | 1830 | `(x, y)` | bhCRR sourced fns (`fit_ssl_psdh`, `tune_ssl_psdh`), **`dplyr`** (`%>%`, `filter`) | **Yes — `%>%` pipe, bare `filter`** |
| `model_sim_func` | 1877 | `(sim)` | **`survival`** (`Surv`), **`dplyr`** (`mutate`, `filter`, `select`, `full_join`, `case_when`), calls all modeling wrappers above | **Yes — heavy** (see detail below) |

---

## Tidyverse-specific code — complete catalogue

### `%>%` magrittr pipes (from `dplyr`)

Appears in: `cox_snell_func` (line 1012), `outcome_simulation_diagnostic_plots` table block (lines 1052–1072), `psdh_ssl_func` (line 1838), and all six summary-table / plot code blocks after each scenario (lines 1157–1615).

> **Note:** `cs_oracle_func`, `psdh_oracle_func`, and `model_sim_func` use the native `|>` pipe instead of `%>%`.

### `dplyr` verbs

| Verb | Locations |
|------|-----------|
| `select` | `cox_snell_func` (line 1013); `outcome_simulation_diagnostic_plots` table (line 1062); scenario summaries (lines 1165, 1257, 1349, 1439, 1519, 1601); `model_sim_func` (lines 1951, 1969) |
| `mutate` | `outcome_simulation_diagnostic_plots` table (lines 1053, 1061, 1064); `cs_oracle_func` (line 1687); `psdh_oracle_func` (line 1725); `model_sim_func` (lines 1949, 1967) |
| `filter` | `cs_oracle_func` (line 1686); `psdh_oracle_func` (line 1724); `psdh_ssl_func` (line 1839); `model_sim_func` (lines 1958, 1963, 1976, 1979); scenario summary plots (line 1170) |
| `group_by` | `outcome_simulation_diagnostic_plots` table (line 1056); scenario summary chunks (lines 1158, 1248, 1341, 1432, 1512, 1596) |
| `summarize` / `summarise` | `outcome_simulation_diagnostic_plots` table (line 1057); scenario summary chunks |
| `ungroup` | `outcome_simulation_diagnostic_plots` table (line 1063); scenario summary chunks (lines 1169, 1261, 1353, etc.) |
| `full_join` | `model_sim_func` (lines 1952, 1954, 1956, 1960, 1970, 1972, 1974, 1978) |
| `case_when` | `outcome_simulation_diagnostic_plots` table (lines 1064, 1067, 1070); `model_sim_func` (lines 1949, 1967) |
| `dplyr::intersect` | `sim_generate_truth` (line 498) — **in an otherwise base-only simulation function** |
| `dplyr::starts_with` | Data loading block (line 806) — used inside `dplyr::select()` |

### `tidyr` verbs

| Verb | Locations |
|------|-----------|
| `tidyr::pivot_longer` | Scenario plot blocks: lines 1166, 1257, 1349, 1439, 1519, 1602 — always called with explicit `tidyr::` namespace |

### `ggplot2` (tidyverse-adjacent, not a data-manipulation concern)

Used throughout `outcome_simulation_diagnostic_plots` and all scenario plot blocks. Main verbs: `ggplot`, `aes`, `geom_density`, `geom_step`, `geom_boxplot`, `geom_abline`, `stat_function`, `scale_color_manual`, `labs`, `theme_classic`, `ggtitle`. All ggplot2 calls use `+` composition, not pipes.

---

## Summary of tidyverse exposure by function

| Function | Tidyverse exposure |
|----------|--------------------|
| `` `%||%` `` | None |
| `sim_merge_lists` | None |
| `sim_spec` | None |
| `sim_generate` | None |
| `sim_generate_x` | None |
| `sim_generate_truth` | `dplyr::intersect` (one call) |
| `sim_compute_effects` | None |
| `sim_weibull_ph` | None |
| `sim_generate_censor_time` | None |
| `permute_rows` | None |
| `outcome_simulation_diagnostic_plots` | Heavy: `%>%`, `select`, `mutate`, `group_by`, `summarize`, `ungroup`, `case_when`, `tidyr::pivot_longer`, `ggplot2` |
| ↳ `cox_snell_func` (nested) | `%>%`, bare `select`; `ggplot2` |
| `cs_oracle_func` | `dplyr::mutate`, `dplyr::select`, bare `filter`; native `\|>` pipe |
| `psdh_oracle_func` | `dplyr::mutate`, `dplyr::select`, bare `filter`; native `\|>` pipe |
| `cs_ssl_func` | None |
| `psdh_ssl_func` | `%>%`, bare `filter` |
| `model_sim_func` | `dplyr::mutate`, `dplyr::filter`, bare `select`, `full_join`, `case_when`; native `\|>` pipe |

---

## Potential refactoring notes (not implemented here)

- `sim_generate_truth` uses `dplyr::intersect` but `base::intersect` is identical in behavior — easy to swap.
- The `%>%` pipe in `psdh_ssl_func` and `cox_snell_func` (nested) could be replaced with native `|>` for consistency with `cs_oracle_func` / `psdh_oracle_func` / `model_sim_func`.
- `dplyr` verbs in `outcome_simulation_diagnostic_plots`'s table block are the most complex tidyverse dependency in the simulation infrastructure; they could be replaced with base R if the diagnostic function ever needs to be vignette-independent.