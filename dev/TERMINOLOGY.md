# mcp terminology
Use one word per concept in documentation, messages, and user-facing output, so users are not confused.


## Formulas and segments
`vignettes/formulas.Rmd#how-segments-connect` defines these:

- **included** / **not included** (terms): what terms are written in a formula. This matches base R wording. Prefer "not included" over "omitted", "dropped", or "excluded". Avoid: declared, active, written.
- **x-terms**: formula terms involving `par_x`, i.e., the change point dimension, e.g., `x`, `I(x^2)`, `x:z`. They are *measured from the change point*. Avoid: segment-local, "count from".
- **joined**: a segment without an intercept (`0 + ...`). Avoid: carry, persist, retain.
- **disjoined**: a segment with an intercept. It starts afresh. Avoid: reset.
- **predictor**: a variable included in a formula, e.g., "continuous predictor" or "categorical predictor".
- **linear predictor**: the segmented linear predictor η of a dpar. Always write "linear predictor" in full.
- **change-point part** / **predictor part**: the two parts of a segment formula, i.e., `cp()` and the rest. Matches the `part` column of `mcp_pars()`.
- **reuse** and **same**: include an earlier coefficient in a segment with `same()`, e.g., "reuse a coefficient". "same" is an OK word too. Avoid "share" except for priors that share a value, e.g., `Intercept_2 = "Intercept_1"`.
- **plateau**: a segment or dpar without x-terms, e.g., `~ 1` or a joined `~ 0`.
- **slope**: the coefficient of `x`. Informal; use "x-terms" for the general case.


## Parameters
`mcp_pars(fit)` classifies every parameter by `role` and `scope`. Its column names and values define these terms.

- **parameter**: anything with a prior and posterior, i.e., everything `mcp_pars(fit)` lists: change points, coefficients, group-level SDs, and group-level deviations.
- **coefficient**: the coefficient of an included term in a dpar's linear predictor (`role = "coefficient"`). This includes intercepts and group-level deviations.
- **change point**: `cp_1`, `cp_2`, ... (`role = "change_point"`). Not a dpar: its `dpar` is `NA`.
- **distributional parameter** (`dpar`): anything with its own segmented linear predictor: `mu`, `sigma`, `shape`, and each AR/MA lag (`ar1`, `ma1`, ...). Defined in `vignettes/dpar.Rmd`. `sigma` is a dpar, whereas `sigma_1` is a coefficient.
- **group-level effect**: a group-level term in a formula, e.g., `(1 | id)` or `cp((1 | id))`. Avoid: group effect.
- **group-level deviation**: a parameter of a group-level effect, e.g., `Intercept_1_id` or `cp_1_id`.
- **group-specific**: the resulting value for each group, e.g., "group-specific change point" for `cp_1 + cp_1_id[g]`.
- **population-level change point**: e.g., `cp_1` in a model with group-level change points. Avoid: population change point.
- **realized**: the result of the bounds of group-level change points, e.g., "realized SD" of the group-specific change points, as opposed to `cp_1_id_sd`, which is the SD *before truncation*.
- **population-level** / **group-level**: the `scope` of a parameter. Population-level includes the coefficients of all dpars, including AR/MA. It is OK to paranthetically mention fixed/random effects for lme4 users.
- **variable**: a column in the data, e.g., the change-point variable `par_x`.
- **level**: a factor level, or in "group-level". Avoid for anything else!
- **par_x** should consistently be called "the change-point variable"

## Priors
See `vignettes/priors.Rmd#cp_prior`.

- **specified prior**: the prior as set by the user or by default, i.e., what `fit$prior` and `prior_summary()` show.
- **marginal prior**: the prior of one parameter after all restrictions, i.e., what `plot_pars(fit, prior = TRUE)` shows. OK to write "resulting marginal prior" when explaining effect of restrictions.
- **data-dependent terms**: `min(x)`, `sd(y)`, `n_cp()`, etc. in prior strings.
- **resolved**: with data-dependent terms replaced by numbers.
- **bounds** / **boundary**: the lower and upper bounds of a parameter: the arguments of `dunif()` or `T()`, or the support of its prior distribution (e.g., `[0, 1]` for `dbeta()`, `[0, Inf]` for `sigma`). Avoid: limits, range, constraints, support.
- **restriction** / **restricted to**: anything that removes prior mass, e.g., ordering.
- **truncation**: JAGS's renormalizing `T()`. For change points, prefer "restriction", because the ordering and minimum segment lengths are not renormalizing.


## Draws and predictions
- **draws**: prior or posterior draws, following the `posterior` package. Avoid: samples (except for "sample" as a verb, e.g., `mcp(..., sample = "both")`).


## Words reserved for other concepts
Do not use "bounds" for these.

- **observed range** (of x): the data from `min(x)` to `max(x)`. Its **ends** are `min(x)` and `max(x)`.
- **lower and upper ends**: of a central interval, e.g., `.lower`/`.upper` and `lower`/`upper`.
- **limits**: of plots, e.g., `limits_y` (ggplot convention).
