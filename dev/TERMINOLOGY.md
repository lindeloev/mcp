# mcp terminology
Use one word per concept, so users are not confused:

- Documentation, messages, user-facing output, and code comments use the full terms below.
- Code may use the listed abbreviations (`pars`, `coefs`, `cp`), but nothing else.


## Formulas and segments
`vignettes/formulas.Rmd#how-segments-connect` defines these:

- **included** / **not included** (terms): what terms are written in a formula. This matches base R wording. Prefer "not included" over "omitted", "dropped", or "excluded". Avoid: declared, active, written.
- **x-terms**: formula terms involving `par_x`, i.e., the change point dimension, e.g., `x`, `I(x^2)`, `x:z`. They are *measured from the change point*. Avoid: segment-local, "count from".
- **joined**: a segment without an intercept (`0 + ...`). Avoid: carry, persist, retain.
- **disjoined**: a segment with an intercept. It starts afresh. Applies separately to each dpar, each grouping column of its group-level terms, and each AR/MA lag, e.g., "joined mu, disjoined sigma". As a verb, the intercept *disjoins* the earlier joined terms. Avoid: reset.
- **turns off**: `(0 | group)`, `ar(0)`/`ma(0)`, or not including an earlier AR/MA lag. The earlier group-level terms or AR/MA terms no longer apply from that segment, i.e., they are disjoined from zero. Use only for these; use "disjoined" for intercepts, and "not included" for other terms.
- **linear predictor**: the segmented linear predictor η of a dpar. Always write "linear predictor" in full.
- **change-point part** / **predictor part**: the two parts of a segment formula, i.e., `cp()` and the rest. Matches the `part` column of `mcp_pars()`.
- **`same()`**: include an earlier coefficient in a segment. Name `same()` where it applies, e.g., "terms in `same()`" or "`same(x)` includes `x_1` in segment 2". "Reuse" is OK only as the verb for what `same()` does, e.g., "`same(x)` reuses `x_1`". Avoid "share" except for priors that share a value, e.g., `Intercept_2 = "Intercept_1"`.
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


## Data
- **variable**: a column in the data, e.g., the change-point variable `par_x`.
- **predictor**: a variable included in a formula, e.g., "continuous predictor" or "categorical predictor". One predictor can give several coefficients, e.g., a factor or `x + I(x^2)`, so do not use "predictor" for a coefficient.
- **level**: a factor level, or in "group-level". Avoid for anything else!
- **par_x** should consistently be called "the change-point variable"


## Priors
See `vignettes/priors.Rmd#cp_prior`.

- **specified prior**: the prior as set by the user or by default, i.e., what `fit$prior` and `prior_summary()` show.
- **marginal prior**: the prior of one parameter after all restrictions, i.e., what `plot_pars(fit, prior = TRUE)` shows. OK to write "resulting marginal prior" when explaining effect of restrictions.
- **constant** / **constant prior**: a numerical prior, e.g., `cp_2 = 80`, as in `brms`'s `constant()` prior and `kind = "constant"` in `prior_summary(fit, verbose = TRUE)`. The parameter is not estimated. "Fix" is OK as a verb, e.g., "fix `cp_2` to 80". Avoid: fixed value, fixed parameter, since "fixed" can be confused with fixed (population-level) effects and `fixef()`.
- **data-dependent terms**: `min(x)`, `sd(y)`, `n_cp()`, etc. in prior strings.
- **resolved**: with data-dependent terms replaced by numbers.
- **bounds** / **boundary**: the lower and upper bounds of a parameter: the arguments of `dunif()` or `T()`, or the support of its prior distribution (e.g., `[0, 1]` for `dbeta()`, `[0, Inf]` for `sigma`). Avoid: limits, range, constraints, support.
- **restriction** / **restricted to**: anything that removes prior mass, e.g., ordering.
- **truncation**: JAGS's renormalizing `T()`. For change points, prefer "restriction", because the ordering and minimum segment lengths are not renormalizing.


## Draws and predictions
- **draws**: prior or posterior draws, following the `posterior` package. Avoid: samples (except for "sample" as a verb, e.g., `mcp(..., sample = "both")`).
- **empty fit**: an `mcpfit` without draws, from `mcp(..., sample = FALSE)`, e.g., to get `fit$simulate()`, `fit$prior`, or `fit$jags_code`. Conventionally named `empty` in examples. Avoid: empty model.


## Words reserved for other concepts
Do not use "bounds" for these.

- **observed range** (of x): the data from `min(x)` to `max(x)`. Its **ends** are `min(x)` and `max(x)`.
- **lower and upper ends**: of a central interval, e.g., `.lower`/`.upper` and `lower`/`upper`.
- **limits**: of plots, e.g., `limits_y` (ggplot convention).


## Code only
- **definition / occurrence**: a coefficient is defined in the first segment that includes it (`definition_segment`) and occurs there and in every later segment with it in `same()`. In docs, refer to `same()` instead.
- **par / pars**: Shorthand for "parameter".
- **coef / coefs**: Shorthand for "coefficient".
- **cp**: Shorthand for "change point".
- **cols / \*_cols**: Shorthand for a character vector of column names, e.g., `group_cols` and `required_cols`.
- **data_columns**: The named list of column roles returned by `mcp_columns()` (`par_x`, `response`, `series`, `trials`, ...). Not shortened, to match `mcp_columns()` and to keep it distinct from `*_cols`.
- **x_range / range_ends / range_nodes**: The observed range of x as numbers (`range(x)`), as strings (`"min(x)"`, `"max(x)"`), and as the JAGS nodes `cp_0` and `cp_K`.
- **name**: In model tables, the parameter name (e.g., `cp_1`, `x_1`); `mcp_pars()` shows it as parameter.
