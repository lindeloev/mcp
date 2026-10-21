# Compare cp_1 mixing for the Bernoulli fit-recovery models across seeds. Run under JAGS 4 and JAGS 5.
pkgload::load_all(quiet = TRUE)
model = list(y ~ 1 + x, ~ 1 + x)
designs = list(
  logit = list(family = bernoulli(link = "logit"), chains = 2, simulated = list(cp_1 = 100, Intercept_1 = -0.4, x_1 = 0.006, Intercept_2 = 2.2, x_2 = -0.004)),
  probit = list(family = bernoulli(link = "probit"), chains = 2, simulated = list(cp_1 = 93, Intercept_1 = -0.3, x_1 = 0.005, Intercept_2 = 1.4, x_2 = -0.003)),
  identity = list(family = bernoulli(link = "identity"), chains = 3, simulated = list(cp_1 = 93, Intercept_1 = 0.2, x_1 = 0.001, Intercept_2 = 0.8, x_2 = -0.001))
)

results = do.call(rbind, lapply(names(designs), function(link) {
  d = designs[[link]]
  newdata = data.frame(x = seq(1, 200, length.out = 400), y = 0)
  empty = mcp(model, data = newdata, family = d$family, sample = FALSE, par_x = "x")
  set.seed(42)
  newdata$y = suppressMessages(do.call(empty$simulate, c(list(fit = empty, newdata = newdata), d$simulated)))
  rows = parallel::mclapply(1:12, function(seed) {
    fit = mcp(model, newdata, family = d$family, par_x = "x", chains = d$chains, warmup = 500, iter = 1000, seed = seed, diagnostics = FALSE, quiet = TRUE)
    s = get_summary(fit, width = 0.97, scope = "population")
    data.frame(link = link, seed = seed, cp_rhat = s$rhat[1], cp_ess = s$ess_bulk[1], min_ess = min(s$ess_bulk))
  }, mc.cores = max(1, parallel::detectCores() - 1))
  do.call(rbind, rows)
}))

cat("JAGS", format(tryCatch(rjags::jags.version(), error = function(e) conditionMessage(e))), "\n")
print(results, row.names = FALSE)
print(aggregate(cbind(cp_rhat, cp_ess, min_ess) ~ link, results, median))
cat("Fraction with min_ess <= 50:\n")
print(tapply(results$min_ess <= 50, results$link, mean))
