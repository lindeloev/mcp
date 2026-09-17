test_that("shared group blocks define one hierarchy and simulate one deviation per group", {
  data = data.frame(x = 1:12, y = 0, id = rep(c("a", "b"), 6))
  fit = mcp(list(y ~ 1 + (1 | id), ~ 0 + same((1 || id))), data, par_x = "x", sample = FALSE)
  expect_equal(mcp_pars(fit, scope = "group")$name, "Intercept_1_id")
  expect_equal(mcp_pars(fit, role = "group_sd")$name, "Intercept_1_id_sd")
  for (ids in list(data$id, rep(c("new-a", "new-b"), 6))) {
    newdata = data
    newdata$id = ids
    y = fit$simulate(fit, newdata, cp_1 = 6, Intercept_1 = 2, Intercept_1_id_sd = 1,
      sigma_1 = 1, .type = "fitted")
    expect_equal(length(unique(y)), 2L)
    expect_equal(as.numeric(y[1:2]), as.numeric(y[7:8]))
  }
  expect_error(mcp(fit$model, data, par_x = "x", sample = FALSE,
    prior = list(Intercept_2_id_sd = "dnorm(0,1)")), "Set prior on `Intercept_1_id_sd` instead of `Intercept_2_id_sd`")
  bare = mcp(list(y ~ 1 + (1 | id), ~ 0 + (1 | id)), data, par_x = "x", sample = FALSE)
  expect_equal(nrow(mcp_pars(bare, role = "group_sd")), 2L)
})


test_that("whole group slopes reuse local coordinates, endpoints, and explicit sources", {
  data = data.frame(x = 1:12, y = 0, id = rep(c("a", "b"), 6))
  fit = mcp(list(y ~ 1 + (0 + x || id), ~ 0 + same((0 + x || id)), ~ 0),
    data, par_x = "x", sample = FALSE)
  values = fit$simulate(fit, data, cp_1 = 4, cp_2 = 8, Intercept_1 = 1,
    x_1_id_sd = 1, sigma_1 = 1, .type = "fitted")
  b = attr(values, "simulated")$x_1_id
  expect_equal(as.numeric(values), 1 + b * pmin(data$x, 8))
  expect_error(mcp(list(y ~ 1 + (1 + x || id), ~ 0 + same((0 + x || id))),
    data, par_x = "x", sample = FALSE), "complete active group block")
  expect_error(mcp(list(y ~ 1 + (1 | id), ~ 0, ~ same((1 | id))),
    data, par_x = "x", sample = FALSE), "complete active group block")
  expect_no_error(mcp(list(y ~ 1 + (1 | id), ~ 0, ~ same((1 | id), as = 1)),
    data, par_x = "x", sample = FALSE))
  expect_error(mcp(list(y ~ 1 + (1 | id), ~ (0 + x || id) + same((1 | id))),
    data, par_x = "x", sample = FALSE), "Only one predictor group-level term")
})


test_that("posterior prediction and likelihood reuse fitted and new group deviations", {
  data = data.frame(x = 1:18, id = rep(c("a", "b", "c"), 6), y = sin(1:18))
  fit = suppressWarnings(mcp(list(y ~ 1 + (1 | id), ~ 0 + same((1 | id))),
    data, par_x = "x", iter = 40, warmup = 30, chains = 1, diagnostics = FALSE, quiet = TRUE))
  draws = as.matrix(coda::as.mcmc(fit))
  expected = vapply(seq_len(nrow(data)), function(i) stats::dnorm(data$y[i],
    draws[, "Intercept_1"] + draws[, paste0("Intercept_1_id[", data$id[i], "]")],
    draws[, "sigma_1"], log = TRUE), numeric(nrow(draws)))
  expect_equal(unname(log_lik(fit, summary = FALSE)), expected)
  for (id in c("a", "new")) {
    nd = data.frame(x = c(2, 17), id = id)
    pred = fitted(fit, newdata = nd, summary = FALSE)
    expect_equal(pred$.epred[pred$data_row == 1], pred$.epred[pred$data_row == 2])
  }
})
