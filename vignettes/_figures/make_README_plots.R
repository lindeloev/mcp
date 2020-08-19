`%>%` = magrittr::`%>%`

theme_it = function(x, title) {
  x +
    ggplot2::ggtitle(title) +
    ggplot2::theme_gray(13) +
    ggplot2::theme(legend.position = "none")  # +
  # theme(axis.title = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank())
}

save_it = function(filename) {
  ggplot2::ggsave(paste0("./vignettes/_figures/", filename), width=6, height=3, dpi = 100, type = "cairo")
}



#############
# Example 1 #
#############
library(mcp)
options(mc.cores = 3)  # Run in parallel

ex_demo = mcp_example("demo")
fit_demo = mcp(model_demo$model, data = ex_demo$data, adapt = 3000)  # dataset included in mcp
theme_it(plot(fit_demo), "")
save_it("ex_demo.png")

plot_pars(fit_demo, regex_pars = "cp_")
save_it("ex_demo_combo.png")

# LOO
# Fit the model
model_null = list(
  response ~ 1 + time,
  ~ 1 + time
)
fit_null = mcp(model_null, ex_demo$data, adapt = 3000)

# Compare loos:
fit_demo$loo = loo(fit_demo)
fit_null$loo = loo(fit_null)
loo::loo_compare(fit_demo$loo, fit_null$loo)


################
# Two plateaus #
################
ex_intercepts = mcp_demo("intercepts")
fit_plateaus = mcp(ex_intercepts$model, ex_intercepts$data, par_x = "x", adapt = 3000)
theme_it(plot(fit_plateaus, lines = 25), "Two plateaus")
save_it("ex_plateaus.png")



########################
# VARYING SLOPE CHANGE #
########################
ex_varying = mcp_example("varying")
fit_varying = mcp(ex_varying$model, ex_varying$data, adapt = 3000)
theme_it(plot(fit_varying, facet_by = "id"), "Varying slope change")
save_it("ex_varying.png")


############
# BINOMIAL #
############
ex_binomial = mcp_example("binomial")
fit_binomial = mcp(ex_binomial$model, ex_binomial$data, family = binomial(), adapt = 5000)
theme_it(plot(fit_binomial, q_fit = TRUE), "Binomial")
save_it("ex_binomial.png")


##########################
# FIXED, RELATIVE, PRIOR #
##########################
ex_rel = mcp_example("rel_prior")
prior_rel = list(
  int_1 = 10,  # fixed value
  x_3 = "x_1",  # shared slope in segment 1 and 3
  int_2 = "dnorm(0, 20)",
  cp_1 = "dunif(20, 50)"  # has to occur in this interval
)

fit_rel = mcp(ex_rel$model, ex_rel$data, prior_rel, iter = 10000)
theme_it(plot(fit_rel, cp_dens = FALSE), "rel() and prior")
save_it("ex_fix_rel.png")


#############
# QUADRATIC #
#############
ex_quadratic = mcp_example("quadratic")
fit_quadratic = mcp(model_quadratic, ex_quadratic, adapt = 3000)
theme_it(plot(fit_quadratic), "Quadratic and other exponentiations")
save_it("ex_quadratic.png")



#################
# TRIGONOMETRIC #
#################
ex_trig = mcp_example("trigonometric")
fit_trig = mcp(ex_trig$model, ex_trig$data, adapt = 3000)
theme_it(plot(fit_trig), "Trigonometric for periodic trends")
save_it("ex_trig.png")



############
# VARIANCE #
############
ex_variance = mcp_example("variance")
fit_variance = mcp(ex_variance$model, ex_variance$data, iter = 10000, adapt = 10000)
theme_it(plot(fit_variance, q_predict = TRUE), "Variance and prediction intervals")
save_it("ex_variance.png")



#########
# AR(N) #
#########
ex_ar = mcp_example("ar")
fit_ar = mcp(ex_ar$model, ex_ar$data, adapt = 3000)
theme_it(plot(fit_ar), "Time series with autoregressive residuals")
save_it("ex_ar.png")
