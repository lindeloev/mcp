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
  ggplot2::ggsave(paste0("./man/figures/", filename), width=6, height=3, dpi = 100, type = "cairo")
}



#############
# Example 1 #
#############
library(mcp)
options(mc.cores = 3)  # Run in parallel

model_demo = list(
  response ~ 1,  # plateau (int_1)
  ~ 0 + time,  # joined slope (time_2) at cp_1
  ~ 1 + time  # disjoined slope (int_1, time_2) at cp_2
)
fit_demo = mcp(model_demo, data = ex_demo)  # dataset included in mcp
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
fit_null = mcp(model_null, ex_demo)

# Compare loos:
fit_demo$loo = loo(fit_demo)
fit_null$loo = loo(fit_null)
loo::loo_compare(fit_demo$loo, fit_null$loo)


################
# Two plateaus #
################
model_plateaus = list(
  y ~ 1,  # plateau (int_1)
  ~ 1  # plateau (int_2)
)
fit_plateaus = mcp(model_plateaus, ex_plateaus, par_x = "x")
theme_it(plot(fit_plateaus, lines = 25), "Two plateaus")
save_it("ex_plateaus.png")



########################
# VARYING SLOPE CHANGE #
########################

model_varying = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x + sigma(1)  # joined slope, varying by id
)
fit_varying = mcp(model_varying, ex_varying)
theme_it(plot(fit_varying, facet_by = "id"), "Varying slope change")
save_it("ex_varying.png")


############
# BINOMIAL #
############
model_binomial = list(
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,  # joined changing rate
  ~ 1 + x  # disjoined changing rate
)
fit_binomial = mcp(model_binomial, ex_binomial, family = binomial(), adapt = 5000)
theme_it(plot(fit_binomial, q_fit = TRUE), "Binomial")
save_it("ex_binomial.png")


##########################
# FIXED, RELATIVE, PRIOR #
##########################
model_rel = list(
  y ~ 1 + x,
  ~ rel(1) + rel(x),
  rel(1) ~ 0 + x
)
prior_rel = list(
  int_1 = 10,  # fixed value
  x_3 = "x_1",  # shared slope in segment 1 and 3
  int_2 = "dnorm(0, 20)",
  cp_1 = "dunif(20, 50)"  # has to occur in this interval
)

fit_rel = mcp(model_rel, ex_rel_prior, prior_rel, iter = 10000)
theme_it(plot(fit_rel, cp_dens = FALSE), "rel() and prior")
save_it("ex_fix_rel.png")


#############
# QUADRATIC #
#############
model_quadratic = list(
  y ~ 1,
  1 ~ 0 + x + I(x^2)
)
fit_quadratic = mcp(model_quadratic, ex_quadratic)
theme_it(plot(fit_quadratic), "Quadratic and other exponentiations")
save_it("ex_quadratic.png")



#################
# TRIGONOMETRIC #
#################
# Define the model
model_trig = list(
  y ~ 1 + sin(x),
  ~ 0 + cos(x) + x
)
fit_trig = mcp(model_trig, ex_trig)
theme_it(plot(fit_trig), "Trigonometric for periodic trends")
save_it("ex_trig.png")



############
# VARIANCE #
############
model_variance = list(
  y ~ 1,
  ~ 0 + sigma(1 + x),
  ~ 0 + x
)

fit_variance = mcp(model_variance, ex_variance, iter = 10000, adapt = 10000)
theme_it(plot(fit_variance, q_predict = TRUE), "Variance and prediction intervals")
save_it("ex_variance.png")



#########
# AR(N) #
#########
model_ar = list(
  price ~ 1 + ar(2),
  ~ 0 + time + ar(1)
)

fit_ar = mcp(model_ar, ex_ar)
plot_ar = plot(fit_ar)
theme_it(plot_ar, "Time series with autoregressive residuals")
save_it("ex_ar.png")
