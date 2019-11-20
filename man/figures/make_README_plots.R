`%>%` = magrittr::`%>%`

theme_it = function(x, title) {
  x +
    ggplot2::ggtitle(title) +
    ggplot2::theme_gray(13) # +
  # theme(axis.title = element_blank(),
  #       axis.text = element_blank(),
  #       axis.ticks = element_blank())
}

save_it = function(filename) {
  ggplot2::ggsave(paste0("./man/figures/", filename), width=6, height=3, dpi = 100)
}



#############
# Example 1 #
#############
library(mcp)
segments_demo = list(
  response ~ 1,  # plateau (int_1)
  ~ 0 + time,  # joined slope (time_2) at cp_1
  ~ 1 + time  # disjoined slope (int_1, time_2) at cp_2
)
fit_demo = mcp(segments_demo, data = ex_demo)  # dataset included in mcp
theme_it(plot(fit_demo), "")
save_it("ex_demo.png")

plot(fit_demo, "combo", regex_pars = "cp_")
save_it("ex_demo_combo.png")

# LOO
# Fit the model
segments_null = list(
  response ~ 1 + time,
  ~ 1 + time
)
fit_null = mcp(segments_null, ex_demo)

# Compare loos:
fit$loo = loo(fit)
fit_null$loo = loo(fit_null)
loo::loo_compare(fit$loo, fit_null$loo)


################
# Two plateaus #
################
segments_plateaus = list(
  y ~ 1,  # plateau (int_1)
  ~ 1  # plateau (int_2)
)
fit_plateaus = mcp::mcp(segments_plateaus, ex_plateaus, par_x = "x")
theme_it(plot(fit_plateaus, lines = 25), "Two plateaus")
save_it("ex_plateaus.png")



########################
# VARYING SLOPE CHANGE #
########################

segments_varying = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope, varying by id
)
fit_varying = mcp::mcp(segments_varying, ex_varying)
theme_it(plot(fit_varying, facet_by = "id"), "Varying slope change")
save_it("ex_varying.png")


############
# BINOMIAL #
############
segments_binomial = list(
  y | trials(N) ~ 1,  # constant rate
  ~ 0 + x,  # joined changing rate
  ~ 1 + x  # disjoined changing rate
)
fit_binomial = mcp::mcp(segments_binomial, ex_binomial, family = binomial())
theme_it(plot(fit_binomial), "Binomial")
save_it("ex_binomial.png")


##########################
# FIXED, RELATIVE, PRIOR #
##########################
segments_rel = list(
  y ~ 1 + x,
  ~ rel(1) + rel(x),
  rel(1) ~ 0
)
prior_rel = list(
  int_1 = 10,  # fixed value
  x_3 = "x_1",  # shared slope in segment 1 and 3
  int_2 = "dnorm(0, 20)",
  cp_1 = "dunif(20, 50)"  # has to occur in this interval
)

fit_rel = mcp::mcp(segments_rel, ex_rel_prior, prior_rel)
theme_it(plot(fit_rel), "rel() and prior")
save_it("ex_fix_rel.png")


#############
# QUADRATIC #
#############
segments_quadratic = list(
  y ~ 1,
  1 ~ 0 + x + I(x^2)
)
fit_quadratic = mcp(segments_quadratic, ex_quadratic)
theme_it(plot(fit_quadratic), "Quadratic and other exponentiations")
save_it("ex_quadratic.png")



#################
# TRIGONOMETRIC #
#################
# Define the model
segments_trig = list(
  y ~ 1 + sin(x),
  ~ 0 + cos(x) + x
)
fit_trig = mcp(segments_trig, ex_trig)
theme_it(plot(fit_trig), "Trigonometric for periodic trends")
save_it("ex_trig.png")






#############
# FOR TWEET #
#############
library(mcp)
library(ggplot2)
segments = list(
  y ~ 1 + x,
  ~ 0 + x
)
empty = mcp(segments, sample = FALSE)
ex_tweet = tibble::tibble(
  x = 1:100,
  y = empty$func_y(x, int_1 = 10, x_1 = 1, x_2 = -0.5, cp_1 = 30, sigma = 5)
)
fit = mcp(segments, ex_tweet)
plot(fit)

# Binomial
segments_bin = list(
  score | trials(N) ~ 1,  # plateau
  1 + (1 | id) ~ 0 + difficulty  # joined slope
)
empty_bin = mcp(segments_bin, family = binomial(), sample = FALSE)
ex_tweet_bin = tibble::tibble(id = 1:5) %>%
  tidyr::expand_grid(difficulty = rep(1:10, each = 3)) %>%
  dplyr::mutate(
    N = 10,
    score = empty_bin$func_y(difficulty, N, int_1 = 2, difficulty_2 = -0.8, cp_1 = 5, cp_1_id = 1 * (id - mean(id)))
  )
fit_bin = mcp(segments_bin, ex_tweet_bin, family = binomial())
#plot(fit_bin, facet_by="id")
hypothesis(fit_bin, "`cp_1_id[1]` < `cp_1_id[2]`")
