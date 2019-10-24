library(mcp)
library(tidyverse)

theme_it = function(x, title) {
  x +
    ggtitle(title) +
    theme_gray(15) # +
    # theme(axis.title = element_blank(),
    #       axis.text = element_blank(),
    #       axis.ticks = element_blank())
}

save_it = function(filename) {
  ggsave(paste0("doc/", filename), width=7, height=4, dpi = 100)
}

################
# Two plateaus #
################

segments = list(
  y ~ 1,
  1 ~ 1
)
empty = mcp(segments, sample=FALSE, par_x = "x")
data = tibble(x = 1:100,
              y = empty$func_y(x, 20, 50, 20, 40))
fit = mcp(segments, data, par_x = "x")
theme_it(plot(fit), "Two plateaus")
save_it("ex_plateaus.png")




################
# SLOPE CHANGE #
################

segments = list(
  y ~ 1 + x,  # intercept + slope
  1 ~ 0 + x  # joined slope
)
empty = mcp(segments, sample=FALSE)
data = tibble(x = 1:100,
              y = empty$func_y(x, 15, 20, 3, -2, 30))
fit = mcp(segments, data)
theme_it(plot(fit), "Slope change")
save_it("ex_slopes.png")




########################
# VARYING SLOPE CHANGE #
########################

segments = list(
  y ~ 1 + x,  # intercept + slope
  1 + (1|id) ~ 0 + x  # joined slope
)
empty = mcp(segments, sample=FALSE)
data = tibble(id = c("John", "Benny", "Rose", "Cath", "Bill", "Erin")) %>%
  expand_grid(x = seq(1, 100, by=4)) %>%
  mutate(
    id2 = as.numeric(as.factor(id)),
    y = empty$func_y(x, 25, 15, 3, -2, 40, 7*(id2 - mean(id2)))
  )
fit = mcp(segments, data)
theme_it(plot(fit, facet_by = "id"), "(1|group): Varying slope change")
save_it("ex_slopes_varying.png")




##########################
# FIXED, RELATIVE, PRIOR #
##########################
segments = list(
  y ~ 1 + x,
  1 ~ rel(1) + rel(x),
  rel(1) ~ 0
)
prior = list(
  int_1 = "dnorm(0, 20)",
  x_1 = -3,
  cp_1 = "dunif(20, 50)"
)
empty = mcp(segments, sample=FALSE)
data = tibble(x = 1:100,
              y = empty$func_y(x, 30, 20, -3, 20, 5, 30, 45))
fit = mcp(segments, data, prior, n.adapt = 2000, n.update = 2000)
theme_it(plot(fit), "rel() and prior")
save_it("ex_fix_rel.png")
