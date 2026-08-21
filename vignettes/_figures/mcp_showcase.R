# ABOUT 
# This plots various mcp models along with their model formulas.
# It is a condensed overview of models supported by mcp. It is used in README
# and other communication.

library(ggplot2)
library(patchwork)
devtools::load_all()

# Insert title and model in faceted plot
style_gg = function(gg, model, title, left, bottom, right, top) {
  # Create model box
  gg_label = ggplot() + 
    geom_label(aes(label = paste0(model, collapse = "\n   "), x = 0, y = 1), hjust = 0, vjust = 1, size = 4, fill = "grey") +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE, clip = "off") +
    theme_void()

  gg + 
    labs(title = title) +
    ggplot2::theme(
      strip.text = ggplot2::element_blank(),
      strip.background = ggplot2::element_blank()
    ) + 

    # Insert model box
    patchwork::inset_element(
      gg_label,
      left = left, bottom = bottom,
      right = right, top = top,
      align_to = "full",
      clip = FALSE
    )
}


# Make the base plots
future::plan(future::multisession, workers = 6)
fit_intercepts = mcp_example("intercepts")
gg_intercepts_org = force(last_plot())  # force to avoid lazy evaluation issues with patchwork
intercepts_model = list(y ~ 1, ~ 1)

fit_binomial = mcp_example("binomial")
gg_binomial_org = force(last_plot())
binomial_model = list(y | trials(N) ~ 1, ~ 0 + x, ~ 1 + x)

fit_sigma = mcp_example("sigma")
gg_sigma_org = force(last_plot()[[1]])
sigma_model = list(y ~ 1, ~ 0 + sigma(1 + x), ~ 0 + x)

fit_ar = mcp_example("ar")
gg_ar_org = force(last_plot()[[1]])
ar_model = list(price ~ 1 + ar(2), ~ 0 + time + ar(1))

fit_group_cp = mcp_example("group_cp")
gg_group_cp_org = force(last_plot()) + scale_x_continuous(breaks = seq(0, 200, by = 25))# + facet_wrap(~id, nrow = 2) + scle_x_continuous()
group_cp_model = list(y ~ 1 + x, ~ 1 + (1|id) ~ 0 + x)

fit_group_mu = mcp_example("group_mu")
gg_group_mu_org = force(last_plot()) + scale_x_continuous(breaks = seq(0, 200, by = 25))# + facet_wrap(~id, nrow = 2)
group_mu_model = list(y ~ 1 + condition + (condition || id), ~ 1 + condition)

fit_multiple = mcp_example("multiple")
gg_multiple_org = force(last_plot())
multiple_model = list(y ~ 1 + x:group + z, ~ 1 + x + group, ~ 0 + I(x^2))


# Style them
gg_intercepts = style_gg(gg_intercepts_org, intercepts_model, "Intercept change", 0.2, 0.75, 0.3, 0.85)
gg_binomial = style_gg(gg_binomial_org, binomial_model, "Binomial changes and other GLM", 0.4, 0.75, 0.5, 0.85)
gg_sigma = style_gg(gg_sigma_org, sigma_model, "Sigma changes w/prediction", 0.5, 0.75, 0.6, 0.85)
gg_ar = style_gg(gg_ar_org, ar_model, "AR(N) change", 0.2, 0.75, 0.3, 0.85)
gg_group_cp = style_gg(gg_group_cp_org, group_cp_model, "Group-level change point", 0.5, 0.75, 0.6, 0.85)
gg_group_mu = style_gg(gg_group_mu_org, group_mu_model, "Group-level mean", 0.3, 0.75, 0.4, 0.85)
gg_multiple = style_gg(gg_multiple_org, multiple_model, "Changes w/multiple regression with interactions", 0.1, 0.75, 0.2, 0.85)


# All together via patchwork
layout = "
AAABBB
CCCDDD
EEEFFF
GGGGGG
"

mcp_showcase = gg_intercepts + gg_binomial +
  gg_sigma + gg_ar +
  patchwork::wrap_elements(full = gg_group_cp, clip = FALSE) +
  patchwork::wrap_elements(full = gg_group_mu, clip = FALSE) +
  gg_multiple + 
  plot_layout(design = layout, heights = c(0.9, 1, 1.7, 1.4)) + 
  plot_annotation(
    title = "Regression with change points using {mcp}",
    theme = theme(plot.title = element_text(size = 19))
  )

# Save to .gitignored location; will be published 
dir.create("pkgdown/assets", recursive = TRUE, showWarnings = FALSE)  # Ensure folder exists
ggsave("pkgdown/assets/mcp_showcase.png", plot = mcp_showcase, scale = 1.25, width = 7, height = 10)
ggsave("pkgdown/assets/mcp_showcase_small.png", plot = mcp_showcase, scale = 1.25, width = 7, height = 10, dpi = 85)
