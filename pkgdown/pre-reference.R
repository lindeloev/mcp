# Executed by pkgdown before running reference examples
Sys.setenv(IN_PKGDOWN = "true")
options(width = 9999)  # Wide tables with a horizontal scroll bar, like R.options in vignettes

if (requireNamespace("bayesplot", quietly = TRUE)) {
  bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = 11, base_family = "sans"))
}
