# Executed by pkgdown before running reference examples
Sys.setenv(IN_PKGDOWN = "true")

if (requireNamespace("bayesplot", quietly = TRUE)) {
  bayesplot::bayesplot_theme_set(bayesplot::theme_default(base_size = 11, base_family = "sans"))
}
