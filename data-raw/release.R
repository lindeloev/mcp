build_release_stuff = function() {
  # Reproducible state
  previous_plan = future::plan(future::sequential)
  on.exit(future::plan(previous_plan), add = TRUE)

  # Latest version of mcp
  devtools::install()
  devtools::load_all()

  # Update documentation
  devtools::document()

  # Build demo fit
  source("data-raw/make_demo_fit.R")
  make_demo_fit()

  # User-facing docs: README and showcase figure
  rmarkdown::render("README.Rmd", output_format = "github_document")
  source("vignettes/_figures/mcp_showcase.R")
}