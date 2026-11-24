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

  # Build manual (to test it works)
  devtools::build_manual(path = 'dev')

  # Test leftovers are .gitignore'd but not .Rbuildignore'd, so they end up in the tarball
  leftovers = c(
    list.files("tests", pattern = "\\.Rout\\.fail$", full.names = TRUE),
    list.files("tests/testthat/_snaps", pattern = "\\.new\\.", recursive = TRUE, full.names = TRUE)
  )
  if (length(leftovers) > 0)
    warning("Delete these test leftovers before building for CRAN: ", paste(leftovers, collapse = ", "))
}