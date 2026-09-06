# Ask reminder questions for CRAN export
release_questions = function() {
  c(
    # Do before merging into dev-branch:
    "TEST: Have you run the release-level fit recovery tests? Sys.setenv(MCP_TEST_LEVEL = 'release', TESTTHAT_CPUS = 11); devtools::test()",
    "TEST: Have you run LLM checks using all dev/prompt-*.R and thought about whether the findings should be addressed or added to dev/DECISIONS.md?",
    "TEST: Have you manually reviewed all mcp_example() plots?",

    "DOC: Have you re-built the site using pkgdown::build_site(lazy=FALSE) and checked it locally?",
    

    # Do before merging into main-branch:
    "TEST: Have you run `revdepcheck::revdep_check()` and notified?",
    "TEST: have you run devtools::check_win_devel(); devtools::check_win_release()?",
    "TEST: have you run urlchecker::url_check()?",
    "TEST: have you run devtools::build_manual(path = tempdir())? Requires tinytex",

    "BUILD: Have you run data-raw/release.R?",

    "DOC: Have you pushed to dev-branch and checked all pages in https://lindeloev.github.io/mcp/dev/ render correctly?"
  )
}