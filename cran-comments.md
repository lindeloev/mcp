# mcp 0.3

## Notes for the reviewer
 * This release adds support for `dplyr` 1.0+ and other packages which caused the prior `mcp` to be taken down from CRAN. Sorry it took so long.
 
 * Please see "Expected NOTEs and ERRORs" in the bottom of this document.

## Test environments
* local Windows 10, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci): release
* Ubuntu 16.04.6 LTS (on travis-ci): devel
* Mac OS X 10.13.3 (on travis-ci): release

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
`mcp` has no downstream dependencies.


# Resubmission 3

 * Deleted call to `options(mc.cores = 3)`.
 
 * See the section "Expected NOTEs and ERRORs" below for anticipated ERRORs and NOTEs.
 


# Resubmission 2

 * Fixed grammatical error in DESCRIPTION.
 
 * mcp now spawns at most 2 cores on CRAN.



# Resubmission
This is a resubmission. I believe I have solved all the points raised in the initial review. All tests pass. In this version I have:

* Added single quotes around 'mcp' in DESCRIPTION.

* Added literature to DESCRIPTION with the theoretical foundation for the computations done in mcp.

* mcp no longer copies code from other packages so no attribution/ctb is required.

* `print()` and `cat()` now only reside within `print()` and `summary()` functions.

* All functions have a \value specified now. This has led me to do many other improvements in the documentation too.

* All examples run now. Some have been enclosed in \donttest() to reduce runtimes.

* I have taken the liberty to add a few API-breaking updates to `mcp` in this resubmission, so that the API is as stable as possible from the initial CRAN release. These are: (1) changed plotting of time-series, (2) the function name to simulate data, and (3) changed the `summary()` output for simulated data.



# mcp 0.2

## Notes for the reviewer
* This is the first submission of the `mcp` package and my first CRAN submission personally. I have done my best to adhere to all standards. Extensive documentation for `mcp` is available at https://lindeloev.github.io/mcp/.
* The package `patchwork` is a dependency which just arrived on CRAN a few days ago. It seems it has not rolled out to all test servers. Travis and devtools::check_rhub() install it just fine.
* My email address is long-term. I have had it for 15 years and I co-own the domain.
* `mcp` uses the GPL-2 license. The only code copied (verbatim) from other packages is in R/lme4_utils.R, is GPL (>=2), and has been given proper attribution via the @authors Roxygen tag.
* `mcp` does not make any external changes (files, options, communication, etc.)

## Test environments
* local Windows 10, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci): release
* Ubuntu 16.04.6 LTS (on travis-ci): devel
* Mac OS X 10.13.3 (on travis-ci): release

## R CMD check results
There were no ERRORs or WARNINGs.

## Downstream dependencies
This is the first submission so there are no downstream dependencies.



# Expected NOTEs and ERRORs

* INSTALL ERROR: `mcp` uses JAGS (an external binary) for sampling through the `rjags` package. rjags will fail to install without JAGS on the system. This happens when I run `devtools::check_win_release()` and `check_win_devel()`. Travis install JAGS prior to installing packages, and all tests pass. They do too on my Windows PC. Automatic installation for Mac and OSX is set up in the .travis.yml (https://github.com/lindeloev/mcp/blob/master/.travis.yml). Windows binaries for JAGS are here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/

* TEST ERROR: One or two of the 1500+ unit tests fail around one in ten times with an error like "loo... tails... missing value where TRUE/FALSE needed". This is an artefact of quick unit testing (few MCMC iterations) and the end user will not experience this.

* DESCRIPTION NOTE: rhub says that the DESCRIPTION DOIs return a HTTP 403 error (forbidden). But the DOI work just fine, e.g., http://doi.org/10.2307/2986119.

* The examples and vignettes call `options(mc.cores = 3)` since this is highly recommended in practice, but the number of cores will be capped at two cores on CRAN (if `Sys.getenv("_R_CHECK_LIMIT_CORES_") == "TRUE"`).
