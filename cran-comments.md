
# mcp 0.2

## Notes for the reviewer
* This is the first submission of the `mcp` package and my first CRAN submission personally. I have done my best to adhere to all standards. Extensive documentation for `mcp` is available at https://lindeloev.github.io/mcp/.
* I include 'patchwork' as a dependency even though it is not currently used. It will be used in the near future and users should not install new packages then. Background: mcp worked with dev version of patchwork, but fails for the version just uploaded to CRAN. It is not self-evident whether the bug resides with `patchwork` or `bayesplot`, but I will (have) raise appropriate issues to track down this issue.
* My email will remain for a long time. I co-own the domain.
* `mcp` uses the GPL-2 license. The only code copied (verbatim) from other packages are R/lme4_utils.R which is GPL (>=2).
* `mcp` does not make any external changes (files, options, communication, etc.)
* `mcp` uses JAGS for sampling and the tests will fail without it. This is true for other CRAN packages, including rjags, R2JAGS, etc. I'm unsure how CRAN submission works when a package relies on an external binary, so I'd appreciate feedback if there's an issue. JAGS is available for all platforms and automatic installation for Mac and OSX is set up in the .travis.yml (https://github.com/lindeloev/mcp/blob/master/.travis.yml). Windows binaries for JAGS are here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/

## Test environments
* local Windows 10, R 3.6.1
* Ubuntu 16.04.6 LTS (on travis-ci): release
* Ubuntu 16.04.6 LTS (on travis-ci): devel
* Mac OS X 10.13.3 (on travis-ci): release
* devtools::check_win_release()
* devtools::check_win_devel()

## R CMD check results
There were no ERRORs or WARNINGs.

## Downstream dependencies
This is the first submission so there are no downstream dependencies.
