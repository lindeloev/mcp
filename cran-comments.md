
# mcp 0.2

## Notes for the reviewer
* This is the first submission of the `mcp` package and my first CRAN submission personally. I have done my best to adhere to all standards. Extensive documentation for `mcp` is available at https://lindeloev.github.io/mcp/.
* The package `patchwork` is a dependency which just arrived on CRAN a few days ago. It seems it has not rolled out to all test servers. Travis and devtools::check_rhub() install it just fine.
* `mcp` uses JAGS (an external binary) for sampling through the `rjags` package. rjags will fail to install without JAGS on the system. This happens when I run `devtools::check_win_release()` and `check_win_devel()`. Travis install JAGS prior to installing packages, and all tests pass. They do too on my Windows PC. Automatic installation for Mac and OSX is set up in the .travis.yml (https://github.com/lindeloev/mcp/blob/master/.travis.yml). Windows binaries for JAGS are here: https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/
* My email address is long-term. I have had it for 15 years and I co-own the domain.
* `mcp` uses the GPL-2 license. The only code copied (verbatim) from other packages are R/lme4_utils.R which is GPL (>=2).
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
