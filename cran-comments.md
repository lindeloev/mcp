# mcp 0.4.0

## Notes for the reviewer
* This is a major update to `mcp`, which has been on CRAN since 2020.
* This release introduces API changes to align `mcp` with R-native conventions (generics, `{posterior}` and `{tidybayes}` and `{brms}` standards). Parameter estimates and core modeling remain consistent with v0.3.4. Breaking changes include informative deprecation messages and backwards-compatible soft fallbacks where appropriate.
* The license changes from GPL-2 to GPL (>= 2). I am the sole author and copyright holder.

## Test environments
* Local: Ubuntu Linux, R 4.6.1
* win-builder: Windows Server (release, devel)
* GitHub Actions: Ubuntu (oldrel, release, devel), macOS (release), Windows (release, devel)

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.

## Downstream dependencies
I checked all 3 reverse dependencies with `revdepcheck`. None had new problems. I notified all maintainers on 2026-08-09 and sent fixes where relevant:
* `glossa` (Imports): calls `mcp::phi()`, which is now deprecated with a warning but still works. GitHub pull request sent.
* `smoothbp` (Suggests): used only in a vignette. The maintainer has replied and will update it when 0.4.0 is released.
* `ggchangepoint` (Suggests; one of ~60 suggested packages): GitHub pull request sent.