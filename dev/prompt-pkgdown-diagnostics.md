# Detect errors and warnings in pkgdown and examples
Review all output blocks in the pkgdown site vignettes and function examples. Report back the following:

 - Any stop()-like errors
 - Any warning()-like warnings
 - Any diagnostics firing from {mcp}, e.g., `mcp()` or mcp-summary saying parameters had low ESS or high rhat.
 - Any diagnostics firing from {loo}, e.g., divergent transitions.
 - Any diagnostics firing from other functions/packages.
 - If simulated values are not recovered. This can be OK if it is just one or two out of many, but otherwise flag it.
 - Any messages that does not seem to be purposeful.


List them all. For each, propose a solution. In general, suppressing warnings is the last resort - the site and examples is a good test in their own respect and new warnings or diagnostics is a signal that should be captured. Do not implement anything.

For chunks currently using chunks using warning=FALSE/message=FALSE, re-knit into the scratchpad with warnings and messages forced on, to see what those chunks hide. If some of these flags are unnecessary, mention that too.

Rank by your guessed severity. Make a separate section for findings where the surrounding text (code comments; plain text) reference this diagnostic firing. This probably needs no action from the user, but it is nice to have it included in your "report".