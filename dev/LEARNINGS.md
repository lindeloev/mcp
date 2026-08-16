# Learnings
This files contains learnings that took a lot of work to reach and some of which meant ideas were abandoned, so it is not "documented" in code/articles. 

See DECISIONS.md for similar considerations.

## JAGS sensitivity to changepoint-parameter values
Even small rounding-to-7-decimal-points can catastrophically collapse JAGS sampling of cp-parameters.
See commit 9b887baacc018b5dd55016366e31e84ee548ed24 and docs for register_jags_constant() and jagsify_constants().

## Poor sampling unless mean-centering group-level effects
Mean-centering individual realized effects (participants; IDs; etc.) is necessary because otherwise ESS/sec drops ~10-fold. I think this may be a deviation from brms and similar.