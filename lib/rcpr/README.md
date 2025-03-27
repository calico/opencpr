# Calico Pathway Resources (CPR) in R!

An R package which is part of the [Calico Pathway Resources (CPR)](https://github.com/calico/opencpr) project.

## Repo Contents

- Functions for working with bioconductor's organism-specific packages.
- Functions for querying Open Target's graphQL database to summarize targets and their top indications.
- Shiny modules and related functions for exploring networks interactively.

## Setup directions

This package can be installed in several ways:

1.  Install from GitHub with a PAT:
    1.  Add "GITHUB_PAT" as a secret in your .Renviron. This key can be obtained through your user settings in GitHub. (This instruction can be dropped once rcpr is open-sourced.)
    2.  Install `rcpr` in an interactive R session using remotes (or devtools) from GitHub:

```r
remotes::install_github("calicolabs/opencpr",  subdir = "lib/rcpr)
```

2.  Clone and build locally:

```bash
!bash
git clone git@github.com:calico/opencpr.git
R CMD build lib/rcpr
```
