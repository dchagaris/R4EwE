---
title: "Building environmental response shapes from AquaMaps"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Building environmental response shapes from AquaMaps}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  eval = FALSE          # AquaMaps DB (~2 GB) is not available on the build machine
)
```

## Overview

`fn.make_aquamaps_shapes()` turns a list of functional groups and their
representative species into EwE/Ecospace **environmental response shapes**. For
every species it can match in the [AquaMaps](https://www.aquamaps.org) database,
it reads the species' environmental envelope (HSPEN) and emits a **trapezoid
(EwE function type 9)** shape for each of five drivers:

* depth
* temperature
* salinity
* distance to shore
* dissolved oxygen (DO)

The four trapezoid parameters for each shape are taken from the species'
`Min`, `PrefMin`, `PrefMax`, and `Max` envelope values for that driver.

The result is a `shapes.out` data frame ready to be combined with shapes from
other sources (e.g. GAM-derived shapes) and written to an EwE response-shape
parameter file.

## Prerequisites

The AquaMaps workflow depends on two packages that are **not** hard
dependencies of R4EwE:

* [`aquamapsdata`](https://github.com/raquamaps/aquamapsdata) (GitHub-only)
* `dplyr` (CRAN)

The first time you call `fn.make_aquamaps_shapes()`, it checks for these and, if
they are missing, **prompts you to install them**. Because `aquamapsdata` is
built from source, it also verifies that a compiler toolchain is present
(**Rtools** on Windows) before installing, and points you to the correct
installer if it is not.

To install everything ahead of time:

```{r}
# Windows only: first install Rtools matching your R version from
#   https://cran.r-project.org/bin/windows/Rtools/
install.packages("remotes")
remotes::install_github("raquamaps/aquamapsdata", dependencies = TRUE)
install.packages("dplyr")
```

## The AquaMaps database

`aquamapsdata` ships **no data** — the species database (`am.db`, ~2 GB) is
downloaded separately and, by design, stored at a fixed location managed by the
package. Download it **once**:

```{r}
library(aquamapsdata)
download_db()        # ~2 GB, one-time download to the package's default location
```

### Telling R4EwE where `am.db` lives

`fn.make_aquamaps_shapes()` (and the helper `fn.connect_aquamaps_db()`) accept a
`db_path` argument:

* **`db_path = NULL` (default)** — use the package's standard location, the same
  place `download_db()` writes to. If no database is found there, you get an
  error pointing you to `download_db()`.
* **`db_path = "…/am.db"`** — use a copy stored anywhere you like (for example a
  shared/network drive), so you don't have to re-download the 2 GB file on every
  machine. The file is opened read-only.

## Input file

The species list is a CSV with (at minimum) these columns:

| number | name              | species                |
|-------:|-------------------|------------------------|
| 1      | offshore dolphins | Stenella longirostris  |
| 2      | coastal dolphins  | Tursiops truncatus     |
| 3      | seabirds          |                        |
| 4      | pelagic sharks    | Carcharhinus plumbeus  |

* `number` — functional group number
* `name` — functional group name
* `species` — scientific name used to query AquaMaps (blank rows are skipped)

Any additional columns are ignored.

## Usage

```{r}
library(R4EwE)

# default database location
shapes.out <- fn.make_aquamaps_shapes("species list for aquamaps query.csv")

# or a database kept on a shared drive
shapes.out <- fn.make_aquamaps_shapes(
  "species list for aquamaps query.csv",
  db_path = "D:/data/aquamaps/am.db"
)

head(shapes.out)
```

Species in the list that AquaMaps does not recognize are reported via a message
and omitted from the output.

## Output

`shapes.out` has 11 columns:

| column          | meaning                                            |
|-----------------|----------------------------------------------------|
| `Function name` | `<fg num>_<fg name>_<var>_aqm`                      |
| `Function type` | `9` (trapezoid)                                     |
| `Param 1`–`4`   | trapezoid points: Min, PrefMin, PrefMax, Max       |
| `Param 5`–`6`   | unused (`NA`) for trapezoids                        |
| `var`           | driver: depth, temp, salinity, dist2shore, DO      |
| `fg num`        | functional group number                            |
| `source`        | `"aqm"`                                             |

The parameter columns are numeric and rounded to 3 decimal places, so the data
frame can be written straight to CSV:

```{r}
write.csv(shapes.out, "response shape pars - aquamaps.csv", row.names = FALSE)
```
