# Objective function, weights and noise

## fn.ecospace_objfxn

Per series: `cv = 1/weight`, `sd = sqrt(log(1 + cv^2))`,
`cost = sum(log(pred/obs)^2) / (2 sd^2)` (relative series: `pred` scaled by the MLE q
first). Weight-0 series are recorded but cost nothing. `LL.total` = timeseries cost +
spatial cost (`spatial.weight` × the GFISHER/sdmTMB map term, if `spatial.obs` is set).

- `return = "vector"`: `LL.total`, `LL.<group>`, `LL.spatial.<map>`.
- `return = "detail"`: one row per series (`Title`, `Type`, `Weight`, `loglik`, ...) - use
  this to attribute gains.
- `return = "both"`.
- `spatial.obs = NULL` scores the timeseries only (fast; the spatial term does not depend on
  timeseries weights).
- Supported obs types: 0, 1 (biomass, relative/absolute), 6 (catch, always absolute),
  61 (catch, relative), 12, 13 (discards: **total** discards), 19, 20, 104 (F = catch/biomass).

## Absolute vs relative

Decided in setup, not by the type code alone:
- `settings$relative.types` (e.g. 12/13/19) forced relative;
- type 104 is absolute only for the **lead stanza** of a multi-stanza group; other stanzas'
  F are relative (a weight on them constrains shape, not level);
- type 6 absolute regardless of `fit.abs.catch`; 61 is its relative twin.

Always print `obs.ts$ts.head[, c("Title","Type","Poolcode","Absolute","Weight")]` for the
series you are reasoning about.

## Weights as precisions

`w` | 1 | 2 | 3 | 5 | 7 | 10 | 15 | 20
---|---|---|---|---|---|---|---|---
cost multiplier vs w=1 | 1 | 3.1 | 7.3 | 18 | 34 | 70 | 157 | 278

Multiplier = `log(1+1) / log(1+1/w^2)`. Because cost is linear in `1/sd^2` and q does not
depend on the weight, **any weight set can be evaluated from one unit-weight scoring**:
score every run once with all weights = 1 (`detail`), then
`cost(w) = cost(1) * log(2) / log(1 + 1/w^2)`. Verified exactly against direct
scoring. This makes weight scans over hundreds of runs take seconds.

How to choose weights for a purpose: [reweighting.md](reweighting.md).

## Why lognormal fits prefer flat predictions

For one series, `E[Σ (log p − log o)²] = n [var(log p) + var(log o) − 2 cov + bias²]`.
If the model's year-to-year timing is uncorrelated with the data (cov ≈ 0), any predicted
variance only adds cost, so the fit improves by **damping** the prediction toward a flat
line. Variability is rewarded only to the extent the timing matches (cov > 0). Report the
SD ratio and correlation for indices whose variability matters; a weight change can't fix
this, only better timing can.

## Run-to-run noise

Ecospace IBM movement is stochastic: the same cmd file gives different LLs. Measure it
(≥10 base replicates) per phase - it changes with weights and model.

- `min.dLL` (sensitivity cutoff) ≈ 2 SD.
- GA stop tolerance ≈ 1–2 SD; the GA uses `stop.rtol × |best|`, so check what that is in LL
  units at the expected best LL.
- A GUI-vs-console gap is "within noise" if |z| < 2 against the replicate mean and SD.
- Look at the per-group SD too: the noise usually sits in one or two groups (often a
  juvenile stanza with strong IBM movement).
- Noise changes when weights change: remeasure (or rescore the replicates) after a
  reweight.
