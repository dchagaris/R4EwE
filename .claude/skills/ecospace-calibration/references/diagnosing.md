# Playbook: diagnosing what a parameter change did

Use this for any LL change worth explaining: a top sensitivity run, a GA best, a weight
change, a GUI edit. The order is fixed: **attribute → mechanism → escapes → verdict.**
Don't explain a gain before step 1.

All of this reads existing run output. Score in parallel (`parallel::makeCluster`,
`library(R4EwE)` on workers, export `obs.ts`, `group.names`, `fleet.names`, `ph`,
`spatial.obs`). Run directories may move from scratch to the archive drive mid-job, so
look in both per run.

## 1. Attribute the gain
Score the reference (base) and the run with `return = "detail"`, `spatial.obs = NULL`
(score the spatial term separately); difference per series; sort. Sum by group and by obs
type. Report:
- the top contributors and the share of the total they carry;
- timeseries vs spatial;
- the series that got **worse** (often the clue to the mechanism).
A gain concentrated in a few series of one group is common. It is not wrong in itself, but
it means the result rests on those series' weights and on what they can see.

## 2. Find the mechanism
Ask what the changed parameters do to the predictions behind the top series:
- **Which parameters moved?** Compare with base (ratio, or in bound widths). For a GA
  best, list every parameter far from base and every one at a bound.
- **What did the predictions do?** For the key groups: biomass trend, catch and F level
  and trend, index SD and correlation. Use the residual check below - the *sign* of each
  residual says which way the run moved.
- **Link the two** with the symptom table. Confirm a suspected mechanism with the one-at-
  a-time sensitivity run of that parameter where one exists.

### Symptom → likely cause

| what you see | likely cause | check |
|---|---|---|
| Gain concentrated in one group of indices | fit rewards smoothness of noisy indices | flat-run test; SD ratio and correlation |
| A flat run scores close to (or better than) the fit | those series reward removing dynamics | which series give the flat run its score |
| Predicted SD ≪ observed, correlation ≈ 0 | damping, not fitting | correlation by lag (−2…+2 yr); timing drivers |
| Catch too high, F matched | biomass held up (productivity, mortality or vul changes) | biomass trend vs observed |
| Catch matched, F too low | fleet / effort parameters lowering F | fleet parameters vs base |
| Ranking barely changes after a reweight | reweighted series have no signal | unit-weight signal; absolute vs relative; forced prediction |
| Big gain from switching a driver off | driver timing or shape wrong, or driver unrealistic | the driver's own data; partial (weakened) runs |
| Parameter at a bound | an escape, or a bound that is too tight | what it does to the protected behaviour; its sensitivity run |
| Many parameters at the same bound | structural mismatch the GA is pushing against | which data series pull them there |
| Spatial term flat while timeseries gains a lot | gain is temporal only | expected if maps are weak; check map residuals if not |
| Few groups off from t = 0 vs GUI | unsaved GUI edits | save and rerun |
| Large sudden LL drop in a GA | truncated output scored low | region file count in the run dir |

## 3. Check for escapes
An escape is a route to a better LL around a protected behaviour. After every sensitivity
and GA:
- report each protected-behaviour metric for the best run vs base and observed;
- for each parameter at a bound or far from base, ask: is the value plausible on its own
  evidence, and does it work by weakening something the data say happens (a mortality
  event, a fishery, a recruitment pulse)?
- when the previous escape was closed (by reweighting), expect the GA to find the
  next-cheapest one; look at what *new* parameters moved compared with the last phase.

## 4. Verdict
One paragraph: what improved, by how much, through which mechanism, whether any protected
behaviour was lost, and which parameters need the user's judgement (bounds to widen or
hold, drivers to revisit). Say plainly when a gain is an escape.

## Recipes

**Flat run.** Copy a base run's `Region_*_(Biomass|Catch|Landings).csv` with every year set
to the first year's values; score it. It shows what "no dynamics" is worth under the
current weights.

**Trend check.** Per run: end/start ratio and log-linear slope of predicted biomass, mean
predicted F (catch / biomass from Region_0 files), the same for observed. Spearman
correlation of ΔLL with the ratio across sensitivity runs: negative = the objective
rewards removing the trend.

**Residual check.** For catch and F (or any pair) of the key group: mean log(pred/obs),
RMS and trend, for base and the runs being explained. Shows *how* a run escapes.

**Index variability.** Per index: SD of log pred / SD of log obs (on overlapping years) and
correlation, at lags −2…+2 years. Low SD ratio with low correlation = damping; good
correlation with low SD ratio = right timing, too little amplitude.

**Parameters vs bounds.** For a GA best: each parameter's value, base value, bounds, and
position within the bounds (0–1). List those within 1% of a bound and those beyond ×3 / ÷3
of base.

**GUI vs console.** Score both under the same obs; z = (GUI − replicate mean) / SD, total
and per group.
