# Playbook: reweighting timeseries to protect a behaviour

Use this when a sensitivity or GA rewards runs that lose something the model must keep, or
before any GA when weights have changed. Everything here rescores **existing** runs: run
output does not depend on weights, so no new Ecospace runs are needed.

## 1. State the behaviour and its metric
One behaviour at a time, with a number that can be computed from run output:
- a trend → end/start biomass ratio and log-linear slope (%/yr) of the predicted group;
- recruitment or index variability → SD ratio (pred/obs, log scale) and correlation per
  index;
- a spatial pattern → the relevant map term or a regional ratio.
Compute it for the observed data too, so "kept" has a target.

## 2. Find the runs that break it
From the sensitivity, past GA bests and quick GAs: runs whose metric moves away from the
target (e.g. the stock no longer declines) and still improve LL. Add a **flat run** (base
output with every year held at the first year's value; diagnosing.md) - it is the
cleanest "behaviour removed" case. If no run breaks the behaviour, the weights are not the
problem; stop here.

Check the rank correlation of ΔLL with the metric across all sensitivity runs. If runs that
remove the behaviour tend to score better, the objective is rewarding its removal.

## 3. Ask which series can see it
For each candidate series:
- **Does its prediction move at all?** A series whose prediction is set by a forcing can't
  discriminate. In Ecospace, F (type 104) follows the effort timeseries, so F series barely
  differ between parameter sets; catch differs, because catch = F × B.
- **Absolute or relative?** A relative series (q estimated) sees shape, not level. Check
  `obs.ts$ts.head[, c("Title","Type","Poolcode","Absolute","Weight")]` after setup: type
  104 is absolute only for the lead stanza, and setup options can force types relative.
- **Which way does it pull?** Score base, the breakers and the flat run at **unit
  weights** (`return = "detail"`); per series `cost(breaker) − cost(base)`. Positive = it
  penalises the break. The size of that unit signal tells you what weight it needs to
  matter against the series pulling the other way.
- **Combinations pin what single series can't.** If F is forced and catch is absolute,
  then fitting both pins biomass. Look for pairs like this when no single series sees the
  behaviour.

Also look at the *other side*: which series reward the break, and by how much. If they are
weighted from data CVs, they stay; the protecting series must outweigh them.

## 4. Scan weight sets on the existing runs
1. Matrix U = every run × series cost at unit weights (one parallel scoring pass).
2. For weights w: `ts(w) = U %*% (log(2) / log(1 + 1/w^2))` (0 for w = 0), plus each run's
   spatial term (weight-independent). First check that `ts(current w)` reproduces the real
   LL.
3. For each candidate set report: base LL; ΔLL of each breaker, the flat run and past GA
   bests; how many runs improve by > 2 noise SD; the rank correlation of ΔLL with the
   metric.
4. Choose the **lightest** set at which breakers no longer gain and the correlation turns
   the right way. Heavier weights make the behaviour series dominate everything else.

## 5. Sanity-check the chosen weights
- Translate to CV (`cv = 1/w`): landings ~5–10% → weight 10–20; assessment-derived
  biomass/F ~10–20% → 5–10; survey indices from their own CVs. Flag anything outside that
  and say why it is still chosen.
- Don't touch weights derived from data CVs to make the fit behave; protect the behaviour
  by adding weight where it is visible.
- Expect the base LL and the noise SD to change; remeasure noise and reset `min.dLL` and
  the GA stop tolerance.

## 6. Apply and verify
Rescore the full sensitivity under the new weights (keep the original `LL_<tag>.rds`),
check the ranking (runs that break the behaviour should drop; runs moving toward the
target should rise), update noise files and `run_info.csv`. During and after the GA,
report the metric for every best: reweighting worked only if the GA best keeps the
behaviour.

## What doesn't work
- Raising the weight of a series that can't see the behaviour (forced or relative).
- A trend or penalty term added to the objective, when the data already contain the
  behaviour: fix the weights first.
- Expecting weights to keep variability whose **timing** the model gets wrong: with
  uncorrelated timing the lognormal cost always prefers a flatter prediction
  (likelihood.md). Fix timing through parameters or drivers.
