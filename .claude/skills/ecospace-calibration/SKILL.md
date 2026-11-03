---
name: ecospace-calibration
description: Calibrate an Ecopath with Ecosim (EwE) Ecospace model with R4EwE - setting up a calibration phase, checking a GUI base against the console, measuring run-to-run noise, running and reading the parameter sensitivity, choosing and testing timeseries weights, configuring/launching/monitoring the genetic algorithm (fn.GA), and diagnosing what a parameter change or a GA result did to the fits. Use whenever the work involves fn.ecospace_objfxn, fn.GA, a sensitivity runlist, Ecospace console runs, or timeseries weights for an EwE model, for any model.
---

# Ecospace calibration with R4EwE

Calibration is a loop of **phases**. Each phase is one set of inputs (model, cmd base, obs
timeseries with weights, env-response pars, settings), one sensitivity, and one or more GAs.
A phase ends when its fits are reviewed and the next phase's changes are decided.

Two playbooks carry most of the judgement. Use them every phase:

| playbook | use it when |
|---|---|
| [references/reweighting.md](references/reweighting.md) | choosing or changing timeseries weights; a fit improves by losing a behaviour the model must keep |
| [references/diagnosing.md](references/diagnosing.md) | explaining a sensitivity ranking, a GA best, or any LL change; deciding whether a gain is real |

Supporting references - read the one that matches the step, not all up front:

| topic | reference |
|---|---|
| objective function, weights as precisions, noise | [references/likelihood.md](references/likelihood.md) |
| EwE / Ecospace behaviour that surprises people | [references/ewe_mechanics.md](references/ewe_mechanics.md) |
| GA setup, stopping, monitoring, failure modes | [references/ga.md](references/ga.md) |
| both playbooks applied to one real problem | [references/example_wfs_mice_gag.md](references/example_wfs_mice_gag.md) |
| general lessons with their evidence and status | [references/lessons.md](references/lessons.md) |
| ad hoc analyses that should become R4EwE functions | [references/r4ewe_roadmap.md](references/r4ewe_roadmap.md) |

Model-specific facts (group numbers, which series matter, past LLs, chosen weights) belong
in that model's project memory, **not** in this skill.

## The phase loop

1. **Set up the phase.** `ph <- fn.calib_phase(phase, inputs=, refs=, settings=)`; inputs are
   edited in `ga_results/<phase>/inputs/`. Every run gets a label:
   `fn.calib_paths(ph, label)` → `fn.calib_start()` (copies inputs, writes `run_info.csv` with
   md5s) → … → `fn.calib_finish()`. Scripts read `ph$settings$...`, never literals.
2. **Name the behaviours to protect.** With the user, list what the model must reproduce
   whatever the LL says: a documented stock decline, recruitment variability, a spatial
   gradient, a collapse and recovery. Each one gets a metric (e.g. end/start biomass ratio,
   SD and correlation of an index) that is reported for every best run from here on.
3. **Check the GUI base against the console.** Remind the user to **save the model in the
   GUI first**. Score both with `fn.ecospace_objfxn(..., return="vector")` and compare group
   by group against measured noise (step 4), not by eye.
4. **Measure noise.** Run the unchanged base ≥10 times. The SD of `LL.total` sets `min.dLL`
   (≈2 SD) and the GA stop tolerance (≈1–2 SD). Note which groups carry the noise.
5. **Sensitivity.** One-at-a-time runs: each vul pair low/high, dispersal ±50%, fleet
   effective power / effort ±50%, each env response off. Rank by ΔLL vs base, then read it
   with **diagnosing.md**: what the data want, and which top runs break a protected
   behaviour.
6. **Check the weights before any GA** with **reweighting.md**. Run outputs don't depend
   on weights, so this is a rescore of existing runs, not a rerun. If the weights change,
   rescore the sensitivity, swap in the rescored `LL.rds` (keep the original as
   `LL_<tag>.rds`), and update the sens `inputs/` copy, noise files and `run_info.csv`.
7. **GA.** Quick runs first (`gaq01`, … small popSize, ~6–10 gens) to test bounds; then the
   full run (`ga01`). The user launches multi-hour jobs from their own R session: prepare
   the script, give the `source()` line. See ga.md.
8. **Monitor.** Every generation: table, interpretation, verdict, push notification, and
   the protected-behaviour metrics of the current best. Keep every report complete.
9. **Review the result** with **diagnosing.md**: attribute, find the mechanism, check for
   escapes. Plot fits (`fn.ecospace_plot_fits`, `fn.ecospace_plot_series_fits`). Then
   archive (`archive_ga`), slim-copy the final gen, move scratch to the archive drive.
10. **Close the phase** - see below.

## Decision rules that hold across models

- **The LL is not the only judge.** A better LL that comes from removing a protected
  behaviour is a weighting or structure problem, not a better model. Say so.
- **Decompose before explaining.** Rescore per series and attribute the gain before
  proposing a cause. Guesses about "what drove the gain" are often wrong.
- **Ask which series can see the behaviour** - and whether they are absolute or relative,
  and whether a forcing fixes their prediction. Weight on a blind series does nothing.
- **Weights are precisions.** weight = 1/CV; cost scales as `1/log(1+1/w²)`. Translate a
  proposed weight to a CV and ask whether it is defensible for that data source. Weights
  derived from data CVs are left alone unless the user says otherwise.
- **Test before launching.** Any weight change is checked against the runs that break a
  protected behaviour, by rescoring, before a GA spends days on it.
- **When one escape is closed, look for the next.** After every GA, check parameters at
  bounds and far from base for a new route around the protected behaviour.
- **Prefer simple, direct script code**; shared logic goes into R4EwE as an exported
  function (lift and export, never `:::`), tested with the package *installed*.
- **Don't commit or push unless asked.**

## Closing a phase

1. Ask: **would this help on a different model?** Only those items go into the skill:
   - a better way to reason (a step, a check, a symptom → cause row) → the playbook it
     belongs to;
   - a fact about EwE, R4EwE or the likelihood → `lessons.md` (with one line of evidence),
     and into the mechanics/likelihood reference once confirmed;
   - a recipe done by hand twice → `r4ewe_roadmap.md`.
   Everything about this model (its groups, numbers, weights, decisions) → project memory.
2. Draft the edits and show them to the user; apply only what they approve.
3. A lesson seen in one model is `tentative`; `confirmed` when it holds in a second model
   or is proven from EwE source, the maths, or a dedicated experiment. Contradicted lessons
   are marked `retracted` with the reason.
4. Every few models, prune: merge duplicates, drop what R4EwE now does, keep this file short.
   Changes are reviewed as git commits in the R4EwE repo.
