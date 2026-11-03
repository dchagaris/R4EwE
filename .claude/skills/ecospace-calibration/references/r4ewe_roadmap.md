# Analyses to lift into R4EwE

Recipes done ad hoc in scratch scripts more than once. Each should become an exported,
tested function; then the playbooks call the function instead of describing the recipe.
Remove an entry once it ships (note the commit).

| proposed function | does | playbook step |
|---|---|---|
| `fn.score_runs(dirs, obs.ts, weights = 1, ncores)` | parallel `detail` scoring of many run dirs (scratch → archive fallback); runs × series matrix | diagnosing 1, reweighting 3–4 |
| `fn.attribute_gain(ref, run, obs.ts, groups)` | per-series and per-group ΔLL, sorted, with shares | diagnosing 1 |
| `fn.reweight_LL(U, weights, spatial)` | apply a weight vector to a unit-weight matrix; LL.total per run | reweighting 4 |
| `fn.weight_scan(U, grid, test.runs, metric)` | the scan table | reweighting 4 |
| `fn.rescore_sens(label, obs.ts)` | rescore a sensitivity, write `LL.rds` (keep `LL_<tag>.rds`), update noise files and run_info | reweighting 6 |
| `fn.flat_run(dir.base, dir.out)` | build the flat-first-year pseudo-run | reweighting 2, diagnosing |
| `fn.trend_check(dirs, groups, obs.ts)` | end/start ratio, %/yr slope, mean F, vs observed | reweighting 1, diagnosing 2 |
| `fn.index_variability(dirs, obs.ts, series, lags)` | SD ratio and correlation per index, by lag | diagnosing 2 |
| `fn.par_bounds_report(best, bounds, base)` | position in bounds, bound hits, far-from-base | diagnosing 2–3 |
| `fn.ga_report(label)` | per-gen monitoring table + stop-rule status + protected metrics | GA monitoring |
| `fn.base_noise(LL)` | replicate SD, per-group SD, GUI z | phase loop step 4 |

GA fixes:
- `timeout` argument through `fn.GA` → `safe_runEwE` (default to the sensitivity's 1800 s).
- Retry the same parvec once before a fresh random draw (covers timeouts and first-wave
  failures); log the replacement count per gen.
- Local mutation using `mutate.margin`.
- `ga_results` header misalignment when `myconfig` has NULL entries.
- Record `myconfig` in `run_info.csv`.

Templates: generic setup / sensitivity / GA / noise scripts derived from the WFS MICE
calibration scripts, with model-specific values pulled into the USER INPUT blocks.
