# GA: setup, stopping, monitoring, failure modes

## Parameter selection (from the sensitivity)
- Rank sensitivity runs by |ΔLL| vs base; keep those ≥ `min.dLL` (≈ 2 noise SD), capped at
  `n_parsens` = number of scored data sources (weight > 0 series + spatial maps).
- Keep parameters the user wants searched regardless (e.g. drivers of a protected behaviour).
- Deduplicate by (prey, pred, tag.type, fleet, fxn) - the low and high test of one pair are
  one parameter.
- Sensitivity replicate rows (`tag.type == 'base replicate'`) must not enter selection.
- After choosing, dry-run `fn.makeparvec` and check: parameter count, groups, that every
  hand-set `L.bounds[...]`/`U.bounds[...]` name exists (a missing name silently appends an
  element), and no NA.

## Configuration (fn.GA / myconfig)
- popSize = 5 × number of parameter groups (a response shape counts as one group).
- elitism ≈ 5% of popSize; tournament selection, k = 3; pmutation ~0.03.
- Bounds: vuls relative (`vul.mult` × base), log-uniform draws; env shape multipliers
  [0.5, 2] (can't flatten a response - known limitation); env locations bounded in curve
  widths (`env.loc.sd`); red-tide slope bounds per response (steep slopes drive
  catastrophic fits).
- Plot the initial population against its bounds before launching
  (`fn.plot_ga_pop_responses`; linear axes on vul histograms).
- Stop rule: from `stop.min_gen`, stop when best improves < `stop.rtol × |best|` over
  `stop.patience` gens, or population SD < that tolerance. Set rtol so the tolerance is
  1–2 noise SD at the expected best LL. `n_gen` is a hard cap.
- Resume a dead run with `myconfig$resume.gen` = last completed gen, same label and script.

## Run layout
`run_dir/gen<NN>/` per generation; elites move forward with their names unchanged, the
previous gen folder is deleted. Elite fitness is the original evaluation (no reruns -
reruns are a fresh stochastic draw). `ga_evals_<label>.csv` has every evaluated individual
(gen, idx, role, kept, fitness, run, parameters) - use it to attribute bad fits to
parameter ranges.

## Monitoring (per generation)
Report: best, Δ vs previous and vs base, median, mean, max, SD, NA count, % of offspring
below base, offspring median, new bests, stop-rule status (gain over the patience window vs
tolerance), pace/ETA, and the protected behaviour of the current best (e.g. its trend
ratio). Interpret (convergence, diversity, anomalies) and give a verdict. Push-notify on
each new generation only.

Watch for:
- a **sudden large LL drop** → suspect a truncated run (only Region_0 files) scored low;
- the `base-in-pop` row not matching base → the first wave failed and was replaced by fresh
  random draws (see below);
- a stall (no new run output for ~45 min) → hung console processes.

## Known failure modes
- **First wave fails.** When all workers start the console at once, the whole first batch
  can fail with no output; the retry loop substitutes **fresh random parvecs**, so row 1
  (the base individual) is lost. Harmless for the search, but the in-population base
  reference is gone.
- **Console timeout.** `safe_runEwE` calls `fn.runEwE` with its default 300 s timeout
  (the sensitivity uses 1800 s). Slow parameter sets fail and are replaced by fresh random
  draws - the GA selects against slow runs and wastes evaluations (hundreds per gen late in
  a run). Set the timeout from the slowest sensitivity run; count replacements per gen.
- **Mutation is global.** `mutate()` draws over the full bounds; `mutate.margin` is unused.
  Late generations therefore produce a tail of poor offspring; this is expected, not a sign
  of a problem.
- **Stop rule and the cap.** The stall test compares best[g − patience] with best[g], so
  one improvement larger than the tolerance anywhere in the window prevents a stop; a run
  that keeps making occasional jumps creeps to `n_gen`.
- **Cluster exports.** Workers need R4EwE internals via `fn.cluster_export(cl)` (namespace
  export); a bare `clusterExport` of package names fails under `library()`. Test worker
  changes with the package installed.
- **ga_results header.** A `NULL` entry in `myconfig` (e.g. `resume.gen = NULL`) shifts the
  name=value pairs in the results header; cosmetic, the run uses `myconfig` directly.
- **Truncated output** - see ewe_mechanics.md; requires `expected_regions` exported.
