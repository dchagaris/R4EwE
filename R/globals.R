# This suppresses “no visible binding for global variable …” NOTES.
# It’s the accepted mechanism when NSE or global state is intentional. 
# Over time, consider refactoring the most important ones (e.g., pass pop_size, n_pars, cl 
# as explicit function arguments) for robustness.
# R/globals.R
utils::globalVariables(c(
  "L.bounds", "U.bounds", "cl", "depth", "depth.5min", "df.names", "dir.cefi",
  "dir.glorys", "dir.main", "dir.pred", "disp.par.idx", "disp_pars", "do.penalty",
  "elitism", "env.par.idx", "est_par_vec", "file.console", "files.cmd",
  "fleet.num", "gapop.dist", "group.names", "group.num", "i", "lk.ts",
  "makeSOCKcluster", "map", "med.par.idx", "month",
  "mutation_rate", "myconfig", "n_generations", "n_pars", "obs.map", "obs.ts",
  "redtide.par.idx", "selection.method", "tournament.k",
  "par", "par.groups", "par.labels", "par_cv_vec", "pars1", "pars2", "pars3",
  "pars4", "pen.wt.mult", "pop_size", "predprey_pars", "res", "respfxn_num",
  "respfxn_type", "run_dir", "styear", "timestamp", "type", "values", "viridis",
  "vul.par.idx", "vul_pars", "winProgressBar", "year"
))
