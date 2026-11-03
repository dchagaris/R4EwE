# EwE / Ecospace mechanics that affect calibration

Only items proven from EwE source or a dedicated experiment are here. Source line numbers
refer to the Feb 2025 EwE source snapshot; build/source match is not formally verified.

## Fishing mortality follows effort, not biomass
Predicted F (catch / biomass) is set by the effort timeseries and the spatial overlap of
effort and biomass. In one model's 578-run sensitivity, adult F of the key group moved only
0.24–0.27 whatever the stock did. Consequence: **an F series cannot tell a declining stock
from a flat one**; only catch (absolute) together with F pins biomass (B = C / F). A
parameter set can also match catch with a flat stock by lowering F through fleet
parameters - only an absolute, well-weighted F series catches that.

## Effort relaxation is a constant (1 − w) multiplier
`cEcoSpace.vb` `PredictEffortDistributionThreaded` zeroes `EffortSpace` before the relaxation
reads it, so `EffortSpace = (1 − w) · EffPred` every timestep. With the default
"Predicted effort relaxation weight" w = 0.9, map effort is cut 10× (confirmed to 5
significant figures). The no-fish attraction (`NoFishWeight`) is a second effort sink.
`UsePenaltySearch = FALSE` conserves effort and avoids both.

## Habitat capacity
HabCap = base map × habitat-use × **product** of env responses, normalised by its max over
water cells **at t = 0** (later timesteps keep that divisor). A response that is uniformly
low cancels out; only its variation across cells and over time matters. Sigmoid, beta,
exponential responses are rescaled to peak = their scalar parameter; normal peaks at its
max parameter; logistic4 and trapezoid are native 0–1. A sigmoid midpoint at the axis edge is
effectively a slope-only response.

## Initial distribution
Initial biomass per cell = HabCap × nWaterCells / TotHabCap × Ecopath B: total is preserved,
so more restrictive capacity means higher density where fish are, and patchier biomass.
Patchy biomass breaks gravity-model effort allocation (catch can be far below Ecosim at
month 1 with biomass unchanged).

## Console and cmd files
- The console applies the cmd file top to bottom; the **last tag wins**. A cmd base already
  encodes the hand-tuned baseline; appending a full parameter block silently overrides it. The
  GA base run must be the cmd base with only the output dir changed.
- The console writes to the model db on every run (imports the effort ts), so the db mtime is
  not a "saved" signal. Unsaved GUI edits are the usual cause of a console/GUI mismatch.
- Discard forcing (types 10, 11) can't be loaded by the console.
- An Ecospace run that dies mid-write leaves only Region_0 files; without the
  `expected_regions` guard it scores a spuriously **low** LL and becomes a permanent elite.

## Observed catch definitions
Ecospace `Catch` = landings + dead discards. Type 13 obs are **total** discards
(DT = (C − L) / Dmort), so compare them with the total-discard prediction.
