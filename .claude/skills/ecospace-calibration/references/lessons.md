# Lessons log

General claims only - true of EwE, R4EwE or the likelihood, useful on any model. Each has
one line of evidence. Status: `tentative` (seen in one model), `confirmed` (second model,
EwE source, the maths, or a dedicated experiment), `retracted` (keep, with the reason).
Confirmed lessons also live in the reference named; this log keeps the evidence.

Format:
```
## YYYY-MM-DD - claim  [status]
Evidence: model, run, key numbers (one line)
Where it lives: file/section
```

---

## 2026-08-24 - Effort relaxation is a constant (1 − w) multiplier  [confirmed]
Evidence: WFS MICE, 5 dedicated runs: map effort exactly 0.5× / 0.1× at w = 0.5 / 0.9; EwE source.
Where it lives: ewe_mechanics.md.

## 2026-09-01 - Truncated runs score low and become permanent elites  [confirmed]
Evidence: WFS MICE phase-4 GAs: elites with 10–43 of 170 region files.
Where it lives: ewe_mechanics.md, ga.md.

## 2026-10-23 - Unsaved GUI edits cause console/GUI mismatch  [confirmed]
Evidence: WFS MICE: one group at 2.6× GUI from t = 0; after saving, within 1–3%.
Where it lives: SKILL.md step 3, ewe_mechanics.md.

## 2026-10-29 - Large LL gains can come from removing a behaviour  [tentative]
Evidence: WFS MICE phase 5: −39.6% LL with the stock decline gone; a flat run beat base.
Where it lives: SKILL.md decision rules, diagnosing.md.

## 2026-10-30 - Ecospace F series can't discriminate stock trends  [tentative]
Evidence: WFS MICE phase 6: F 0.24–0.27 across 578 runs; F weight 7 left the ranking unchanged.
Where it lives: ewe_mechanics.md, reweighting.md step 3.

## 2026-11-01 - Absolute catch + absolute lead-stanza F pins a trend through a GA  [tentative]
Evidence: WFS MICE phase 6 ga01: decline kept (×0.37) at −41.8% LL; catch and F within 4%.
Where it lives: reweighting.md step 3.

## 2026-10-30 - Weight changes can be tested by rescoring, not rerunning  [confirmed]
Evidence: unit-weight scaling reproduced direct scoring exactly; maths (q independent of weight).
Where it lives: likelihood.md, reweighting.md step 4.

## 2026-11-01 - Lognormal fits damp variability whose timing is wrong  [confirmed]
Evidence: maths (cov term); WFS MICE sens: every improving run damped age-0 indices.
Where it lives: likelihood.md, reweighting.md "What doesn't work".

## 2026-11-01 - Closing one escape exposes the next  [tentative]
Evidence: WFS MICE phase 6: with F and catch pinned, the only bound hits were adult mortality-driver responses.
Where it lives: SKILL.md decision rules, diagnosing.md step 3.

## 2026-10-30 - Noise ≈ 0.2% of LL, but changes with weights  [tentative]
Evidence: WFS MICE: SD 15.9 → 21.2 after reweighting, both ≈ 0.19% of base LL.
Where it lives: likelihood.md.

## 2026-11-01 - GA console timeout (300 s) selects against slow runs  [tentative]
Evidence: WFS MICE phase 6 ga01: ≥730 fresh random replacements, 100–144 per gen late in the run.
Where it lives: ga.md failure modes.

## 2026-10-30 - First-wave console failures lose the base individual  [tentative]
Evidence: WFS MICE phase 5 gen 1 and phase 6 gen 0 (i0001–i0024 no output).
Where it lives: ga.md failure modes.
