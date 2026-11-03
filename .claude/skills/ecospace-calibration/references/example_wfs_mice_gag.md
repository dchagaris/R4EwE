# Worked example: keeping a stock decline (WFS MICE gag, 2026)

Shows both playbooks on one problem. The numbers are this model's; the sequence of
questions is the transferable part.

## The problem
Protected behaviour: the assessed decline of adult gag (metric: gag 5+ biomass end/start
ratio; base ×0.37, assessment ×0.17). The phase-5 GA cut LL 39.6%, but its best ended at
×1.13: the decline was gone.

## Diagnosing playbook on the phase-5 best
1. **Attribute.** ~90% of the gain came from five age-0 regional indices (weights 5.5–6.4,
   from index CVs). Gag F series (weight 1) contributed −6.5 LL. Assessment biomass series
   had weight 0.
2. **Mechanism.** Residual check: observed catch −3.2%/yr and F +1.7%/yr imply biomass
   −4.9%/yr. Base matched F and overpredicted catch +43%. The phase-5 best matched catch
   with F 36% low - fleet parameters had cut F (symptom table: "catch matched, F too low").
3. **Flat run** beat base: the age-0 indices rewarded smoothness (SD ratio < 1,
   correlation ≈ 0 - damping, not fitting).
4. **Verdict:** an escape; the weights did not let any series see the decline.

## Reweighting playbook (phase 6)
1. **Breakers:** the phase-5 best, a flat run, and the top sensitivity run (one red-tide
   response off, gag 5+ ×1.06).
2. **First attempt failed, and why:** adult F weight 7 left the ranking almost unchanged.
   F follows effort in Ecospace (0.24–0.27 in every run), and one "adult" F series mapped
   to a non-lead stanza, so it was relative. Lesson: check what each series can see before
   weighting it.
3. **Which series see it:** at unit weights, catch penalised the driver-off breaker and
   absolute F penalised the F-cutting one. Neither alone covered both; together they pin
   biomass (B = C/F).
4. **Scan** (589 runs + breakers, from one unit-weight scoring):

   catch | F (lead) | biomass | breaker ΔLL | old GA best ΔLL | flat ΔLL | ρ(ΔLL, ratio)
   ---|---|---|---|---|---|---
   3 | 7 | 1 | −2568 | −3239 | −2498 | −0.11
   15 | 20 | 5 | −624 | −1457 | −864 | +0.17
   20 | 20 | 5 | +888 | −1275 | −204 | +0.27

   Chosen: the lightest row at which the breaker stops gaining. CV check: catch weight 20 =
   CV 5% (landings, defensible); F 20 for an assessment-derived F is tight, accepted to pin
   the trend. The age-0 index weights were left alone (data CVs).
5. **Verify:** rescored sensitivity - runs raising gag 5+ >10% that improved fell from 7
   of 12 to 1. Noise SD 15.9 → 21.2; stop tolerance reset.

## Result (phase-6 GA, 30 gens)
- LL −41.8% with gag 5+ at ×0.37 and catch and F each within 4% of observed: the decline
  was kept.
- **Next escape found** (diagnosing step 3): the only 4 of 98 parameters at bounds were
  the inflection and slope of two adult red-tide mortality responses, both at the
  "weakest" end - the GA reduced a mortality driver instead of the fishery. That became a
  question for the user (bounds, or the driver's data), not something to fix by weights.
- Age-0 timing partly matched (correlation up to 0.57 from ≈ 0), but variability still
  damped (SD ratio 0.31–0.62), as the likelihood maths predicts.

## Transferable points
Decompose before explaining. Ask which series can see the behaviour, and whether they are
absolute and free to move. Test weights on the breaking runs before launching. Translate
weights to CVs. After the fix, look for the next escape.
