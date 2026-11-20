# DarkSieve results (2026-11-19)

## Accuracy on fake lenses (J0946 noise, PSF and source; frozen rule chosen on calibration half)
Rule: d > 5, mirror < max(2, 0.3 d), offset < 0.3", Hide-and-Predict > 0, split-exposure > 0, lopsided adversary > 0.

Targeted mode (verification at the true / a random position), 329 fakes:
| | result | 95% interval |
|---|---|---|
| completeness M >= 1e9 | 60/67 = 89.6% | 80-96% |
| completeness M >= 10^9.25 | 44/47 = 93.6% | 82-99% |
| completeness M >= 10^9.5 | 21/22 = 95.5% | 77-100% |
| purity | 66/68 = 97.1% | 90-100% |
Held-out half only: completeness M >= 1e9 93%, purity 97%. Both false accepts were lopsided (m=1) null lenses (2/48);
0/158 other nulls passed.

Blind mode (full search), first 31 fakes: 4 false candidates (1 lopsided null, 2 echoes of strong real subhalos,
1 weak). Fixed in commit eb541e4: stronger lopsided adversary (m=1 grid scan + polish) and peeling of accepted clumps.
Re-test running (targeted seed 11000, blind seed 12000).

## Correction (after advisor review)
- Held-out only (odd-numbered fakes): completeness M >= 1e9 27/29 = 93.1% (77-99%); false-positive rate per null lens
  1/106 = 0.9% (0.0-5.1%); purity 28/29. All fakes: 60/67, 2/220 = 0.9% (0.1-3.2%), 66/68.
- The three validators change none of these numbers: detection + mirror test did all the work so far.
- The sigma and global-p columns below extrapolate a Gaussian tail from 25-40 noise-only simulations and are WITHDRAWN.
  What the simulations support: no null exceeded the observed score, local p < ~1/(N+1). A tail model will be
  calibrated on the hard-null benchmark. Noise-only tests do not cover lens-model errors.

## Blind search of 29 real lenses (796 arc positions, no earlier candidates used)
Positive control J0946: known subhalo (Vegetti et al. 2010, (-1.08,-0.65)) recovered blind at (-1.03,-0.67), passes all validators.

Candidates passing every gate, with per-lens null test (lens re-simulated from its own smooth model + fresh noise):
| lens | position | d | HP | adversary | noise-only sigma | global p (796 trials) |
|---|---|---|---|---|---|---|
| J0956+5100 | (-0.75,-0.39) | 32.9 | 44.3 | 38.3 | 13.4 | 3e-38 |
| J0946 (known) | (-1.03,-0.67) | 19.3 | 20.0 | 18.0 | 12.5 | 5e-33 |
| J0946 | (0.80,-1.30) | 15.7 | 69.1 | 15.2 | 8.0 | 4e-13 |
| J1430+4105 | (0.18,1.10) | 12.7 | 58.7 | 13.3 | 11.8 | 2e-29 |
| E102020532 | (0.20,0.81) | 8.8 | 5.9 | 8.8 | 7.2 | 2e-10 |
| E102019128 | (-0.10,-0.41) | 10.8 | 4.5 | 13.3 | 4.8 | 5e-4 |
| J1430+4105 | (0.79,1.12) | 5.0 | 8.1 | 3.5 | 5.2 | 8e-5 |
| E102044825 | (0.20,-0.40) | 8.7 | 9.7 | 28.9 | 2.3 | not significant |

Caveats: the null test covers noise only, not model misspecification; sigma values use a Gaussian fit to 25-40
null simulations. J0946 has a second source plane not in the model, which could affect (0.80,-1.30).
Rejected by validators despite a large signal: J0956 (1.14,1.17) d=96.7 (Hide-and-Predict -16.3).

## Update 2026-11-20 (v2 code: stronger adversary, peeling; verdict bug fixed)
Targeted, lopsided-heavy benchmark:
| fit budget | lopsided nulls accepted | other nulls | completeness >= 1e9 |
|---|---|---|---|
| 150 steps (seed 11000) | 3/78 | 0/90 | 28/38 = 74% |
| 400 steps (seed 13000) | 6/56 (detection only 7/56) | 0/54 | 26/32 = 81% |
Longer, equal fits make lopsided false positives MORE common: the m=1 / subhalo confusion is a real evidence-level
degeneracy, not mainly a fitter artefact. First case of a validator removing a false accept (7 -> 6).

Blind benchmark (115 fakes, ~134 positions searched each): at d > 5, 11/71 null lenses falsely accepted, purity ~61%.
d > 12 (chosen on calibration half): held-out 2/35 null lenses falsely accepted, completeness 5/9. The targeted
benchmark (one position per fake) understated false positives by ignoring the trials factor.

Linear false-positive predictor, prototype v1 (152 lopsided nulls, 5 false positives): position-specific predicted gain
AUC 0.36 (fails); leftover chi2 after nuisance absorption AUC 0.74.

Finder, 9 new HST lenses (verdicts recomputed with fixed gates): pass at d > 12: J1020+1122 (-1.30, 0.80) d 13.3,
J1416+5136 (0.50, 1.72) d 24.7. Pass only at d > 5: J1218, J1250, J1416 (-1.57,-0.08), J1525 x2.
