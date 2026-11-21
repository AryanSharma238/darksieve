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

## Independent-lens test of the false-positive emulator (2026-11-20) - FAILED to transfer
Pre-registered (commits bac371b, 866d937) with frozen J0946 calibration.
- J0252 (144 blind fakes): per-fake AUC 0.28-0.53 (J0946 held-out was 0.80-0.94). Rates: predicted 22% vs observed 39% at
  d>5; predicted 9% vs observed 2% at d>12. Lens-specific threshold d>16 kept false detections at 1/109 but conservatively.
- J1430 (66 blind fakes): 72% of null lenses give a false detection even at d>12; 1/10 injected subhalos found.
Diagnosis:
- J0252 false detections cluster at one ring position in fake after fake; they survive long fits, fits from the truth
  and white noise. Noise-free fake at the true model: gain 0.2; at the parameters a noisy fit settles on: 4.2 ->
  smooth-model drift into a near-degenerate solution that leaves a local mismatch (not captured by a linear emulator).
- J1430: noise-free fake at the TRUE lens model already gains 4-5 (30x30 source grid) and 12-16 (50x50) at two
  hotspots -> source-regularization bias: a regularized pixel source cannot reproduce sharp source structure and a
  subhalo is preferred. Real-data J1430 candidates are therefore not credible (consistent with Nightingale+2024 FP).
- Hotspot map from half the null fakes did not beat a fixed d>12 bar (J0252 0-5% vs 0-3%; J0946 8-9% vs 3-6%).
Conclusion so far: a fixed d>12 bar gives 0-6% false detections per lens on J0946 and J0252, but fails completely on
J1430. Per-lens null audits are required before trusting any candidate; audits running for J0956, J1416, J1020.

## Asimov audit (noise-free source-regularization check), 14 HST lenses
Noise-free image of each lens's smooth model, searched with the pipeline's own 30x30 pixel source at the true lens.
Validated against full blind benchmarks: J1430 hotspot-prone (benchmark 72% false at d>12), J0252 and J0946 clean (2-6%).
| lens | max noise-free gain | positions with gain > 4 |
|---|---|---|
| J1416+5136 | 16.5 | 20 (around (0.2..0.5, 1.4..1.7)) |
| J1430+4105 | 6.1 | 9 |
| J1525+3327 | 4.7 | 1 (at (-0.1,-1.6)) |
| J1250+0523 | 4.3 | 2 |
| J1205, J0252, J1020, J1218, J0956, J2300, J1627, J0903, J0946, J1029 | -1.5 .. 2.8 | 0 |
4 of 14 lenses are hotspot-prone. Pre-registered (commit 9b27c76): J1416 candidate (0.50,1.72) likely an artefact;
J1020 and J0956 clean. J1525 candidate (-0.08,-1.86) is 0.26" from its hotspot -> suspect.

## Per-lens null audits (J0956, J1416, J1020) vs pre-registered Asimov predictions (commit 9b27c76)
| lens | Asimov prediction | false detections d>5 / d>12 / d>20 | false clusters | completeness >=1e9 at d>5 |
|---|---|---|---|---|
| J1416 | hotspot-prone near (0.2..0.5,1.4..1.7), >20% at d>12 | 35% / 14% / 2% | (0.2,1.7) | 0/9 |
| J0956 | clean | 12% / 5% / 2% | scattered | 2/8 |
| J1020 | clean | 23% / 2% / 0% | scattered | 2/14 |
Qualitatively 3/3 correct (J1416 worst, cluster where predicted); the J1416 rate was overstated (14% vs >20%).
Real candidates against their own lens audit:
- J0956 (-0.81,-0.45) d=115: no null fake comes close (0/42), no null detection within 0.4". Strongest candidate. May be
  the luminous satellite included by Vegetti+2014 or Nightingale+2024's Category A candidate - position match pending.
- J1020 (-1.30,0.80) d=13.3: 1/53 nulls exceed it; nulls give weak detections nearby (max d 6.0). Moderate.
- J1416 (0.50,1.72) d=24.7: above all nulls overall, but 7 null detections within 0.4" up to d=18.6 (confirmed
  hotspot). Suspect.
Blind completeness is low in these lenses (0-2 of 8-14 injected >=1e9 subhalos found).

## External cross-checks (literature)
- Nightingale+2024 report their J1430 false positive at (x,y) ~ (0.2", -1.2") (PyAutoLens convention, as here). Our
  Asimov audit's 2nd-strongest J1430 hotspot is (-1.6, 0.2) [y,x], 0.40" away. Consistent, but J1430 has 9 hotspots
  covering much of the arc, so a chance match this close is plausible (~30-40%); not counted as confirmation.
- Nightingale+2024: in J0956 a line-of-sight galaxy lies inside the Einstein radius; Vegetti+2014 modelled a luminous
  satellite in J0956. Our J0956 candidate (-0.81,-0.45) may be one of these; positions not yet matched.
- 8 more lenses' Asimov predictions pre-registered (commit ca54e08); their null audits are running (seeds 30000-37000).

## 8-lens pre-registered Asimov test (commit ca54e08) - MIXED
Clean predictions: 6/6 correct (false detections at d>12: 0-4%). Hotspot-prone predictions: 0/2 correct (J1525 gain 4.7
and J1250 gain 4.3 gave 2% and 0% at d>12). All 14 lenses: Spearman(Asimov max gain, false rate) = 0.58 (p=0.03) at
d>8, 0.49 (p=0.08) at d>12; pooled prone 21% vs clean 1.8% at d>12, driven by J1430 (gain 6.1, 72%) and J1416 (16.5, 14%).
The cut-off of 4 was too low; a revised rule (>= 6) is pre-registered for a fresh Euclid test (configs/PREREG_asimov_v2_euclid.json).
Post hoc (to be tested): false detections in clean lenses often sit near the Asimov top positions (J2300, J1205, J1020).

## Toy inverse problem with exact ground truth (toy/toy_inverse.py, toy/run_toy.py, toy/run_toy_fair.py)
Linear-Gaussian model (exact Bayesian evidence, no fitter): 16x16 source of smooth + sharp blobs, local radial deflection
("perturber"), Gaussian PSF, white noise (5% of peak), gradient (smoothness) prior with regularization searched.
150 random instances; per instance one noise-free Asimov search and 8 noisy null searches; rules on even, scored on odd.
| held-out | Asimov from TRUE source (oracle) | FAIR Asimov (source estimated from one noisy dataset) | data-only chi2 baseline |
|---|---|---|---|
| AUC, flag instances with >=50% false-detection rate (d>10 / d>12 / d>14) | 0.95 / 0.97 / 0.98 | 0.89 / 0.93 / 0.97 | 0.53 / 0.46 / 0.75 |
| Spearman with false-detection rate | 0.79 / 0.88 / 0.58 | 0.69 / 0.75 / 0.41 | 0.03 / 0.04 / 0.25 |
False detections (d>10) sit at the FAIR Asimov peak position 81% of the time (chance 4%).
Unexpected: sharp source features did NOT drive the bias here; false detections arise from general mismatch between the
smoothness prior and the source shape (lenses with 0 sharp features were as affected).
Caveats: toy is linear and exact; it has no smooth-model drift (the J0252 mechanism) and a high overall false-detection
regime; it has not been compared with a parametric-bootstrap baseline (simulating nulls from the fitted model).

## Euclid pre-registered Asimov test (configs/PREREG_asimov_v2_euclid.json, audit set PREREG_euclid_auditset.json) - FAILED
10 Euclid lenses, 1,800+ null fakes. Rate rule: the only hotspot-prone lens (E102042916, Asimov gain 18.5) gave 3/115
false detections at d>8, vs 25/1735 (1.4%) for clean lenses; one-sided Fisher p=0.25 (registered set alone p=0.39).
Spearman(Asimov gain, false rate) = 0.29 (p=0.42). Location rule: 12/28 false detections near the Asimov top-3
positions vs 14.9 expected by chance (p=0.90). Neither pre-registered rule held on Euclid.
Overall record of the Asimov check on real-pipeline lenses: severe HST hotspots (J1430, J1416) correctly flagged; clean
HST predictions 8/8 correct; marginal HST hotspots (J1525, J1250) 0/2; Euclid 0/1 prone and no location signal.
In the exact-evidence toy it works well (held-out AUC 0.89-0.97). The gap between toy and real pipeline is unexplained.
