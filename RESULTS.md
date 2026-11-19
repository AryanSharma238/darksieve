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
