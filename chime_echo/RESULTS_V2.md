# Detector version 2: results (not blind; see RULES.md)

Winner by the locked development-half rule: hp9 (removes trends smoother than a 9-sample running median).

| | v1 | hp9 |
|---|---|---|
| Planted-lens recovery, held-out half (thresholds from development repeats) | 6.0% | 16.5% |
| Planted-lens recovery, full search | 5.9% | 13.9% |
| Threshold (highest repeat-burst score) | 22.93 | 8.00 |
| Held-out repeat bursts above development threshold | 0/428 | 1/428 |
| Non-repeating bursts above full threshold | 4 | 4 |
| 95% limit on f at 100/300/1000/3000 Msun | 80/8.2/5.8/12.6 | 12.8/2.0/2.0/6.3 |

No constraint at any mass (all limits above 1). Scores of different versions are not on the same scale.

## Bursts above the hp9 threshold
If one-off and repeat bursts behaved alike, about 3.5 of 3,476 would exceed the maximum of 981 repeats (Poisson chance of 4 or more: 0.47).
That assumption is imperfect: at catalog S/N>=50 the repeat maximum is only 6.1, so bright sharp bursts are poorly represented.

Checks chosen after seeing the candidates:
- FRB20210113C, FRB20190805E, FRB20200302C: catalog widths of the two components differ by 2.3 to 3.0 times
  (several times the quoted errors). A point lens should give equal widths. Strong evidence against, assuming the catalog fits are right.
- FRB20190131D (S=16.9, 8.85 ms; ranks first): equal catalog widths (0.273 ms). Catalog spectral index and running differ by about 2.4 sigma.
  Frequency-ratio test, calibrated on 129 planted model copies: 1.5% (4 bands), 5.4% (8), 36.9% (16) look at least as uneven.
  Inconclusive. Planted copies are model-based, and the band choice was not fixed in advance.
  This burst was already published as a lens candidate (arXiv:2605.19653). It is not a new object and this run does not confirm it.
If all four were lenses, f would be about 0.9, which conflicts with other published bounds.

AI use: code and report produced with AI assistance (Claude).
