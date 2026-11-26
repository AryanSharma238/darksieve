# Full CHIME/FRB Catalog 2 echo search: results

Rules and amendments: RULES.md. Numbers: out_full/summary.json. Per-burst scores and injections: out_full/scores.json.gz.

## Outcome
- 4536 spectrum files processed. 3476/3478 non-repeating and 981/981 repeat bursts scored.
- Old scores reproduced: 322/322.
- Threshold 23 (largest repeat-burst score 22.93). 4 bursts exceed it.
- **No constraint on f at any tested mass (all 95% limits are above 1).**

| Mass (Msun) | locked 95% limit on f (4 candidates counted) | after-the-fact, if all 4 are removed (not a valid headline) |
|---|---|---|
| 100 | 81.9 | 26.8 |
| 300 | 8.2 | 2.7 |
| 1000 | 5.9 | 1.9 |
| 3000 | 12.7 | 4.1 |

## The four bursts above threshold
FRB20200306B, FRB20220219A, FRB20220720B, FRB20211006D. Each has one catalog component, no visible second peak at the fitted delay
(out_full/candidates.png), and a catalog scattering time of 43 to 560 ms. A post-hoc per-band test (vet_candidates.py) finds the
echo ratio changes with frequency (p about 2e-4 or smaller with 16 bands; approximate errors). Interpretation: the detector compared the start of a long tail
with its later part. They are not evidence of lensing. These checks were chosen after seeing them.
23 of the 26 non-repeating bursts with S>15 have a delay within 2 samples of the smallest allowed, which suggests the same weakness more widely
(not individually inspected).

## Problems found in the earlier 340-burst result (withdraw these)
1. Its optical depth counted 3,564 sources, but only 322 bursts were scored. Using only the scored sources, the old 0.93 limit at 1000 Msun becomes about 19.
2. Planted echoes copied the host burst's noise. On six bursts: recovery 16.1% (old) vs 1.8% (noise-free model). Whole catalog: 21.8% vs 5.8%.
   So the 17.8% efficiency and the 3.5-11x advantage over the published method must not be quoted. The published method has not been re-tested here.
3. The top flux-ratio bin applied R=5 results to all R>4; the integration overcounted by about 3%.

## What this does not show
It does not show the Zhou et al. candidates are false. FRB20190131D (S=18) is below the threshold; 5 of 981 repeat bursts score as high or higher.
Not searched: echoes overlapping the main pulse, delays above about 69 ms. 79 files have no catalog model (2 in the search set).
Model-only injections may be somewhat pessimistic. The DM-to-distance assumption changes the optical depth by roughly -22% to +20%
in the alternatives a reviewer tested. A limited web search by an assistant found no earlier full-catalog search of this kind.

## AI use
Code, review and this report were produced with AI assistance (Claude). The earlier code was also reviewed and corrected by AI agents.
