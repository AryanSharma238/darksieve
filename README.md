# DarkSieve

**Finding and validating dark-matter subhalos in strong gravitational lenses.**

A dark-matter subhalo bends the light of a lensed background galaxy slightly. DarkSieve searches lens images for those
bends, then puts every candidate through validation tests designed so that each known type of false alarm has a test aimed at it.

## Pipeline

| Stage | What it does | Kills which false alarm |
|---|---|---|
| D1 Smooth model | Power-law + shear + m=3,4 multipoles, adaptive pixel source (PyAutoLens) | — |
| D2 Clump scan | Clump starts small at each arc grid point; mass, lens shape and source smoothing re-adjust | — |
| D3 Verification | Joint fit vs equal-budget smooth fit; **mirror control** at (−y,−x); position check | Fit improvements that a clump anywhere would give |
| V1 **Hide-and-Predict** *(new)* | Hide pixels near the candidate, refit on the rest, compare predictions of the hidden pixels with/without the clump | Over-fitting, source-driven detections |
| V2 **Split test** *(new for subhalos)* | Two independent halves of the data must both prefer the clump at the same place | Noise |
| V3 Lopsided adversary | Smooth model gains an m=1 multipole; clump must still win | Lopsided lens galaxies |
| Calibration | Thresholds fixed on one half of the simulated lenses, accuracy measured on the other half | Self-tuning |

Novelty checks (ADS full-text, 19 Nov 2026): no prior use found of held-out-pixel prediction, split-exposure replication,
trials-factor / false-discovery control, or curl null tests for subhalo validation in extended arcs. Details in `NOVELTY.md`.

## Accuracy benchmark

`darksieve/bench.py` builds realistic and adversarial simulated lenses from a real lens (SDSS J0946+1006 by default):
no clump, lopsided lens (m=1, 1–2%), extra boxy/triangular structure (m=3/4, 1–2%), clumpy star-forming source,
and injected clumps of 10^8.5–10^9.5 solar masses. `darksieve/bench_eval.py` scores completeness (per mass) and purity.

Targets: completeness ≥90% (goal 95%) above the measured mass limit, purity ≥90% (goal 95%).

## Running

Locally:

```
pip install -r requirements.txt
DS_DATA=data python3 darksieve/bench.py configs/bench_J0946.json <worker_id> <n_fakes> results/bench.jsonl
python3 darksieve/bench_eval.py results/bench.jsonl 9.0
```

At scale on GitHub Actions: Actions → *darksieve-benchmark* → Run workflow (default: 20 parallel jobs × 4 workers × 5.5 h).
The *score* job merges all results and prints the accuracy table.

Real lens:

```
DS_DATA=data python3 darksieve/run_real.py configs/<lens>.json results/all_<run>_thresholds.json out.json
```

## Status

Pipeline and benchmark are built and run end to end. Accuracy numbers come from the benchmark runs in `results/`.
On a Euclid Q1 lens (E102044825), a 10^10 solar-mass clump on the arc gave only Δlog-evidence ≈ 2 after a full refit,
so individual-subhalo detection in Euclid VIS is limited to the most massive clumps. The accuracy benchmark therefore
uses HST data; Euclid is best suited to population-level statistics.
