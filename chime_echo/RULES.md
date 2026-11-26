# Locked rules: full CHIME/FRB Catalog 2 echo search
Committed before any score outside the earlier 340-burst run existed. The earlier 322 scores
(147 non-repeaters) HAD been seen, so this is blind only for the other bursts.

1. Data: CHIME/FRB Catalog 2 Stokes-I spectra (doi:10.11570/25.0066), files in manifest.txt; catalog table chimefrbcat2.csv, sub_num=0 row.
2. Detector: echo.py unchanged. cmax=1.35, delays 2..70 samples (and longer than the main-pulse window).
3. Reproduction gate: at least 99% of the previously scored bursts must give the same S within 0.01.
4. Calibration set: every burst with a repeater_name. Search set: every burst without one, with snr_fitb>0.
5. Threshold T = max(23.0, largest finite S in the calibration set). T=23 is also reported, labelled secondary.
6. k = number of search-set bursts with S>T. Every one is only a CANDIDATE, not a lens detection.
7. Optical depth: summed ONLY over search-set bursts that were actually scored and injection-tested. Repeating sources are excluded.
8. Efficiency: echoes planted in every scored search-set burst (11 delays x 5 flux ratios). Recovered means |d-d_true|<=1 sample and S>T (chi<1.35 enforced by detector). Same bins as the earlier run; empty bins count as 0.
9. Limit: f95 = N95(k)/SumTau_eff, N95 = Poisson 95% one-sided upper limit (k=0: 3.00), no background subtraction. f95>=1 means no constraint.
10. If more than 5% of search-set files fail, the limit is reported as incomplete. All failures are listed.
11. The result is reported whatever it is.

Changes from the earlier run: (a) optical depth no longer counts unsearched bursts; (b) injections use all searched bursts, not 60 of the 340;
(c) injections use cmax=1.35 like the real search (earlier: 99, then cut on chi); (d) repeaters are used only for calibration.

Known limits: echoes overlapping the main pulse are not searched; delays above about 69 ms are not searched; simplified white-noise model;
single lens mass; approximate DM-redshift relation; the catalog model defines the injected signal region.
