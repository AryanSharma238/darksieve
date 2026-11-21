"""Score the toy experiment. Prediction target per instance: false-detection rate (fraction of K noisy null searches with
max d > T). Predictors: Asimov max gain (one noise-free search) vs smooth-fit chi2 on one noisy dataset (data-only baseline).
Rules are fitted on even instances, scored on odd instances."""
import json, sys, numpy as np
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score
R = [json.loads(l) for f in sys.argv[1:] for l in open(f)]
print('instances', len(R))
for T in (8, 10, 12):
    fr = np.array([np.mean([n['max'] > T for n in r['nulls']]) for r in R]); A = np.array([r['asimov_max'] for r in R]); C = np.array([r['chi2'] for r in R])
    odd = np.array([r['inst'] % 2 == 1 for r in R]); bad = fr >= 0.5
    print('T=%d: mean false rate %.2f | Spearman(Asimov, rate) %.2f | Spearman(chi2, rate) %.2f' % (T, fr.mean(), spearmanr(A, fr)[0], spearmanr(C, fr)[0]))
    if 0 < bad[odd].sum() < odd.sum():
        print('      held-out AUC for flagging instances with >=50%% false rate: Asimov %.2f, chi2 baseline %.2f' % (roc_auc_score(bad[odd], A[odd]), roc_auc_score(bad[odd], C[odd])))
# location: does the noisy false detection sit at the Asimov peak?
hit = [n['argmax'] == r['asimov_argmax'] for r in R for n in r['nulls'] if n['max'] > 10]
print('false detections (d>10) located at the Asimov peak position: %d/%d = %.0f%% (chance = 1/25 = 4%%)' % (sum(hit), len(hit), 100 * np.mean(hit) if hit else 0))
