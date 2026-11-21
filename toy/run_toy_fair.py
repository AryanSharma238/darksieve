"""FAIR Asimov: the analyst never sees the true source. From ONE noisy dataset, estimate the source with weak
regularization (lam=0.05, keeps sharp structure), make the noise-free image of that estimate, and run the Asimov search
on it. Same instances and noise seeds as run_toy.py. usage: run_toy_fair.py first last out.jsonl"""
import sys, json, numpy as np, toy_inverse as T
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]; SIG = 0.05
pos = [(py, px) for py in range(3, 13, 2) for px in range(3, 13, 2)]
Wc = {(i, a): T.warp_matrix(p, a) for i, p in enumerate(pos) for a in T.AMPS}
for inst in range(lo, hi):
    rng = np.random.default_rng(1000 + inst); s = T.random_source(rng, n_smooth=3, n_sharp=int(inst % 4), sharp_w=0.55)
    y0 = T.B @ s; y = y0 + rng.normal(0, SIG, len(y0))                     # = the first noisy null dataset of run_toy.py
    A = 0.05 * T.DTD + T.EPS * T.EYE + T.B.T @ T.B / SIG ** 2; m = np.linalg.solve(A, T.B.T @ y / SIG ** 2)
    a = T.search(T.B @ m, SIG, pos, Wc)
    open(out, 'a').write(json.dumps(dict(inst=inst, fair_max=float(a.max()), fair_argmax=int(a.argmax()))) + '\n'); print(inst, round(float(a.max()), 1), flush=True)
