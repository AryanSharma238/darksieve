"""Toy experiment: does the noise-free (Asimov) search predict which instances give false detections?
For each random instance: Asimov map on noise-free data; K noisy null searches; smooth-fit chi2 on one noisy dataset
(a data-only baseline an analyst could compute). usage: run_toy.py first last out.jsonl"""
import sys, json, numpy as np, toy_inverse as T
lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
SIG = 0.05; K = 8
pos = [(py, px) for py in range(3, 13, 2) for px in range(3, 13, 2)]
Wc = {(i, a): T.warp_matrix(p, a) for i, p in enumerate(pos) for a in T.AMPS}
for inst in range(lo, hi):
    rng = np.random.default_rng(1000 + inst)
    s = T.random_source(rng, n_smooth=3, n_sharp=int(inst % 4), sharp_w=0.55)
    y0 = T.B @ s
    asim = T.search(y0, SIG, pos, Wc)
    nulls = []
    for k in range(K):
        y = y0 + rng.normal(0, SIG, len(y0)); d = T.search(y, SIG, pos, Wc)
        nulls.append(dict(max=float(d.max()), argmax=int(d.argmax())))
        if k == 0:   # data-only baseline: residual chi2 of the best smooth reconstruction (MAP) on this noisy dataset
            best = None
            for lam in T.LAMS:
                P = lam * T.DTD + T.EPS * T.EYE; A = P + T.B.T @ T.B / SIG**2; m = np.linalg.solve(A, T.B.T @ y / SIG**2)
                c2 = float(((y - T.B @ m) / SIG) ** 2 @ np.ones(len(y)))
                best = c2 if best is None else min(best, c2)
            chi2 = best
    rec = dict(inst=inst, n_sharp=inst % 4, asimov=asim.tolist(), asimov_max=float(asim.max()), asimov_argmax=int(asim.argmax()), nulls=nulls, chi2=chi2)
    open(out, 'a').write(json.dumps(rec) + '\n'); print(inst, 'asimov %.1f' % asim.max(), 'null max', [round(n['max'], 1) for n in nulls], flush=True)
