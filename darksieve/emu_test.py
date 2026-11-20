"""Test the emulator against full nonlinear blind searches, fake by fake.
For each blind null fake: its true lens error (m1 / m3 / m4 amplitudes) is known -> emulator predicts
P(blind search returns a candidate with d > T) by Monte Carlo over noise. One number (evidence offset c) is fitted on
even-numbered fakes only; scored on odd-numbered fakes (log-loss, AUC, predicted vs observed rates).
usage: emu_test.py emu.npz blind.jsonl..."""
import sys, json, numpy as np
from sklearn.metrics import roc_auc_score
Z=np.load(sys.argv[1]); G=Z['G']; half=Z['half']; keys=list(Z['keys']); NN=len(G)
Gc=(G-G.mean(0))/np.sqrt(NN-1)
rng=np.random.default_rng(1); ND=4000; Zs=rng.standard_normal((ND,NN)); noise=Zs@Gc       # ND emulated noise realizations, template space
def err_vec(r):
    v=np.zeros_like(half)
    if r.get('m1'): v+=Z['E_m1c']*r['m1'][0]/0.01+Z['E_m1s']*r['m1'][1]/0.01
    if r.get('m34'):
        m,c,s=r['m34']; v+=Z['E_m%dc'%m]*c/0.01+Z['E_m%ds'%m]*s/0.01
    return v
def maxgain(r): return (noise+err_vec(r)-half).max(1)          # ND samples of the blind max statistic
recs=[json.loads(l) for f in sys.argv[2:] for l in open(f) if l.strip()]
recs=[r for r in recs if 'error' not in r and r['kind']!='inj']
det=lambda c: c['mirror']<max(2,.3*c['d']) and c['off']<0.3
obs=lambda r,T: any(det(c) and c['d']>T for c in r.get('cands',[]))
MG={id(r):maxgain(r) for r in recs}
def prob(r,T,c): return np.clip((MG[id(r)]-c>T).mean(),1e-4,1-1e-4)
cal=[r for r in recs if r['k']%2==0]; tst=[r for r in recs if r['k']%2==1]
def ll(rs,c,Ts=(5,8,12)): return -sum(np.log(prob(r,T,c) if obs(r,T) else 1-prob(r,T,c)) for r in rs for T in Ts)
cs=np.arange(-5,60,0.5); best=cs[int(np.argmin([ll(cal,c) for c in cs]))]
print('null fakes: cal %d test %d | fitted evidence offset c = %.1f'%(len(cal),len(tst),best))
for T in (5,8,10,12,15):
    p=np.array([prob(r,T,best) for r in tst]); o=np.array([obs(r,T) for r in tst])
    auc=roc_auc_score(o,p) if 0<o.sum()<len(o) else float('nan')
    print('T=%2d  test: predicted false-accept rate %.3f  observed %d/%d = %.3f  | per-fake AUC %.2f'%(T,p.mean(),o.sum(),len(o),o.mean(),auc))
for k in ('null','null_m1','null_m34','null_knots'):
    rs=[r for r in tst if r['kind']==k]
    if rs: print('  %-10s T=8 predicted %.2f observed %d/%d'%(k,np.mean([prob(r,8,best) for r in rs]),sum(obs(r,8) for r in rs),len(rs)))
json.dump(dict(c=float(best)),open(sys.argv[1].replace('.npz','_calib.json'),'w'))
