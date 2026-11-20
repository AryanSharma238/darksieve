"""Emulator score -> calibrated false-detection probability.
score(fake) = median over emulated noise of the blind max gain (lens error + noise, no fitting).
P(false detection at threshold T) = sigmoid(a + b*log(1+score) - g*T): 3 numbers, fitted on even-numbered fakes only,
scored on odd-numbered fakes (AUC, Brier vs base rate, predicted vs observed rate).
usage: emu_calib.py emu.npz blind.jsonl... [--freeze out.json] [--use frozen.json]"""
import sys, json, numpy as np
from scipy.optimize import minimize
from sklearn.metrics import roc_auc_score
args=[a for a in sys.argv[1:] if not a.startswith('--')]; Z=np.load(args[0]); G=Z['G']; half=Z['half']; NN=len(G)
frz=sys.argv[sys.argv.index('--freeze')+1] if '--freeze' in sys.argv else None
use=sys.argv[sys.argv.index('--use')+1] if '--use' in sys.argv else None
if frz: args=[a for a in args if a!=frz]
if use: args=[a for a in args if a!=use]
Gc=(G-G.mean(0))/np.sqrt(NN-1); rng=np.random.default_rng(1); noise=rng.standard_normal((2000,NN))@Gc
def err_vec(r):
    v=np.zeros_like(half)
    if r.get('m1'): v=v+Z['E_m1c']*r['m1'][0]/0.01+Z['E_m1s']*r['m1'][1]/0.01
    if r.get('m34'): m,c,s=r['m34']; v=v+Z['E_m%dc'%m]*c/0.01+Z['E_m%ds'%m]*s/0.01
    return v
score=lambda r: float(np.median((noise+err_vec(r)-half).max(1)))
recs=[json.loads(l) for f in args[1:] for l in open(f) if l.strip()]; recs=[r for r in recs if 'error' not in r and r['kind']!='inj']
det=lambda c: c['mirror']<max(2,.3*c['d']) and c['off']<0.3
obs=lambda r,T: any(det(c) and c['d']>T for c in r.get('cands',[]))
S={id(r):np.log1p(max(score(r),0)) for r in recs}; TS=(5,8,10,12,15)
sig=lambda x: 1/(1+np.exp(-x))
def nll(w,rs): return -sum(np.log(np.clip(sig(w[0]+w[1]*S[id(r)]-w[2]*T),1e-6,1-1e-6) if obs(r,T) else np.clip(1-sig(w[0]+w[1]*S[id(r)]-w[2]*T),1e-6,1)) for r in rs for T in TS)
cal=[r for r in recs if r['k']%2==0]; tst=[r for r in recs if r['k']%2==1]
if use: w=json.load(open(use))['w']; tst=recs; print('using frozen weights',w,'on ALL',len(tst),'null fakes')
else: w=minimize(nll,[-2,1,0.1],args=(cal,),method='Nelder-Mead').x; print('fitted on %d calibration null fakes: a=%.2f b=%.2f g=%.3f'%(len(cal),*w))
for T in TS:
    p=np.array([sig(w[0]+w[1]*S[id(r)]-w[2]*T) for r in tst]); o=np.array([obs(r,T) for r in tst])
    auc=roc_auc_score(o,p) if 0<o.sum()<len(o) else float('nan')
    print('T=%2d  predicted rate %.3f  observed %d/%d = %.3f | AUC %.2f | Brier %.3f (base-rate %.3f)'%(T,p.mean(),o.sum(),len(o),o.mean(),auc,np.mean((p-o)**2),np.mean((o.mean()-o)**2)))
if frz: json.dump(dict(w=list(map(float,w))),open(frz,'w')); print('frozen ->',frz)
