"""Linearized false-positive predictor (prototype).
Uses the pipeline's own pixelized-source solve, linearizes only the smooth-model parameters:
  data_E = noise-free truth image with the omitted error E (here an m=1 multipole)
  r(Q, clumps) = noise-weighted residual after the pixelized source is re-solved on data_E
  J = d r / d theta (smooth-model parameters), computed once on the error-free truth
  chi2_smooth = |perp_J r(Qt)|^2 ;  chi2_sub(p,M) = |perp_J r(Qt, subhalo at p)|^2
  predicted false gain = max_M (chi2_smooth - chi2_sub) / 2     (no noise: the systematic part only)
usage: fp_predict.py CFG records.jsonl..."""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
cfgp=sys.argv[1]; files=sys.argv[2:]
sys.argv=['bench.py',cfgp,'0','0','/tmp/none.jsonl']
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'bench.py')).read(); exec(src[:src.index('rng=np.random.default_rng(7000+wid)')])
from model import MACRO, MP
mask=~np.array(data0.mask); sig=np.array(data0.noise_map.native)[mask]
def truth_img(extra=None):
    t=ds.tracer(L,Qt,[],extra); g=list(t.galaxies)[:-1]+[srcG]
    return np.array(al.FitImaging(dataset=data0,tracer=al.Tracer(galaxies=g)).model_data.native)
def resid(Q,clumps=()):
    f=al.FitImaging(dataset=L.ds,tracer=ds.tracer(L,Q,list(clumps))); return (np.array(L.ds.data.native)-np.array(f.model_data.native))[mask]/sig
T0=truth_img(); L.ds=data0; L.set_data(T0); r0=resid(Qt)
H=dict(te=0.003,slope=0.01,my=0.003,mx=0.003,e0=0.003,e1=0.003,g1=0.003,g2=0.003,a3c=0.002,a3s=0.002,a4c=0.002,a4s=0.002)
keys=[k for k in MACRO+MP if k in Qt and k!='reg']
J=np.array([(resid(dict(Qt,**{k:Qt[k]+H.get(k,0.003)}))-r0)/H.get(k,0.003) for k in keys]).T
Qj,_=np.linalg.qr(J); perp=lambda v: v-Qj@(Qj.T@v)
print('jacobian',J.shape,flush=True)
out=[]; seen=set()
for fn in files:
    for l in open(fn):
        r=json.loads(l)
        if r.get('kind')!='null_m1' or not r.get('cands') or 'm1' not in r or (r['w'],r['k']) in seen: continue
        seen.add((r['w'],r['k'])); c=r['cands'][0]; p=c['start'][:2]
        ex=dict(mpx1=al.mp.PowerLawMultipole(m=1,centre=(Qt['my'],Qt['mx']),einstein_radius=Qt['te'],slope=2.05,multipole_comps=tuple(r['m1'])))
        L.ds=data0; L.set_data(truth_img(ex))
        rs=perp(resid(Qt)); c2s=float(rs@rs); best=-1e9
        for lm in (8.5,9.0,9.5):
            rc=perp(resid(Qt,[(p[0]+0.0013,p[1]+0.0017,lm)])); best=max(best,(c2s-float(rc@rc))/2)
        out.append(dict(w=r['w'],k=r['k'],amp=float(np.hypot(*r['m1'])),chi2_smooth=c2s,pred=best,d=c['d'],mir=c['mirror'],fp=bool(c['d']>5 and c['mirror']<max(2,.3*c['d']))))
        print(json.dumps(out[-1]),flush=True)
        json.dump(out,open('/tmp/fp_pred.json','w'))
