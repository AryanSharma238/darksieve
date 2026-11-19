"""Per-lens null test for one candidate.  usage: null_at.py LENSCFG y x lM nsim seed out.json
Simulates the lens from its own best smooth model + fresh noise, with NO subhalo ('null') and WITH a subhalo
of mass lM at (y,x) ('inj'), alternately, and runs the same joint-fit detection at (y,x). The candidate's real
score is compared to the null distribution (false-alarm probability) and to the injection distribution (could a
real subhalo of this mass even produce it?)."""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
import autolens as al, pj, hst
from qload import load; from model import Lens, MACRO, MP
import ds
from mock import noise
cfg=json.load(open(sys.argv[1])); y,x,lm=map(float,sys.argv[2:5]); nsim=int(sys.argv[5]); seed=int(sys.argv[6]); out=sys.argv[7]
P,Q=load(cfg['start'],cfg['fit']); fast=cfg.get('fast',False)
L=Lens(cfg['name'],cfg['zl'],cfg['zs'],P,mesh=(30,30) if fast else None,os_=1 if fast else None)
d=cfg['det']; F=float(hst.inputs(cfg['name']).get('F',1.0))
Qb,_,vb=ds.nm_fit(lambda Q_,c:ds.post(L,Q_),Q,MACRO+MP,None,d['nf_macro'])
f=al.FitImaging(dataset=L.ds,tracer=ds.tracer(L,Qb)); smooth=np.array(f.model_data.native)
# subhalo image added on top of the smooth truth (lensing perturbation of the reconstructed source)
fc=al.FitImaging(dataset=L.ds,tracer=ds.tracer(L,Qb,[(y+0.0013,x+0.0017,lm)])); withc=np.array(fc.model_data.native)
data0=L.ds; res=dict(name=cfg['name'],pos=[y,x,lm],null=[],inj=[])
for i in range(nsim):
    kind='null' if i%2==0 else 'inj'
    L.ds=data0; L.set_data((smooth if kind=='null' else withc)+noise(L,F,seed*1000+i))
    Q1,_,v1=ds.nm_fit(lambda Q_,c:ds.post(L,Q_),Qb,MACRO+MP,None,d['nf_macro'])
    Q2,_,v2=ds.nm_fit(lambda Q_,c:ds.post(L,Q_),Q1,MACRO+MP,None,d['nf_verify']); v0=max(v1,v2); Q0=Q2 if v2>v1 else Q1
    _,c,vc=ds.nm_fit(lambda Q_,c_:ds.post(L,Q_,[c_]),Q0,MACRO+MP,(y,x,lm),d['nf_verify'])
    res[kind].append(float(vc-v0)); json.dump(res,open(out,'w')); print(kind,round(vc-v0,2),flush=True)
