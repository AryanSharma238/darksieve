"""Blind screen of one slice of a lens's arc.  usage: screen_shard.py LENSCFG shard nshard out.json
Fits the smooth macro model (deterministic, so every slice gets the same baseline), then screens every
nshard-th arc position starting at `shard`. Writes [(score,y,x,lM),...] after every point."""
import sys, os, json
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qload import load; from model import Lens, MACRO, MP
import ds, stages as S
cfg=json.load(open(sys.argv[1])); k,K,out=int(sys.argv[2]),int(sys.argv[3]),sys.argv[4]
P,Q=load(cfg['start'],cfg['fit']); fast=cfg.get('fast',False)
L=Lens(cfg['name'],cfg['zl'],cfg['zs'],P,mesh=(30,30) if fast else None,os_=1 if fast else None)
ds.BASE=[tuple(c) for c in cfg.get('base_clumps',[])]
d=cfg['det']; Qb,_,vb=ds.nm_fit(lambda Q,c:ds.post(L,Q),Q,MACRO+MP,None,d['nf_macro'])
pts=S.arc_positions(L,Qb,spacing=d['spacing']); mine=pts[k::K]
res=dict(name=cfg['name'],shard=k,nshard=K,n_total=len(pts),v_base=vb,scores=[])
print('lens',cfg['name'],'points total',len(pts),'this slice',len(mine),flush=True)
for p in mine:
    res['scores'].append(list(ds.screen_point(L,Qb,vb,p,d['nf_screen'])))
    json.dump(res,open(out,'w'),default=float); print('PT',[round(float(v),3) for v in res['scores'][-1]],flush=True)
res['done']=True; json.dump(res,open(out,'w'),default=float)
