"""DarkSieve on a REAL lens: detection -> validation -> verdict, with thresholds frozen from the benchmark.
usage: run_real.py LENSCFG.json THRESHOLDS.json out.json [halfA.npy halfB.npy]
LENSCFG: {name, zl, zs, start, fit, base_clumps?, fast?, det{...}, R, nf_hp, nf_adv}
Split test runs only if two independent half-data images (from separate exposures) are given."""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from qload import load; from model import Lens
import ds
cfg=json.load(open(sys.argv[1])); T=json.load(open(sys.argv[2])); out=sys.argv[3]
halves=[np.load(sys.argv[4]),np.load(sys.argv[5])] if len(sys.argv)>5 else None
P,Q=load(cfg['start'],cfg['fit']); fast=cfg.get('fast',False)
L=Lens(cfg['name'],cfg['zl'],cfg['zs'],P,mesh=(30,30) if fast else None,os_=1 if fast else None)
ds.BASE=[tuple(c) for c in cfg.get('base_clumps',[])]
det=ds.detect(L,Q,cfg['det']); res=dict(name=cfg['name'],v0=det['v0'],screen_top=det['screen_top'],cands=[])
for c in det['cands']:
    cc={k:c[k] for k in ('screen','start','clump','d','mirror','off','peeled') if k in c}
    gates=dict(detect=c['d']>T['d'],mirror=c['mirror']<max(T['mir'],T.get('mir_frac',0.3)*c['d']),position=c['off']<0.3)
    if all(gates.values()):
        cc['hp']=ds.hide_and_predict(L,P,c,cfg['R'],cfg['nf_hp']); gates['hide_and_predict']=cc['hp']['score']>T['hp']
        cc['adv']=ds.lopsided_adversary(L,P,c,cfg['nf_adv']); gates['lopsided_adversary']=cc['adv']['d_m1']>T['adv']
        if halves is not None:
            cc['split']=ds.split_test(L,P,c,halves); gates['split']=cc['split']['min']>T['sp']
        else: gates['split']='not run (needs separate exposures)'
    cc['gates']=gates
    cc['verdict']='VERIFIED CANDIDATE' if all(v is True for k,v in gates.items() if k!='split') and gates.get('split') is not False else 'rejected'
    res['cands'].append(cc); print('CAND',np.round(c['clump'],3).tolist(),cc['verdict'],gates,flush=True)
    json.dump(res,open(out,"w"),indent=1,default=float)
json.dump(res,open(out,'w'),indent=1,default=float)
