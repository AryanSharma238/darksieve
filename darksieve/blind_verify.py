"""Stage 2 of the blind search: merge all screen slices for a lens, then verify + validate its peaks.
usage: blind_verify.py LENSCFG SLICE_DIR THRESHOLDS out.json"""
import sys, os, json, glob, subprocess
cfg=json.load(open(sys.argv[1])); sl=sorted(glob.glob(os.path.join(sys.argv[2],'**','*.json'),recursive=True))
scores=[]; done=0; tot=None
for f in sl:
    try: r=json.load(open(f))
    except Exception: continue
    if not isinstance(r,dict) or r.get('name')!=cfg['name']: continue
    scores+=r.get('scores',[]); done+=bool(r.get('done')); tot=r.get('n_total',tot)
nsl=sum(1 for f in sl if os.path.basename(f).startswith('scr_'))
print('slices','complete',done,'points screened',len(scores),'of',tot,flush=True)
cfg['det']['screen_given']=scores; cfg['blind']=dict(slices=done,complete=done,screened=len(scores),total=tot)
json.dump(cfg,open('merged_cfg.json','w'))
r=subprocess.run([sys.executable,'-u',os.path.join(os.path.dirname(os.path.abspath(__file__)),'run_real.py'),'merged_cfg.json',sys.argv[3],sys.argv[4]])
try:
    o=json.load(open(sys.argv[4]))
except Exception: o=dict(name=cfg['name'],cands=[],error='verification did not finish')
o['blind']=cfg['blind']; o['screen_all']=sorted(scores,reverse=True); json.dump(o,open(sys.argv[4],'w'),indent=1,default=float)
