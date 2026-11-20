"""Position-dependent thresholds from null simulations ("hotspot map").
Map: for each 0.3" cell, the null-simulation candidate scores d found there (calibration fakes only).
Rule: accept a candidate at p if d > max(T0, q(cell(p))), q = the highest null score seen in that cell (+margin).
Scored on held-out null and injected fakes.  usage: hotspot.py blind.jsonl [T0] [margin]"""
import sys, json, numpy as np, collections
R=[json.loads(l) for l in open(sys.argv[1]) if l.strip()]; R=[r for r in R if 'error' not in r]
T0=float(sys.argv[2]) if len(sys.argv)>2 else 5.0; mg=float(sys.argv[3]) if len(sys.argv)>3 else 1.0
det=lambda c: c['mirror']<max(2,.3*c['d']) and c['off']<0.3
cell=lambda y,x: (int(np.floor(y/0.3)),int(np.floor(x/0.3)))
nbr=lambda k: [(k[0]+i,k[1]+j) for i in (-1,0,1) for j in (-1,0,1)]
hit=lambda r,c: r['kind']=='inj' and np.hypot(c['clump'][0]-r['clump'][0],c['clump'][1]-r['clump'][1])<0.3
def evaluate(cal,tst,label):
    M=collections.defaultdict(float)
    for r in cal:
        for c in r.get('cands',[]):
            if not hit(r,c) and det(c): k=cell(*c['clump'][:2]); M[k]=max(M[k],c['d'])
    thr=lambda c: max(T0,max(M.get(k,0) for k in nbr(cell(*c['clump'][:2])))+mg)
    out=[]
    for nm,rule in (('fixed d>%g'%T0,lambda c: det(c) and c['d']>T0),('fixed d>12',lambda c: det(c) and c['d']>12),('hotspot map',lambda c: det(c) and c['d']>thr(c))):
        nl=[r for r in tst if r['kind']!='inj']; f=sum(any(rule(c) for c in r.get('cands',[])) for r in nl)
        inj=[r for r in tst if r['kind']=='inj' and r['clump'][2]>=9]; tp=sum(any(rule(c) and hit(r,c) for c in r.get('cands',[])) for r in inj)
        out.append('%-12s null lenses w/ false detection %2d/%2d = %3.0f%% | completeness>=1e9 %2d/%2d'%(nm,f,len(nl),100*f/max(1,len(nl)),tp,len(inj)))
    print(label); [print('   ',o) for o in out]; print('    hotspot cells:',sorted(((round(v,1),k) for k,v in M.items()),reverse=True)[:5])
evaluate([r for r in R if r['k']%2==0],[r for r in R if r['k']%2==1],'map from EVEN fakes, tested on ODD fakes')
evaluate([r for r in R if r['k']%2==1],[r for r in R if r['k']%2==0],'map from ODD fakes, tested on EVEN fakes')
