"""Score the DarkSieve benchmark. Thresholds are fixed on the CALIBRATION half (even k) and accuracy is measured on the TEST half (odd k).
Accepted = passes detection gates (d, mirror, offset) AND validation gates (Hide-and-Predict, split, lopsided adversary).
TP = accepted candidate within 0.3" of the injected clump; FP = any other accepted candidate.
usage: bench_eval.py bench.jsonl [MINMASS]"""
import sys, json, numpy as np, collections
R=[]
for l in open(sys.argv[1]):
    try: R.append(json.loads(l))
    except Exception: pass
R=[r for r in R if 'error' not in r]
MINM=float(sys.argv[2]) if len(sys.argv)>2 else 9.0
def stats(c):
    return dict(d=c['d'],mir=c['mirror'],off=c['off'],hp=c.get('hp',{}).get('score',np.nan),sp=c.get('split',{}).get('min',np.nan),adv=c.get('adv',{}).get('d_m1',np.nan))
def is_true(r,c): return r['kind']=='inj' and np.hypot(c['clump'][0]-r['clump'][0],c['clump'][1]-r['clump'][1])<0.3
def accepted(c,T,use):
    s=stats(c)
    if not (s['d']>T['d'] and s['mir']<T['mir'] and s['off']<0.3): return False
    for k in use:
        if not (s[k]>T[k]): return False
    return True
cal=[r for r in R if r['k']%2==0]; tst=[r for r in R if r['k']%2==1]
# calibration: candidates that are NOT true (false candidates) and true ones
det_ok=lambda s:(s['d']>5.0 and s['mir']<2.0 and s['off']<0.3)
falsec=[s for s in (stats(c) for r in cal for c in r.get('cands',[]) if not is_true(r,c)) if det_ok(s)]
truec=[stats(c) for r in cal for c in r.get('cands',[]) if is_true(r,c)]
print('records',len(R),'cal',len(cal),'test',len(tst),'| cal false cands',len(falsec),'true cands',len(truec))
def thr(key,lo_is_bad=True,margin=0.5):
    """smallest threshold above all calibration false candidates that reached validation (if none: 0)"""
    v=[f[key] for f in falsec if np.isfinite(f[key])]
    return (max(v)+margin) if v else 0.0
T=dict(d=5.0,mir=2.0)
for k in ('hp','sp','adv'): T[k]=thr(k)
# detection threshold d: lowest value with no calibration false accept (given validation)
print('validation thresholds from calibration:',{k:round(v,2) for k,v in T.items()})
json.dump(T,open(sys.argv[1].replace('.jsonl','_thresholds.json'),'w'),indent=1)
def score(rs,T,use,label):
    tp=fp=0; nul=0; nfp_fakes=0; comp=collections.defaultdict(lambda:[0,0])
    for r in rs:
        acc=[c for c in r.get('cands',[]) if accepted(c,T,use)]
        hit=any(is_true(r,c) for c in acc); f=sum(1 for c in acc if not is_true(r,c))
        fp+=f; nfp_fakes+=(f>0)
        if r['kind']=='inj':
            comp[r['clump'][2]][1]+=1; comp[r['clump'][2]][0]+=hit; tp+=hit
    n=len(rs); F=nfp_fakes/max(n,1)
    above=[(v[0],v[1]) for m,v in comp.items() if m>=MINM]; ca=sum(a for a,b in above)/max(sum(b for a,b in above),1)
    pur=tp/max(tp+fp,1)
    print('%-28s fakes %3d | TP %3d FP %3d | purity %.2f | completeness(M>=%.2f) %.2f | per-mass %s | false-alarm rate per fake %.3f'%(label,n,tp,fp,pur,MINM,ca,{m:'%d/%d'%tuple(v) for m,v in sorted(comp.items())},F))
    return pur,ca,F
for use,label in (([],'detection only'),(['hp'],'+ Hide-and-Predict'),(['hp','sp'],'+ split test'),(['hp','sp','adv'],'+ lopsided adversary (full)')):
    print('CAL ',end=''); score(cal,T,use,label)
    print('TEST',end=''); p,c,F=score(tst,T,use,label)
print('by fake class (TEST, full pipeline):')
for kind in sorted(set(r['kind'] for r in tst)):
    rs=[r for r in tst if r['kind']==kind]; score(rs,T,['hp','sp','adv'],kind)
