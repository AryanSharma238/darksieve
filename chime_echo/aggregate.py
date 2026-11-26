import csv,json,glob,sys,numpy as np
from scipy.stats import chi2
from tau import zm,dt_ms,ymin,ymax_R,Dc,E,Oc,H0,c
DT=np.array([2,3,4,6,8,12,16,24,32,48,64],float); RF=np.array([1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0,25.0,50.0])
def weight_points(M,d,wi,s,nz=20,ny=30):
    """optical-depth samples (delay ms, flux ratio, weight) for one burst; midpoint rule in y"""
    zs=zm(d)
    if not np.isfinite(zs) or not s/10>=1: return []
    ym=ymax_R(s/10); P=[]
    for z in np.linspace(zs/nz/2,zs*(1-1/nz/2),nz):
        Dl=Dc(z)/(1+z); Ds=Dc(zs)/(1+zs); Dls=(Dc(zs)-Dc(z))/(1+zs); y0=ymin(z,M,wi)
        if ym<=y0: continue
        pref=1.5*Oc*H0**2/(c*H0*E(z))*Dls*Dl/Ds*(1+z)**2*(zs/nz)
        dy=(ym-y0)/ny
        for y in y0+dy*(np.arange(ny)+0.5):
            R=(y*y+2+y*np.sqrt(y*y+4))/(y*y+2-y*np.sqrt(y*y+4))
            P.append((dt_ms(y,M,z),R,pref*2*y*dy))
    return P

def burst_tau(P,hit):
    """hit: 11x10 array of this burst's own injection outcomes; nearest grid point in log delay and log R.
    Delays outside the searched range (below 1.5 or above 69.5 samples) and R>60 get zero."""
    if not P: return 0.,0.
    t,R,w=np.array(P).T
    j=np.abs(np.log(t)[:,None]-np.log(DT)).argmin(1); k=np.abs(np.log(R)[:,None]-np.log(RF)).argmin(1)
    ok=(t>=1.5*0.983)&(t<=69.5*0.983)&(R<=60)
    return float(w.sum()),float((w*hit[j,k]*ok).sum())

def hits(rows,T):
    return np.array([S is not None and abs(d-s)<=1 and S>T for dt,R,s,S,d in rows],float).reshape(len(DT),len(RF))

if __name__=='__main__':
    catfile=sys.argv[1]
    recs=[json.loads(l) for f in glob.glob('shards/**/*.jsonl',recursive=True) for l in open(f)]
    recs={r['name']:r for r in recs}
    cat={r['tns_name']:r for r in csv.DictReader(open(catfile)) if r['sub_num']=='0'}
    def num(n,k):
        try: return float(cat[n][k])
        except Exception: return float('nan')
    manifest=open('manifest.txt').read().split()
    out=dict(files_listed=len(manifest),records=len(recs),missing=sorted(set(manifest)-set(recs)),
             errors={n:r.get('err') or r.get('inj_err') for n,r in recs.items() if 'err' in r or 'inj_err' in r})
    old={o['name']:o for o in json.load(open('old_scores.json')) if 'err' not in o}
    cmpn=[n for n in old if n in recs and 'err' not in recs[n]]
    def agree(n):  # S None here and -inf in the old file both mean 'no valid echo'
        new,o=recs[n].get('S'),old[n]['S']
        return (new is None and not np.isfinite(o)) or (new is not None and np.isfinite(o) and abs(new-o)<0.01)
    same=sum(agree(n) for n in cmpn)
    out['reproduction']=dict(old_scored=len(old),compared=len(cmpn),matching=int(same),passed=bool(same>=0.99*len(old)),
        mismatches=[n for n in cmpn if not agree(n)][:30])
    cal=[r for r in recs.values() if r['rep'] and 'err' not in r]
    search_all=[n for n in manifest if n in cat and cat[n]['repeater_name'] in ('','-9999') and num(n,'snr_fitb')>0]
    sr=[recs[n] for n in search_all if n in recs and 'err' not in recs[n] and 'inj' in recs[n]]
    out['n_calibration']=len(cal); out['n_search_listed']=len(search_all); out['n_search_scored']=len(sr)
    out['search_fail_fraction']=1-len(sr)/max(len(search_all),1)
    out['limit_complete']=bool(out['search_fail_fraction']<=0.05)
    top=sorted([r for r in cal if r.get('S') is not None],key=lambda r:-r['S'])[:15]
    out['top_calibration']=[dict(name=r['name'],rep=r['rep'],S=r['S'],d=r['d']) for r in top]
    Tp=max(23.0,top[0]['S']) if top else float('nan')
    masses=[100,300,1000,3000]
    pts={M:[weight_points(M,num(r['name'],'dm_exc_ymw16'),num(r['name'],'width_fitb')*1e3,num(r['name'],'snr_fitb')) for r in sr] for M in masses}
    clean_ok=[not any(num(r['name'],k)==1 for k in ('excluded_flag',)) for r in sr]
    cal_listed=[n for n in manifest if n in cat and cat[n]['repeater_name'] not in ('','-9999')]
    out['n_calibration_listed']=len(cal_listed); out['calibration_fail_fraction']=1-len(cal)/max(len(cal_listed),1)
    out['limit_complete']=bool(out['limit_complete'] and out['calibration_fail_fraction']<=0.01 and out['reproduction']['passed'])
    variants=[('primary',Tp,'inj',False),('secondary_T23',23.0,'inj',False),('sensitivity_drop_excluded_flag',Tp,'inj',True),('comparison_old_noisy_injection',Tp,'inj_noisy',False)]
    for label,T,key,drop in variants:
        idx=[i for i in range(len(sr)) if (clean_ok[i] or not drop)]
        H={i:hits(sr[i][key],T) for i in idx}
        det=[]
        for i in idx:
            snr=num(sr[i]['name'],'snr_fitb')
            det+= [H[i][a,b] for a in range(len(DT)) for b in range(len(RF)) if RF[b]<=5 and snr/RF[b]>=10]
        cands=sorted([dict(name=r['name'],S=r['S'],d=r['d'],a=r['a'],chi=r['chi'],snr=num(r['name'],'snr_fitb')) for i,r in enumerate(sr) if i in H and r.get('S') is not None and r['S']>T],key=lambda x:-x['S'])
        k=len(cands); N95=0.5*chi2.ppf(0.95,2*k+2); lim={}
        for M in masses:
            tt=[burst_tau(pts[M][i],H[i]) for i in idx]
            tp=sum(x[0] for x in tt); te=sum(x[1] for x in tt)
            lim[M]=dict(SumTau_if_perfect=tp,SumTau_eff=te,f95=float(N95/te) if te>0 else None)
        out[label]=dict(T=float(T),k=k,N95=float(N95),n_bursts=len(idx),candidates=cands,efficiency_snr2_ge10_R_le5=float(np.mean(det)) if det else None,n_injections=len(det),limits=lim)
    s=[r for r in sr if r.get('S') is not None]
    out['top_search']=[dict(name=r['name'],S=r['S'],d=r['d']) for r in sorted(s,key=lambda r:-r['S'])[:15]]
    os_=__import__('os'); os_.makedirs('results',exist_ok=True)
    json.dump(out,open('results/summary.json','w'),indent=1)
    import gzip
    with gzip.open('results/scores.json.gz','wt') as f: json.dump(list(recs.values()),f)
    print(json.dumps({k:out[k] for k in ['records','n_calibration','n_search_scored','search_fail_fraction','reproduction','primary']},indent=1)[:5000])
