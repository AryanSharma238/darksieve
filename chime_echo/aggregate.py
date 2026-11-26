import csv,json,glob,sys,numpy as np
from scipy.stats import chi2
from tau import zm,dt_ms,ymin,ymax_R,Dc,E,Oc,H0,c
sb=[0,25,40,70,120,1e9]; db=[0,3,5,10,20,40,70]; rb=[1,1.35,1.75,2.5,4,1e9]
def weights(M,src,nz=20,ny=30):
    # same integral as the earlier boot.py
    W=np.zeros((5,6,5))
    for d,wi,s in src:
        zs=zm(d)
        if not np.isfinite(zs) or not s/10>=1: continue
        i=min(max(np.searchsorted(sb,s)-1,0),4); ym=ymax_R(s/10)
        for z in np.linspace(zs/nz/2,zs*(1-1/nz/2),nz):
            Dl=Dc(z)/(1+z); Ds=Dc(zs)/(1+zs); Dls=(Dc(zs)-Dc(z))/(1+zs); y0=ymin(z,M,wi)
            if ym<=y0: continue
            pref=1.5*Oc*H0**2/(c*H0*E(z))*Dls*Dl/Ds*(1+z)**2*(zs/nz)
            ys=np.linspace(y0,ym,ny); dy=ys[1]-ys[0]
            for y in ys:
                R=(y*y+2+y*np.sqrt(y*y+4))/(y*y+2-y*np.sqrt(y*y+4)); t=dt_ms(y,M,z)
                j=np.searchsorted(db,t)-1; k=min(np.searchsorted(rb,R)-1,4)
                if 0<=j<=5 and k>=0: W[i,j,k]+=pref*2*y*dy
    return W

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
    cmpn=[n for n in old if n in recs and recs[n].get('S') is not None]
    same=sum(abs(recs[n]['S']-old[n]['S'])<0.01 for n in cmpn)
    out['reproduction']=dict(old_scored=len(old),compared=len(cmpn),matching=int(same),passed=bool(same>=0.99*len(old)),
        mismatches=[n for n in cmpn if abs(recs[n]['S']-old[n]['S'])>=0.01][:30])
    cal=[r for r in recs.values() if r['rep'] and 'err' not in r and r.get('S') is not None]
    search_all=[n for n in manifest if n in cat and cat[n]['repeater_name'] in ('','-9999') and num(n,'snr_fitb')>0]
    sr=[recs[n] for n in search_all if n in recs and 'err' not in recs[n] and 'inj' in recs[n]]
    out['n_calibration']=len(cal); out['n_search_listed']=len(search_all); out['n_search_scored']=len(sr)
    out['search_fail_fraction']=1-len(sr)/max(len(search_all),1)
    out['limit_complete']=bool(out['search_fail_fraction']<=0.05)
    top=sorted(cal,key=lambda r:-r['S'])[:15]
    out['top_calibration']=[dict(name=r['name'],rep=r['rep'],S=r['S'],d=r['d']) for r in top]
    Tp=max(23.0,top[0]['S']) if top else float('nan')
    src=[(num(r['name'],'dm_exc_ymw16'),num(r['name'],'width_fitb')*1e3,num(r['name'],'snr_fitb')) for r in sr]
    Ws={M:weights(M,src) for M in [100,300,1000,3000]}
    for label,T in [('primary',Tp),('secondary_T23',23.0)]:
        cnt=np.zeros((5,6,5)); rec=np.zeros((5,6,5)); det_all=[]
        for r in sr:
            snr=num(r['name'],'snr_fitb'); i=np.searchsorted(sb,snr)-1
            for dt,R,s,S,d in r['inj']:
                j=np.searchsorted(db,dt)-1; k=np.searchsorted(rb,R)-1
                hit=S is not None and abs(d-s)<=1 and S>T
                cnt[i,j,k]+=1; rec[i,j,k]+=hit
                if snr/R>=10: det_all.append(hit)
        Ef=np.where(cnt>0,rec/np.maximum(cnt,1),0)
        cands=sorted([dict(name=r['name'],S=r['S'],d=r['d'],a=r['a'],chi=r['chi'],snr=r['snr']) for r in sr if r.get('S') is not None and r['S']>T],key=lambda x:-x['S'])
        k=len(cands); N95=0.5*chi2.ppf(0.95,2*k+2); lim={}
        for M,W in Ws.items():
            te=float((W*Ef).sum())
            lim[M]=dict(SumTau_if_perfect=float(W.sum()),SumTau_eff=te,f95=float(N95/te) if te>0 else None)
        out[label]=dict(T=float(T),k=k,N95=float(N95),candidates=cands,efficiency=float(np.mean(det_all)) if det_all else None,n_injections=len(det_all),limits=lim)
    s=[r for r in sr if r.get('S') is not None]
    out['top_search']=[dict(name=r['name'],S=r['S'],d=r['d']) for r in sorted(s,key=lambda r:-r['S'])[:15]]
    os_=__import__('os'); os_.makedirs('results',exist_ok=True)
    json.dump(out,open('results/summary.json','w'),indent=1)
    import gzip
    with gzip.open('results/scores.json.gz','wt') as f: json.dump([{k:v for k,v in r.items() if k!='inj'} for r in recs.values()],f)
    print(json.dumps({k:out[k] for k in ['records','n_calibration','n_search_scored','search_fail_fraction','reproduction','primary']},indent=1)[:5000])
