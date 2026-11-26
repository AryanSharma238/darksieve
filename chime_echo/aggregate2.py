import csv,json,glob,sys,hashlib,gzip,os,numpy as np
from scipy.stats import chi2
from aggregate import weight_points,burst_tau,hits
from echo2 import VARIANTS
MASSES=[100,300,1000,3000]
def half(key): return 'dev' if hashlib.md5(key.encode()).digest()[0]%2==0 else 'test'
if __name__=='__main__':
    cat={r['tns_name']:r for r in csv.DictReader(open(sys.argv[1])) if r['sub_num']=='0'}
    num=lambda n,k: float(cat[n][k] or 'nan')
    recs={}
    for f in glob.glob('shards/**/*.jsonl',recursive=True):
        for l in open(f):
            r=json.loads(l); recs[r['name']]=r
    cal=[r for r in recs.values() if r['rep'] and 'real' in r]
    sr=[r for r in recs.values() if not r['rep'] and 'inj' in r and num(r['name'],'snr_fitb')>0]
    cal_listed=sum(1 for r in recs.values() if r['rep']); s_listed=sum(1 for n,r in recs.items() if not r['rep'] and num(n,'snr_fitb')>0)
    pts={M:{r['name']:weight_points(M,num(r['name'],'dm_exc_ymw16'),num(r['name'],'width_fitb')*1e3,num(r['name'],'snr_fitb')) for r in sr} for M in MASSES}
    S=lambda r,v: -np.inf if r['real'][v]['S'] is None else r['real'][v]['S']
    def evaluate(v,T,rows):
        H={r['name']:hits(r['inj'][v],T) for r in rows}
        tau={M:float(sum(burst_tau(pts[M][r['name']],H[r['name']])[1] for r in rows)) for M in MASSES}
        det=[H[r['name']][a,b] for r in rows for a in range(11) for b in range(5) if num(r['name'],'snr_fitb')/[1.2,1.5,2,3,5][b]>=10]
        return tau,float(np.mean(det))
    out=dict(n_cal=len(cal),cal_listed=cal_listed,n_search=len(sr),search_listed=s_listed)
    # 1. choose the version using ONLY the development half
    dev={}
    for v in VARIANTS:
        Td=max(S(r,v) for r in cal if half(r['rep'])=='dev')
        tau,eff=evaluate(v,Td,[r for r in sr if half(r['name'])=='dev'])
        dev[v]=dict(T=float(Td),efficiency=eff,tau=tau,metric=tau[300]+tau[1000])
    win=max(VARIANTS,key=lambda v:dev[v]['metric']); out['dev']=dev; out['winner']=win
    # 2. held-out half: winner vs v1, thresholds from the development repeaters
    test=[r for r in sr if half(r['name'])=='test']
    out['heldout']={}
    for v in dict.fromkeys(['v1',win]):
        tau,eff=evaluate(v,dev[v]['T'],test)
        fa=[r['name'] for r in cal if half(r['rep'])=='test' and S(r,v)>dev[v]['T']]
        out['heldout'][v]=dict(efficiency=eff,tau=tau,heldout_repeats_above_T=len(fa),n_heldout_repeats=sum(half(r['rep'])=='test' for r in cal),search_above_T=sum(S(r,v)>dev[v]['T'] for r in test))
    # 3. full search, every version reported; only the winner is the headline
    out['full']={}
    for v in VARIANTS:
        top=sorted(cal,key=lambda r:-S(r,v))[:5]; T=S(top[0],v)
        tau,eff=evaluate(v,T,sr)
        c=sorted([dict(name=r['name'],S=S(r,v),**{k:r['real'][v][k] for k in 'dac'.replace('c','')},snr=num(r['name'],'snr_fitb'),scat_ms=num(r['name'],'scat_time')*1e3) for r in sr if S(r,v)>T],key=lambda x:-x['S'])
        N95=0.5*chi2.ppf(0.95,2*len(c)+2)
        out['full'][v]=dict(T=float(T),efficiency=eff,candidates=c,f95={M:(float(N95/tau[M]) if tau[M]>0 else None) for M in MASSES},tau=tau,
            top_calibration=[(r['name'],r['rep'],S(r,v)) for r in top],
            top_search=[(r['name'],S(r,v),r['real'][v]['d']) for r in sorted(sr,key=lambda r:-S(r,v))[:8]])
    os.makedirs('results',exist_ok=True); json.dump(out,open('results/summary2.json','w'),indent=1)
    with gzip.open('results/scores2.json.gz','wt') as f: json.dump(list(recs.values()),f)
    print(json.dumps(dict(winner=win,heldout=out['heldout'],full=out['full'][win]),indent=1)[:4000])
