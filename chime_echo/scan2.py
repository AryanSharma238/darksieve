import os,sys,csv,json,time,tempfile,warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, h5py
from echo import load_sub
from echo2 import scan_all,VARIANTS
from scan import fetch,BASE,DT,RF
warnings.filterwarnings('ignore')

def one(name,rep,tmp):
    rec=dict(name=name,rep=rep); p=os.path.join(tmp,name+'.h5')
    try:
        fetch(f'{BASE}dynamic_spectra/data/{name}_stokesi_dynamic_spectrum.h5',p)
        with h5py.File(p,'r') as f:
            D=np.array(f['data'],float); G=np.array(f['good_freq'])
            if 'calibration/good_freq' in f: G=G&np.array(f['calibration/good_freq'])
            FL=np.array(f['flag']); M=np.array(f['model'],float) if not rep else None
        D[~FL]=np.nan; D[~G]=np.nan
        best,w,snr=scan_all(load_sub(None,D=D)); rec.update(real=best,w=int(w))
        if not rep:
            nt=D.shape[1]; Mz=np.nan_to_num(M); mts=np.abs(Mz).sum(0)
            on=np.convolve(mts>0.02*mts.max(),np.ones(5),'same')>0; sig=np.where(on[None,:],Mz,0.)
            inj={v:[] for v in VARIANTS}
            for dt in DT:
                s=int(round(dt/0.983))
                for R in RF:
                    ec=np.zeros_like(sig)
                    if s<nt: ec[:,s:]=sig[:,:nt-s]/R
                    b=scan_all(load_sub(None,D=D+ec))[0]
                    for v in VARIANTS: inj[v].append([dt,R,s,None if b[v]['S'] is None else round(b[v]['S'],3),b[v]['d']])
            rec['inj']=inj
    except Exception as e:
        rec['err']=repr(e)[:150]
    if os.path.exists(p): os.remove(p)
    return rec

if __name__=='__main__':
    shard,nshard=int(sys.argv[1]),int(sys.argv[2]); tmp=tempfile.mkdtemp()
    cpath=os.path.join(tmp,'cat.csv'); fetch(BASE+'table/chimefrbcat2.csv',cpath)
    meta={r['tns_name']:('' if r['repeater_name'] in ('','-9999') else r['repeater_name']) for r in csv.DictReader(open(cpath)) if r['sub_num']=='0'}
    mine=[n for n in open('manifest.txt').read().split() if n in meta][shard::nshard]; t=time.time()
    with ThreadPoolExecutor(3) as ex, open(f'shard_{shard}.jsonl','w') as out:
        for i,rec in enumerate(ex.map(lambda n: one(n,meta[n],tmp),mine)):
            out.write(json.dumps(rec)+'\n'); out.flush()
            if i%25==0: print(f'{i+1}/{len(mine)} {time.time()-t:.0f}s',flush=True)
