import os,sys,csv,json,time,hashlib,tempfile,urllib.request,warnings
from concurrent.futures import ThreadPoolExecutor
import numpy as np, h5py
from echo import load_sub,echo_scan
warnings.filterwarnings('ignore')
BASE='https://cadc-west-01.canfar.net/vault/files/AstroDataCitationDOI/CISTI.CANFAR/25.0066/data/'
DT=[2,3,4,6,8,12,16,24,32,48,64]; RF=[1.2,1.5,2.0,3.0,5.0,7.0,10.0,15.0,25.0,50.0]

def fetch(url,path,tries=5):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url,timeout=120) as r:
                data=r.read(); n=r.headers.get('Content-Length')
            if n and int(n)!=len(data): raise IOError('short read')
            open(path,'wb').write(data); return
        except Exception:
            if i==tries-1: raise
            time.sleep(3*2**i)

def clean(b): return {k:(None if isinstance(v,(float,np.floating)) and not np.isfinite(v) else (float(v) if isinstance(v,(float,np.floating)) else int(v))) for k,v in b.items()}

def one(name,meta,tmp):
    rec=dict(name=name,rep=meta.get('rep',''))
    p=os.path.join(tmp,name+'.h5')
    try:
        fetch(f'{BASE}dynamic_spectra/data/{name}_stokesi_dynamic_spectrum.h5',p)
        with h5py.File(p,'r') as f:
            D=np.array(f['data'],float); G=np.array(f['good_freq'])&np.array(f['calibration/good_freq']); FL=np.array(f['flag'])
            M=np.array(f['model'],float) if not rec['rep'] else None
        D[~FL]=np.nan; D[~G]=np.nan
        rec.update(clean(echo_scan(load_sub(None,D=D),cmax=1.35)))
    except Exception as e:
        rec['err']=repr(e)[:150]
    if 'err' not in rec and not rec['rep']:
        try:
            nt=D.shape[1]; Mz=np.nan_to_num(M); mts=np.abs(Mz).sum(0)
            on=np.convolve(mts>0.02*mts.max(),np.ones(5),'same')>0
            # 'inj' (headline): noise-free catalog model as echo. 'inj_noisy' (comparison only): old method, copies host noise too.
            inj={}
            for key,src in (('inj',Mz),('inj_noisy',np.nan_to_num(D))):
                sig=np.where(on[None,:],src,0.); rows=[]
                for dt in DT:
                    s=int(round(dt/0.983))
                    for R in RF:
                        ec=np.zeros_like(sig)
                        if s<nt: ec[:,s:]=sig[:,:nt-s]/R
                        b=echo_scan(load_sub(None,D=D+ec),cmax=1.35)
                        S=float(b['S']); rows.append([dt,R,s,S if np.isfinite(S) else None,int(b['d'])])
                inj[key]=rows
            rec.update(inj)
        except Exception as e:
            rec['inj_err']=repr(e)[:150]
    if os.path.exists(p): os.remove(p)
    return rec

if __name__=='__main__':
    shard,nshard,mode=int(sys.argv[1]),int(sys.argv[2]),sys.argv[3]
    tmp=tempfile.mkdtemp()
    cpath=os.path.join(tmp,'cat.csv'); fetch(BASE+'table/chimefrbcat2.csv',cpath)
    meta={}
    for r in csv.DictReader(open(cpath)):
        if r['sub_num']=='0':
            rn=r['repeater_name']; meta[r['tns_name']]=dict(rep='' if rn in ('','-9999') else rn)
    names=open('manifest.txt').read().split()
    if mode=='pilot': names=['FRB20180814B','FRB20180917A','FRB20190131D','FRB20211115A']+[n for n in names[:20]]; names=list(dict.fromkeys(names))
    mine=names[shard::nshard]; t=time.time()
    with ThreadPoolExecutor(3) as ex, open(f'shard_{shard}.jsonl','w') as out:
        # downloads overlap in threads; h5py serialises file reads itself
        for i,rec in enumerate(ex.map(lambda n: one(n,meta.get(n,{'rep':''}),tmp) if n in meta else dict(name=n,rep='',err='not in catalog'),mine)):
            out.write(json.dumps(rec)+'\n'); out.flush()
            if i%25==0: print(f'{i+1}/{len(mine)} {time.time()-t:.0f}s',flush=True)
