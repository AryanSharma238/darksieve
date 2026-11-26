import sys,csv,json,gzip,numpy as np,h5py,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, warnings
from scipy.stats import chi2
from scipy.ndimage import median_filter
sys.path.insert(0,'/home/claude/darksieve/chime_echo'); warnings.filterwarnings('ignore')
from echo import load_sub,whiten,main_window
from echo2 import scan_all
rows=list(csv.DictReader(open('/home/claude/dl/chimefrbcat2.csv')))
recs=json.load(gzip.open('/home/claude/darksieve/chime_echo/out_v2/scores2.json.gz','rt'))
sc=lambda r: r['real']['hp9']['S'] if r['real']['hp9']['S'] is not None else -np.inf
snr={r['tns_name']:float(r['snr_fitb'] or 0) for r in rows if r['sub_num']=='0'}
# how often do repeat bursts and one-off bursts of similar brightness reach high scores?
for lo in (0,30,50):
    for lab,sel in (('repeat',[r for r in recs if r['rep'] and 'real' in r]),('one-off',[r for r in recs if not r['rep'] and 'real' in r and snr[r['name']]>0])):
        x=[sc(r) for r in sel if snr[r['name']]>=lo]
        print(f'snr>={lo:2d} {lab:8s} n={len(x):4d} frac>5: {np.mean(np.array(x)>5):.4f}  >8: {sum(np.array(x)>8)}  max {max(x):.1f}')
names=['FRB20190131D','FRB20210113C','FRB20190805E','FRB20200302C']
fig,ax=plt.subplots(4,3,figsize=(15,13))
for i,n in enumerate(names):
    R=[r for r in rows if r['tns_name']==n]
    with h5py.File(f'/home/claude/dl/bad/{n}.h5') as f:
        D=np.array(f['data'],float); G=np.array(f['good_freq'])
        if 'calibration/good_freq' in f: G=G&np.array(f['calibration/good_freq'])
        FL=np.array(f['flag'])
    D[~FL]=np.nan; D[~G]=np.nan
    X=load_sub(None,D=D); best,w,_=scan_all(X); b=best['hp9']; d=b['d']; nt=X.shape[1]
    pk=int(np.argmax(np.nansum(X,0))); on=np.zeros(nt,bool); on[max(pk-8,0):min(pk+75,nt)]=True
    Y=whiten(X,on); lo,hi,_=main_window(Y)
    Z=Y-median_filter(Y,size=(1,9),mode='nearest'); sd=Z[:,~on].std(1,keepdims=True); sd[sd==0]=1; Z=Z/sd
    print('\n',n,'hp9 S=%.2f delay=%.2f ms R=%.2f chi=%.3f window=%d bins; versions:'%(b['S'],d*.983,1/b['a'],b['chi'],w),{v:(None if x['S'] is None else round(x['S'],1),x['d']) for v,x in best.items()})
    for r in R: print('   catalog sub-burst',r['sub_num'],'t_ms=%.2f'%((float(r['mjd_inf'])-float(R[0]['mjd_inf']))*86400e3),'width_ms=%.2f'%(float(r['width_fitb'])*1e3),'sp_idx=%.1f+-%.1f'%(float(r['sp_idx']),float(r['sp_idx_err'] or 0)),'run=%.1f+-%.1f'%(float(r['sp_run']),float(r['sp_run_err'] or 0)))
    for lab,M in (('raw',Y),('trend removed',Z)):
        A=M[:,lo:hi+1];B=M[:,lo+d:hi+1+d]
        for ng in (4,8,16):
            a=[];e=[]
            for g in np.array_split(np.arange(M.shape[0]),ng):
                sig=max(np.sum(A[g]**2)-A[g].size,1e-9); ai=np.sum(A[g]*B[g])/sig
                a.append(ai); e.append(np.sqrt(sig*(1+ai**2)+A[g].size)/sig)
            a=np.array(a);e=np.array(e);use=e<1; W=1/e[use]**2; m=(W*a[use]).sum()/W.sum(); c=(((a[use]-m)/e[use])**2).sum()
            p=chi2.sf(c,use.sum()-1) if use.sum()>1 else np.nan
            print('   %-13s %2d bands: ratio %.2f, same in every band? p=%.2g (%d usable)'%(lab,ng,m,p,use.sum()))
            if ng==8 and lab=='raw': ax[i,2].errorbar(np.where(use)[0],a[use],e[use],fmt='o'); ax[i,2].axhline(m,c='k',ls='--'); ax[i,2].set_title('second/first ratio by frequency (flat if lensed)')
    t=(np.arange(nt)-pk)*.983; s0=max(pk-20,0); s1=min(pk+50,nt)
    ax[i,0].plot(t[s0:s1],(Y.sum(0)/np.sqrt(len(Y)))[s0:s1]); ax[i,0].axvline(d*.983,c='r',alpha=.4); ax[i,0].set_title('%s  S=%.1f, delay %.1f ms'%(n,b['S'],d*.983))
    ax[i,1].imshow(Y[:,s0:s1],aspect='auto',origin='lower',extent=[t[s0],t[s1-1],0,len(Y)],vmin=-2,vmax=6); ax[i,1].set_title('time vs frequency')
plt.tight_layout(); plt.savefig('/home/claude/work/candidates_v2.png',dpi=105)
