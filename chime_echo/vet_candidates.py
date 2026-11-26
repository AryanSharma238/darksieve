import sys,csv,json,numpy as np,h5py,matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt, warnings
sys.path.insert(0,'/home/claude/darksieve/chime_echo'); warnings.filterwarnings('ignore')
from echo import load_sub,echo_scan,whiten,main_window
rows=list(csv.DictReader(open('/home/claude/dl/chimefrbcat2.csv')))
names=['FRB20200306B','FRB20220219A','FRB20220720B','FRB20211006D']
fig,ax=plt.subplots(4,3,figsize=(15,13))
for i,n in enumerate(names):
    R=[r for r in rows if r['tns_name']==n]
    print('\n',n,'catalog sub-bursts:',len(R),'excluded',R[0]['excluded_flag'],'sidelobe',R[0]['sidelobe_flag'],'intrachan',R[0]['intrachan_flag'],'cat1',R[0]['catalog1_flag'],'DMexc',R[0]['dm_exc_ymw16'][:6])
    for r in R: print('   sub',r['sub_num'],'mjd_ms',round((float(r['mjd_inf'])-float(R[0]['mjd_inf']))*86400e3,2),'width_ms',round(float(r['width_fitb'])*1e3,2),'sp_idx',round(float(r['sp_idx']),1),'run',round(float(r['sp_run']),1),'scat',r['scat_time'][:6])
    with h5py.File(f'/home/claude/dl/bad/{n}.h5') as f:
        D=np.array(f['data'],float); M=np.array(f['model'],float); G=np.array(f['good_freq']); G=G&np.array(f['calibration/good_freq']) if 'calibration/good_freq' in f else G; FL=np.array(f['flag'])
    D[~FL]=np.nan; D[~G]=np.nan
    X=load_sub(None,D=D); b=echo_scan(X,cmax=1.35); d=b['d']; nt=X.shape[1]
    pk=int(np.argmax(np.nansum(X,0))); on=np.zeros(nt,bool); on[max(pk-8,0):min(pk+75,nt)]=True
    Y=whiten(X,on); lo,hi,_=main_window(Y); A=Y[:,lo:hi+1]; B=Y[:,lo+d:hi+1+d]
    print('   S=%.1f delay=%.2f ms window=%d..%d a=%.2f R=%.2f chi=%.3f'%(b['S'],d*0.983,lo,hi,b['a'],1/b['a'],b['chi']))
    # spectral test: does every frequency group need the same echo amplitude?
    for ng in (4,8,16):
        gs=np.array_split(np.arange(Y.shape[0]),ng); a=[];e=[]
        for g in gs:
            aa=np.sum(A[g]**2); ab=np.sum(A[g]*B[g]); ai=ab/aa
            # errors from noise in both windows; debias template power by its noise
            sig=max(aa-A[g].size,1e-9); ai=ab/sig
            e.append(np.sqrt(sig*(1+ai**2)+A[g].size)/sig); a.append(ai)
        a=np.array(a);e=np.array(e); use=e<1; w=1/e[use]**2; m=(w*a[use]).sum()/w.sum(); c=(((a[use]-m)/e[use])**2).sum()
        from scipy.stats import chi2
        print('   %2d bands: mean a=%.3f chi2=%.1f dof=%d p=%.2g'%(ng,m,c,use.sum()-1,chi2.sf(c,use.sum()-1)))
        if ng==8: ax[i,2].errorbar(np.where(use)[0],a[use],e[use],fmt='o'); ax[i,2].axhline(m,c='k',ls='--'); ax[i,2].set_title('echo/main ratio per band (flat if lensed)')
    t=(np.arange(nt)-pk)*0.983; ts=Y.sum(0)/np.sqrt(Y.shape[0])
    ax[i,0].plot(t,ts,label='data'); ax[i,0].set_xlim(-25,70); ax[i,0].axvline(d*0.983,c='r',alpha=.4); ax[i,0].set_title(n+'  S=%.1f, delay %.1f ms'%(b['S'],d*.983))
    ax[i,1].imshow(Y[:,max(pk-25,0):pk+70],aspect='auto',origin='lower',extent=[t[max(pk-25,0)],t[min(pk+69,nt-1)],0,Y.shape[0]],vmin=-2,vmax=6); ax[i,1].set_title('time vs frequency sub-band')
plt.tight_layout(); plt.savefig('/home/claude/work/candidates.png',dpi=110)
