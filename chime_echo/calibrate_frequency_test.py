import numpy as np,h5py,sys,warnings
from scipy.stats import chi2
sys.path.insert(0,'/home/claude/darksieve/chime_echo'); warnings.filterwarnings('ignore')
from echo import load_sub,whiten,main_window
with h5py.File('/home/claude/dl/bad/FRB20190131D.h5') as f:
    D=np.array(f['data'],float); M=np.nan_to_num(np.array(f['model'],float)); G=np.array(f['good_freq'])
    if 'calibration/good_freq' in f: G=G&np.array(f['calibration/good_freq'])
    FL=np.array(f['flag'])
D[~FL]=np.nan; D[~G]=np.nan; nt=D.shape[1]
X0=load_sub(None,D=D); pk=int(np.argmax(np.nansum(X0,0))); on=np.zeros(nt,bool); on[max(pk-8,0):min(pk+75,nt)]=True
Y0=whiten(X0,on); lo,hi,_=main_window(Y0); mu=X0[:,~on].mean(1,keepdims=True); sdv=X0[:,~on].std(1,keepdims=True)
def test(Y,start,ng):
    A=Y[:,lo:hi+1];B=Y[:,start:start+hi-lo+1]; a=[];e=[]
    for g in np.array_split(np.arange(Y.shape[0]),ng):
        sig=max(np.sum(A[g]**2)-A[g].size,1e-9); ai=np.sum(A[g]*B[g])/sig
        a.append(ai); e.append(np.sqrt(sig*(1+ai**2)+A[g].size)/sig)
    a=np.array(a);e=np.array(e);use=e<1;W=1/e[use]**2;m=(W*a[use]).sum()/W.sum()
    return chi2.sf((((a[use]-m)/e[use])**2).sum(),use.sum()-1)
main=np.zeros_like(M); main[:,lo:hi+1]=M[:,lo:hi+1]
starts=[s for s in range(0,nt-(hi-lo+1)) if (s+(hi-lo+1)<lo-4) or s>hi+20]
real={ng:test(Y0,lo+9,ng) for ng in (4,8,16)}
sim={ng:[] for ng in (4,8,16)}
for s in starts:
    ec=np.zeros_like(M); ec[:,s:s+hi-lo+1]=0.41*main[:,lo:hi+1]
    X=load_sub(None,D=D+ec); Y=(X-mu)/sdv
    for ng in sim: sim[ng].append(test(Y,s,ng))
for ng in sim:
    p=np.array(sim[ng]); k=(p<=real[ng]).sum()
    print('%2d bands: real p=%.4f ; %d of %d planted true copies look this uneven or worse (%.1f%%)'%(ng,real[ng],k,len(p),100*(k+1)/(len(p)+1)))
