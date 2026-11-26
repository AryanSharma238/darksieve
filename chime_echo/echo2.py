"""Variants of the echo detector. 'v1' equals echo.py. 'dipN' only accepts a delay if the summed light curve
falls after the main peak and then rises again by more than N noise units at the echo. 'hpW' removes anything
smoother than a running median of W samples from every sub-band before matching."""
import numpy as np
from scipy.ndimage import median_filter, uniform_filter1d
from echo import whiten, main_window
VARIANTS=['v1','dip1','dip2','dip3','hp9','hp21']

def scan_all(X,dmin=2,dmax=70,cmax=1.35):
    nt=X.shape[1]
    ts0=np.nansum(X,0); pk0=int(np.argmax(ts0)); on=np.zeros(nt,bool); on[max(pk0-8,0):min(pk0+75,nt)]=True
    Y=whiten(X,on); lo,hi,snr=main_window(Y); w=hi-lo+1
    H={}
    for W in (9,21):
        h=Y-median_filter(Y,size=(1,W),mode='nearest'); sd=h[:,~on].std(1,keepdims=True); sd[sd==0]=1; H['hp%d'%W]=h/sd
    ts=Y.sum(0); ts=ts/(ts[~on].std()+1e-9); sm=uniform_filter1d(ts,3); sm/=sm[~on].std()+1e-9
    pk=lo+int(np.argmax(sm[lo:hi+1]))
    best={v:dict(S=-np.inf,d=0,a=0.,chi=np.inf) for v in VARIANTS}
    def match(Z,d):
        A=Z[:,lo:hi+1]; B=Z[:,lo+d:hi+1+d]; aa=np.sum(A*A); a=np.sum(A*B)/aa
        r=B-a*A; return a*np.sqrt(aa), a, np.sum(r*r)/(1+a*a)/(r.size-1)
    def keep(v,S,a,chi,d):
        if 0<a<1 and chi<cmax and S>best[v]['S']: best[v]=dict(S=float(S),d=int(d),a=float(a),chi=float(chi))
    for d in range(max(dmin,w),dmax+1):
        if hi+d>=nt: break
        S,a,chi=match(Y,d); keep('v1',S,a,chi,d)
        e0=max(pk+d-2,pk+1); seg=sm[pk+1:pk+d]
        if len(seg):
            j=pk+1+int(np.argmin(seg)); after=sm[max(e0,j+1):pk+d+3]
            rise=after.max()-sm[j] if len(after) else -np.inf
            for n in (1,2,3):
                if rise>n: keep('dip%d'%n,S,a,chi,d)
        for v,Z in H.items(): keep(v,*match(Z,d),d)
    for b in best.values():
        if not np.isfinite(b['S']): b['S']=None; b['chi']=None
    return best,w,snr
