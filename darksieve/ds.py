"""DarkSieve core: detection + validation of subhalo candidates in one strong lens.
Stages
  D1 macro refit (PL+shear+m3,m4, adaptive pixel source)
  D2 fast fixed-macro screen over arc positions (clump masses 8.5..10)
  D3 verification of top screen peaks: joint (macro+clump) vs equal-budget smooth refit, mirror control, position check
Validation of each verified candidate
  V1 Hide-and-Predict: hold out pixels near the candidate, refit with/without the clump, compare chi2 on held-out pixels
  V2 split test: two independent halves of the data (independent noise) must both prefer the clump at the same place
  V3 lopsided adversary: smooth model gets an extra m=1 multipole; clump must still win
Every statistic is returned so thresholds can be calibrated on fakes and frozen before testing."""
import numpy as np, copy, json
import autolens as al, pj
from model import Lens, logprior, MACRO, MP
import stages as S
SIG=0.01; MASSES=(8.5,9.0,9.5,10.0)
BASE=[]   # clumps always present (e.g. already-known subhaloes), fixed

# ---------- generic evaluation with arbitrary clumps / extra profiles ----------
def tracer(L,Q,clumps=(),extra=None):
    g=list(L.tracer(Q,None,extra=extra).galaxies)
    for c in list(BASE)+list(clumps): g.insert(1,pj.pj_galaxy(c[0]+0.0013,c[1]+0.0017,10**c[2],L.zl,L.zs,Q['te']))
    return al.Tracer(galaxies=g)
def logev(L,Q,clumps=(),extra=None,ds=None):
    try:
        v=al.FitImaging(dataset=ds or L.ds,tracer=tracer(L,Q,clumps,extra)).log_evidence
        return v if np.isfinite(v) else -1e9
    except Exception: return -1e9
def post(L,Q,clumps=(),extra=None,ds=None): return logev(L,Q,clumps,extra,ds)+logprior(Q,SIG)
def m1_extra(Q):
    if 'a1c' not in Q: return None
    return dict(mp1=al.mp.PowerLawMultipole(m=1,centre=(Q['my'],Q['mx']),einstein_radius=Q['te'],slope=2.05,multipole_comps=(Q['a1c'],Q['a1s'])))
def nm_fit(f,Q,keys,clump,nfev):
    """Nelder-Mead over macro keys (+ clump y,x,lM)."""
    from scipy.optimize import minimize
    BND=dict(S.BND,a1c=(-.06,.06),a1s=(-.06,.06)); STEP=dict(S.STEP,a1c=.005,a1s=.005)
    ks=list(keys)+(['cy','cx','lM'] if clump is not None else [])
    x0=np.array([Q[k] for k in keys]+(list(clump) if clump is not None else []),float)
    sim=np.array([x0]+[x0+np.eye(len(x0))[i]*STEP[ks[i]] for i in range(len(x0))])
    def g(x):
        if any(not BND[k][0]<=v<=BND[k][1] for k,v in zip(ks,x)): return 1e9
        d=dict(zip(ks,x)); QQ=dict(Q,**{k:d[k] for k in keys}); c=(d['cy'],d['cx'],d['lM']) if clump is not None else None
        return -f(QQ,c)
    r=minimize(g,x0,method='Nelder-Mead',options={'initial_simplex':sim,'maxfev':nfev})
    d=dict(zip(ks,map(float,r.x))); QQ=dict(Q,**{k:d[k] for k in keys}); c=(d['cy'],d['cx'],d['lM']) if clump is not None else None
    return QQ,c,-float(r.fun)

# ---------- detection ----------
def screen_point(L,Qb,vb,p,nfev,lM0=8.3):
    """proven screen (found J0946's clump blind): clump starts small; free lM + te,e0,e1 + multipoles + reg. returns (delta,y,x,lM)"""
    from scipy.optimize import minimize
    K=['te','e0','e1']+MP+['reg']; x0=np.array([lM0]+[Qb[k] for k in K])
    st=np.array([0.3]+[S.STEP[k] for k in K]); sim=np.array([x0]+[x0+np.eye(len(x0))[i]*st[i] for i in range(len(x0))])
    def g(x):
        if not 8.0<=x[0]<=10.7: return 1e9
        return -post(L,dict(Qb,**dict(zip(K,x[1:]))),[(p[0],p[1],x[0])])
    r=minimize(g,x0,method='Nelder-Mead',options={'initial_simplex':sim,'maxfev':nfev})
    return (float(-r.fun-vb),p[0],p[1],float(r.x[0]))
def detect(L,Q0,cfg):
    out=dict()
    Qb,_,vb=nm_fit(lambda Q,c:post(L,Q),Q0,MACRO+MP,None,cfg['nf_macro']); out['v_base']=vb
    if cfg.get('screen_given') is not None:   # reuse an earlier screen: [(score,y,x,lM),...]
        scr=sorted([tuple(s) for s in cfg['screen_given']],reverse=True); pts=[]
    else: pts=S.arc_positions(L,Qb,spacing=cfg['spacing'])
    if cfg.get('presel') and pts:
        from scipy.ndimage import gaussian_filter
        f=al.FitImaging(dataset=L.ds,tracer=tracer(L,Qb)); r=np.array(f.normalized_residual_map.native); m=np.array(L.ds.mask)
        r2=gaussian_filter(np.where(m,0,r**2-1),2.0); n=r2.shape[0]; ps=L.ds.pixel_scales[0]; c=(n-1)/2
        sc=[(r2[int(round(c-p[0]/ps)),int(round(c+p[1]/ps))],p) for p in pts]; sc.sort(key=lambda t:-t[0])
        pts=[p for _,p in sc[:cfg['presel']]]
    if cfg.get('screen_given') is None: scr=[]
    for p in pts:
        scr.append(screen_point(L,Qb,vb,p,cfg['nf_screen']))
    scr.sort(reverse=True); out['screen_top']=scr[:5]; out['n_pts']=len(pts) or len(scr)
    # peaks separated by >0.45"
    peaks=[]
    for s in scr:
        if s[0]<cfg['pre_thr']: break
        if all(np.hypot(s[1]-q[1],s[2]-q[2])>0.45 for q in peaks): peaks.append(s)
        if len(peaks)>=cfg['max_peaks']: break
    # equal-budget smooth continuation (shared baseline for all peaks)
    Qs,_,vs=nm_fit(lambda Q,c:post(L,Q),Qb,MACRO+MP,None,cfg['nf_verify']) if peaks else (Qb,None,vb)
    v0=max(vb,vs); Q0b=Qs if vs>vb else Qb
    out['v0']=v0; cands=[]
    for s in peaks:
        Qc,c,vc=nm_fit(lambda Q,c:post(L,Q,[c]),Q0b,MACRO+MP,(s[1],s[2],s[3]),cfg['nf_verify'])
        Qm,cm,vm=nm_fit(lambda Q,c:post(L,Q,[c]),Q0b,MACRO+MP,(-s[1],-s[2],s[3]),cfg['nf_verify'])
        off=float(np.hypot(c[0]-s[1],c[1]-s[2]))
        cands.append(dict(screen=s[0],start=[s[1],s[2],s[3]],clump=list(c),d=vc-v0,mirror=vm-v0,off=off,Qc={k:Qc[k] for k in MACRO+MP},Q0={k:Q0b[k] for k in MACRO+MP}))
    out['cands']=cands; return out

# ---------- validation ----------
def holdout_ds(L,cy,cx,R):
    n0=np.array(L.ds.noise_map.native); n=n0.shape[0]; ps=L.ds.pixel_scales[0]; c=(n-1)/2
    yy,xx=np.mgrid[0:n,0:n]; hold=(np.hypot((c-yy)*ps-cy,(xx-c)*ps-cx)<R)&(~np.array(L.ds.mask))
    nm=n0.copy(); nm[hold]=1e8; d=copy.copy(L.ds); d.noise_map=al.Array2D(values=nm,mask=L.ds.mask); return d,hold,n0
def hide_and_predict(L,P,cand,R,nfev):
    c=cand['clump']; dsh,hold,n0=holdout_ds(L,c[0],c[1],R); data=np.array(L.ds.data.native)
    res={}
    for tag,Qs,cl in (('noclump',dict(P,**cand['Q0']),[]),('clump',dict(P,**cand['Qc']),[c])):
        QQ,_,_=nm_fit(lambda Q,_c:post(L,Q,cl,ds=dsh),Qs,MACRO+MP,None,nfev)
        f=al.FitImaging(dataset=dsh,tracer=tracer(L,QQ,cl)); mi=np.array(f.model_data.native)
        res[tag]=float((((data-mi)/n0)**2)[hold].sum())
    res['n_hold']=int(hold.sum()); res['score']=res['noclump']-res['clump']; return res
def split_test(L,P,cand,halves):
    """halves: two independent half-data images (noise x sqrt2). Macros fixed from the full-data fits."""
    c=cand['clump']; Q0=dict(P,**cand['Q0']); Qc=dict(P,**cand['Qc']); out=[]
    for img in halves:
        d=copy.copy(L.ds); d.data=al.Array2D(values=np.asarray(img),mask=L.ds.mask); d.noise_map=al.Array2D(values=np.array(L.ds.noise_map.native)*np.sqrt(2),mask=L.ds.mask)
        out.append(logev(L,Qc,[c],ds=d)-logev(L,Q0,ds=d))
    return dict(d_half=out,min=min(out))
def lopsided_adversary(L,P,cand,nfev):
    c=cand['clump']; Qa=dict(P,**cand['Q0'],a1c=0.,a1s=0.); Qac=dict(P,**cand['Qc'],a1c=0.,a1s=0.); K=MACRO+MP+['a1c','a1s']
    _,_,v0=nm_fit(lambda Q,_c:post(L,Q,extra=m1_extra(Q)),Qa,K,None,nfev)
    _,_,v1=nm_fit(lambda Q,_c:post(L,Q,[c],extra=m1_extra(Q)),Qac,K,None,nfev)
    return dict(d_m1=v1-v0)
