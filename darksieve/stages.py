import numpy as np, json
from scipy.optimize import minimize
from model import Lens, logprior, MACRO, MP
BND=dict(my=(-.3,.3),mx=(-.3,.3),e0=(-.4,.4),e1=(-.4,.4),te=(.5,2.5),g1=(-.15,.15),g2=(-.15,.15),reg=(-1,3),a3c=(-.06,.06),a3s=(-.06,.06),a4c=(-.06,.06),a4s=(-.06,.06),cy=(-3,3),cx=(-3,3),lM=(8,10.7))
STEP=dict(my=.02,mx=.02,e0=.03,e1=.03,te=.02,g1=.02,g2=.02,reg=.2,a3c=.005,a3s=.005,a4c=.005,a4s=.005,cy=.05,cx=.05,lM=.25)
def nm(f,x0,keys,nfev):
    x0=np.array(x0,float); sim=np.array([x0]+[x0+np.eye(len(x0))[i]*STEP[keys[i]] for i in range(len(x0))])
    def g(x):
        if any(not BND[k][0]<=v<=BND[k][1] for k,v in zip(keys,x)): return 1e9
        return -f(dict(zip(keys,x)))
    r=minimize(g,x0,method='Nelder-Mead',options={'initial_simplex':sim,'maxfev':nfev}); return dict(zip(keys,r.x)),-r.fun
def post(L,Q,clump,sig): return L.logev(Q,clump)+logprior(Q,sig)
def smooth_fit(L,Q0,sig,nfev=600):
    keys=MACRO+MP; return nm(lambda q:post(L,dict(Q0,**q),None,sig),[Q0[k] for k in keys],keys,nfev)
def screen_point(L,Q,pos,sig,nfev=40,lM0=9.3):
    """partial refit at one position: clump mass + multipoles + reg"""
    keys=['lM']+MP+['reg']; x0=[lM0]+[Q[k] for k in MP]+[Q['reg']]
    def f(q):
        QQ=dict(Q,**{k:q[k] for k in MP+['reg']}); return post(L,QQ,(pos[0],pos[1],q['lM']),sig)
    q,v=nm(f,x0,keys,nfev); return v,q['lM']
def joint_fit(L,Q,clump,sig,nfev=600):
    keys=MACRO+MP+['cy','cx','lM']; x0=[Q[k] for k in MACRO+MP]+list(clump)
    def f(q): return post(L,dict(Q,**{k:q[k] for k in MACRO+MP}),(q['cy'],q['cx'],q['lM']),sig)
    return nm(f,x0,keys,nfev)
def arc_positions(L,Q,spacing=0.2,frac=0.15):
    """grid points inside mask where the lensed-source model is bright"""
    f=L.fit(Q); src=np.array(f.model_images_of_planes_list[-1].native); m=np.array(L.ds.mask)
    n=src.shape[0]; ps=0.05; c=(n-1)/2; pts=[]
    thr=frac*np.nanmax(np.where(m,np.nan,src))
    from scipy.ndimage import maximum_filter
    sm=maximum_filter(np.where(m,0,src),size=int(round(spacing/ps)))
    for y in np.arange(-2.5,2.51,spacing):
        for x in np.arange(-2.5,2.51,spacing):
            i=int(round(c-y/ps)); j=int(round(c+x/ps))
            if 0<=i<n and 0<=j<n and not m[i,j] and sm[i,j]>thr: pts.append((round(y,3),round(x,3)))
    return pts
def screen_point2(L,Q,pos,sig,nfev=80,lM0=9.5):
    """partial refit: clump mass + Einstein radius, ellipticity, shear, multipoles, reg (centre fixed)"""
    free=['te','e0','e1','g1','g2']+MP+['reg']; keys=['lM']+free; x0=[lM0]+[Q[k] for k in free]
    def f(q): return post(L,dict(Q,**{k:q[k] for k in free}),(pos[0],pos[1],q['lM']),sig)
    q,v=nm(f,x0,keys,nfev); return v,q['lM']
def screen_point3(L,Q,pos,sig,nfev=120,lM0=8.3):
    """partial refit starting from (near) the no-clump solution so delta >= ~0"""
    free=['te','e0','e1']+MP+['reg']; keys=['lM']+free; x0=[lM0]+[Q[k] for k in free]
    def f(q): return post(L,dict(Q,**{k:q[k] for k in free}),(pos[0],pos[1],q['lM']),sig)
    q,v=nm(f,x0,keys,nfev); return v,q['lM']
# ---- iterative (multi-clump) support ----
def post_multi(L,Q,clumps,sig):
    t=L.tracer(Q,None); g=list(t.galaxies)
    import pj
    for c in clumps: g.insert(1,pj.pj_galaxy(c[0]+0.0013,c[1]+0.0017,10**c[2],L.zl,L.zs,Q['te']))
    import autolens as al, numpy as np
    try:
        v=al.FitImaging(dataset=L.ds,tracer=al.Tracer(galaxies=g)).log_evidence; v=v if np.isfinite(v) else -1e9
    except Exception: v=-1e9
    return v+logprior(Q,sig)
def joint_fit_multi(L,Q,fixed_clumps,new_clump,sig,nfev=500):
    """refit macro + all clumps (fixed_clumps are refit too); new_clump may be None"""
    keys=MACRO+MP; x0=[Q[k] for k in keys]; allc=list(fixed_clumps)+([new_clump] if new_clump is not None else [])
    for i,c in enumerate(allc): keys=keys+['cy','cx','lM']; x0+=list(c)
    names=MACRO+MP
    def f(q):
        xs=[q[k] for k in names]  # dict keys collide for multiple clumps -> use positional below
        return 0
    import numpy as np
    from scipy.optimize import minimize
    x0=np.array(x0,float); steps=[STEP[k] for k in keys]
    sim=np.array([x0]+[x0+np.eye(len(x0))[i]*steps[i] for i in range(len(x0))])
    def g(x):
        if any(not BND[k][0]<=v<=BND[k][1] for k,v in zip(keys,x)): return 1e9
        QQ=dict(Q,**dict(zip(names,x[:len(names)]))); cl=[tuple(x[len(names)+3*i:len(names)+3*i+3]) for i in range(len(allc))]
        return -post_multi(L,QQ,cl,sig)
    r=minimize(g,x0,method='Nelder-Mead',options={'initial_simplex':sim,'maxfev':nfev})
    QQ=dict(Q,**dict(zip(names,r.x[:len(names)]))); cl=[tuple(map(float,r.x[len(names)+3*i:len(names)+3*i+3])) for i in range(len(allc))]
    return QQ,cl,-r.fun

M2=['m2e0','m2e1','m2f']
BND.update(m2e0=(-.5,.5),m2e1=(-.5,.5),m2f=(0.02,0.9)); STEP.update(m2e0=.05,m2e1=.05,m2f=.1)
def fit_m2(L,Q,clump,sig,nfev):
    Q=dict(Q,use_m2=1); keys=MACRO+MP+M2+(['cy','cx','lM'] if clump else []); x0=[Q.get(k,{'m2e0':Q['e0'],'m2e1':Q['e1'],'m2f':0.3}.get(k)) for k in MACRO+MP+M2]+(list(clump) if clump else [])
    def f(q):
        QQ=dict(Q,**{k:q[k] for k in MACRO+MP+M2}); c=(q['cy'],q['cx'],q['lM']) if clump else None
        return post(L,QQ,c,sig)
    return nm(f,x0,keys,nfev)
