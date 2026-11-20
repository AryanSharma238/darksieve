"""DarkSieve false-positive emulator.
Per lens, precompute (once, ~1500 pixelized-source fits):
  r(data) = noise-weighted residual after re-solving the pixelized source at the smooth model Qt (linear in data)
  J        = d r / d theta, theta = smooth-model parameters the search refits      -> P_perp projects them out
  m[p,M]   = residual change from a subhalo of mass M at blind-grid position p    (subhalo templates)
  e_k      = residual from unit-amplitude lens-model errors (m=1, m=3, m=4; cos/sin)
  G        = m_perp . r(noise) for N noise realizations (correlated HST noise, same generator as the benchmark)
Then any blind search is emulated with linear algebra only:
  gain[p,M] = m_perp.(E a + n) - |m_perp|^2/2 ;  d_pred = max_{p,M} gain - c   (c = evidence offset, the ONE calibrated number)
usage: emulator.py build CFG out.npz [nnoise]"""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
def setup(cfgp):
    g={'__file__':os.path.abspath(__file__)}; sys.argv=['bench.py',cfgp,'0','0','/tmp/none.jsonl']
    src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'bench.py')).read()
    exec(src[:src.index('rng=np.random.default_rng(7000+wid)')],g); return g
if __name__=='__main__' and sys.argv[1]=='build':
    cfgp,outp=sys.argv[2],sys.argv[3]; NN=int(sys.argv[4]) if len(sys.argv)>4 else 200
    g=setup(cfgp); L,Qt,data0,ds,al,S,noise,F,srcG=(g[k] for k in ('L','Qt','data0','ds','al','S','noise','F','srcG'))
    from model import MACRO, MP
    mask=~np.array(data0.mask); sig=np.array(data0.noise_map.native)[mask]
    def truth(extra=None):
        t=ds.tracer(L,Qt,[],extra); gg=list(t.galaxies)[:-1]+[srcG]
        return np.array(al.FitImaging(dataset=data0,tracer=al.Tracer(galaxies=gg)).model_data.native)
    def resid(data,Q=None,cl=()):
        L.ds=data0; L.set_data(data); f=al.FitImaging(dataset=L.ds,tracer=ds.tracer(L,Q or Qt,list(cl)))
        return (np.array(L.ds.data.native)-np.array(f.model_data.native))[mask]/sig
    T0=truth(); r0=resid(T0)
    H=dict(te=0.003,slope=0.01,my=0.003,mx=0.003,e0=0.003,e1=0.003,g1=0.003,g2=0.003,a3c=0.002,a3s=0.002,a4c=0.002,a4s=0.002)
    keys=[k for k in MACRO+MP if k in Qt and k!='reg']
    J=np.array([(resid(T0,dict(Qt,**{k:Qt[k]+H.get(k,0.003)}))-r0)/H.get(k,0.003) for k in keys]).T
    Qj,_=np.linalg.qr(J); P=lambda v: v-Qj@(Qj.T@v); print('J',J.shape,flush=True)
    pts=S.arc_positions(L,Qt,spacing=0.15); pts=[p for p in pts if all(np.hypot(p[0]-c[0],p[1]-c[1])>0.3 for c in ds.BASE)]
    masses=[8.0,8.25,8.5,8.75,9.0,9.25,9.5,9.75,10.0]
    Mp=np.zeros((len(pts),len(masses),len(r0)),np.float32)
    for i,p in enumerate(pts):
        for j,lm in enumerate(masses): Mp[i,j]=P(r0-resid(T0,None,[(p[0]+0.0013,p[1]+0.0017,lm)]))   # model change from the subhalo, projected
        if i%10==0: print('templates',i,len(pts),flush=True)
    E={}
    for m,slope in ((1,2.05),(3,2.0),(4,2.0)):
        for nm,comps in (('c',(0.01,0.)),('s',(0.,0.01))):
            ex=dict(mpx=al.mp.PowerLawMultipole(m=m,centre=(Qt['my'],Qt['mx']),einstein_radius=Qt['te'],slope=slope,multipole_comps=comps))
            E['m%d%s'%(m,nm)]=P(resid(truth(ex))-r0)/0.01
    print('errors',list(E),flush=True)
    Mf=Mp.reshape(-1,len(r0)).astype(np.float64); G=np.zeros((NN,Mf.shape[0]))
    for s in range(NN):
        n=np.asarray(noise(L,F,987654+s)); G[s]=Mf@P(resid(T0+n)-r0)
        if s%25==0: print('noise',s,flush=True)
    np.savez_compressed(outp,pts=np.array(pts),masses=np.array(masses),half=0.5*(Mf**2).sum(1),G=G,
        **{'E_'+k:Mf@v for k,v in E.items()},Eabs=np.array([np.sqrt((v**2).sum()) for v in E.values()]),keys=np.array(list(E)))
    print('saved',outp,'templates',Mf.shape,flush=True)
