"""Shared lens model: PowerLaw+shear (+m3,m4 multipoles), MGE lens light, adaptive pixel source, optional pseudo-Jaffe clump."""
import numpy as np, json, sys
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import autolens as al, hst, pj
MACRO=['my','mx','e0','e1','te','g1','g2','reg']; MP=['a3c','a3s','a4c','a4s']
class Lens:
    def __init__(s,name,zl,zs,P,mask_in=None,mask_out=None,data=None,mesh=None,os_=None):
        mask_in=P.get('mask_in',0.8) if mask_in is None else mask_in; mask_out=P.get('mask_out',1.85) if mask_out is None else mask_out
        hst.set_z(zl,zs); s.name,s.zl,s.zs,s.P=name,zl,zs,dict(P)
        mesh=tuple(P.get('mesh',(40,40))) if mesh is None else mesh; os_=P.get('os',2) if os_ is None else os_
        s.ds=hst.make_dataset(name=name,mask_in=mask_in,mask_out=mask_out,os_lp=os_,os_pix=os_)
        s.mesh=mesh
        if data is not None: s.set_data(data)
    def set_data(s,arr):
        import copy; d=copy.copy(s.ds); d.data=al.Array2D(values=np.asarray(arr),mask=s.ds.mask); s.ds=d
    def tracer(s,Q,clump=None,mp=True,extra=None):
        Q=dict(s.P,**Q); t=hst.tracer(dict(Q,reg=10**Q['reg'])); L=t.galaxies[0]
        kw=dict(redshift=L.redshift,bulge=L.bulge,mass=L.mass,shear=L.shear)
        if mp:
            for m in (3,4): kw['mp%d'%m]=al.mp.PowerLawMultipole(m=m,centre=(Q['my'],Q['mx']),einstein_radius=Q['te'],slope=Q['slope'],multipole_comps=(Q['a%dc'%m],Q['a%ds'%m]))
        if 'm2f' in Q and Q.get('use_m2',0):
            f=Q['m2f']; kw['mass']=al.mp.PowerLaw(centre=(Q['my'],Q['mx']),ell_comps=(Q['e0'],Q['e1']),einstein_radius=Q['te']*max(1-f,1e-3),slope=Q['slope'])
            kw['mass2']=al.mp.Isothermal(centre=(Q['my'],Q['mx']),ell_comps=(Q['m2e0'],Q['m2e1']),einstein_radius=Q['te']*max(f,1e-3))
        if extra: kw.update(extra)
        g=[al.Galaxy(**kw)]
        if clump is not None: g.append(pj.pj_galaxy(clump[0]+0.0013,clump[1]+0.0017,10**clump[2],s.zl,s.zs,Q['te']))
        g.append(al.Galaxy(redshift=s.zs,pixelization=al.Pixelization(mesh=al.mesh.RectangularBilinearAdaptDensity(shape=s.mesh),regularization=al.reg.Constant(coefficient=10**Q['reg']))))
        return al.Tracer(galaxies=g)
    def logev(s,Q,clump=None,mp=True):
        try:
            v=al.FitImaging(dataset=s.ds,tracer=s.tracer(Q,clump,mp)).log_evidence
            return v if np.isfinite(v) else -1e9
        except Exception: return -1e9
    def fit(s,Q,clump=None,mp=True): return al.FitImaging(dataset=s.ds,tracer=s.tracer(Q,clump,mp))
MP_CENTER={}
def logprior(Q,sig):
    return 0. if not sig else -0.5*sum((Q[k]-MP_CENTER.get(k,0.))**2 for k in MP)/sig**2
