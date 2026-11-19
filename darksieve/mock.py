"""Stage 3: per-lens mock calibration.
1) fit a smooth parametric 'truth' source (linear Gaussian basis) under the best no-clump lens model
2) mocks = truth image (+ optional injected clump) + correlated noise matching the drizzled data
3) run the same screen/verify on mocks -> false-alarm threshold & completeness vs mass"""
import sys, json, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import autolens as al
from scipy.optimize import minimize
from scipy.ndimage import gaussian_filter
from qload import load; from model import Lens; import stages as S, pj, hst
def src_basis(cy,cx,e0,e1,n=24,smin=0.01,smax=0.6):
    sig=np.logspace(np.log10(smin),np.log10(smax),n)
    return al.lp_basis.Basis(profile_list=[al.lp_linear.Gaussian(centre=(cy,cx),ell_comps=(e0,e1),sigma=float(s)) for s in sig]+
                                          [al.lp_linear.Gaussian(centre=(cy,cx),ell_comps=(0.,0.),sigma=float(s)) for s in sig[::2]])
def truth_tracer(L,Q,p,clump=None):
    t=L.tracer(Q,clump); g=list(t.galaxies); g[-1]=al.Galaxy(redshift=L.zs,bulge=src_basis(*p)); return al.Tracer(galaxies=g)
def fit_truth(L,Q,nfev=150):
    def f(p):
        try: return -al.FitImaging(dataset=L.ds,tracer=truth_tracer(L,Q,p)).figure_of_merit
        except Exception: return 1e9
    best=None
    for p0 in [(0,0,0,0),(0.1,0.1,0,0),(-0.1,-0.1,0,0),(0.1,-0.1,0,0),(-0.1,0.1,0,0)]:
        r=minimize(f,p0,method='Nelder-Mead',options={'maxfev':nfev//5*2});
        if best is None or r.fun<best.fun: best=r
    r=minimize(f,best.x,method='Nelder-Mead',options={'maxfev':nfev}); fit=al.FitImaging(dataset=L.ds,tracer=truth_tracer(L,Q,r.x))
    fixed=fit.model_obj_linear_light_profiles_to_light_profiles
    return r.x,fixed,float(fit.chi_squared/L.ds.mask.pixels_in_mask)
def model_image(L,fixed,clump=None,Q=None):
    g=list(fixed.galaxies)
    if clump is not None: g=[g[0],pj.pj_galaxy(clump[0],clump[1],10**clump[2],L.zl,L.zs,Q['te']),g[-1]]
    f=al.FitImaging(dataset=L.ds,tracer=al.Tracer(galaxies=g)); return np.array(f.model_data.native)
def noise(L,F,seed):
    rng=np.random.default_rng(seed); s=np.array(L.ds.noise_map.native)/np.sqrt(F)
    w=rng.normal(size=s.shape); k=gaussian_filter(w,0.437); k*=1/np.std(gaussian_filter(rng.normal(size=(400,400)),0.437))
    return k*s
if __name__=='__main__':
    cfg=json.load(open(sys.argv[1])); P,Q=load(cfg['start'],cfg['fit']); L=Lens(cfg['name'],cfg['zl'],cfg['zs'],P)
    p,fixed,c2=fit_truth(L,Q); print('truth source',np.round(p,3),'chi2_red(white)',round(c2,3),flush=True)
    import pickle; pickle.dump(dict(p=list(p),Q=Q),open(cfg['out'].replace('.json','_truth.pkl'),'wb'))
    base=model_image(L,fixed); np.save(cfg['out'].replace('.json','_truth_img.npy'),base)
    F=float(hst.inputs(cfg['name'])['F']); rng=np.random.default_rng(0); out=[]
    plan=[(None,s) for s in range(cfg.get('nnull',4))]
    for lm in cfg.get('masses',[9.0,9.5]):
        for k in range(cfg.get('ninj',4)):
            th=rng.uniform(0,2*np.pi); r=Q['te']*rng.uniform(0.92,1.08); plan.append(((r*np.sin(th),r*np.cos(th),lm),100+len(plan)))
    for i,(cl,seed) in enumerate(plan):
        img=model_image(L,fixed,cl,Q)+noise(L,F,seed); fn=cfg['out'].replace('.json','_mock%02d.npy'%i); np.save(fn,img)
        out.append(dict(i=i,clump=cl,file=fn))
    json.dump(out,open(cfg['out'],'w'),indent=1); print('mocks',len(out),flush=True)
