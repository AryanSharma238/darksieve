"""DarkSieve benchmark worker: build realistic + adversarial fakes of one real lens, run detection + validation, log everything.
usage: bench.py LENSCFG.json worker_id n_trials out.jsonl
Fake classes: null, null_m1 (lopsided 1-2%), null_m34 (extra boxy/triangular 1-2%), null_knots (clumpy source),
              inj (pseudo-Jaffe clump logM 8.5-10 on the arcs; 30% also get knots or m1)
Each fake has independent half-data images (noise x sqrt2) so the split test is exact: full = mean of halves."""
import sys, json, os, pickle, numpy as np
sys.path.insert(0,__import__('os').path.dirname(__import__('os').path.abspath(__file__)))
import autolens as al, pj, hst
from model import Lens, MACRO, MP
import stages as S, ds
from mock import fit_truth, noise
cfg=json.load(open(sys.argv[1])); wid=int(sys.argv[2]); ntr=int(sys.argv[3]); out=sys.argv[4]
if 'verify' in cfg:
    P=json.load(open(cfg['start'])); P.pop('_x',None); V=json.load(open(cfg['verify'])); Qt=dict(P,**V['macro'][cfg.get('macro_key','base_0')])
else:
    from qload import load; P,Qt=load(cfg['start'],cfg['fit'])
fast=cfg.get('fast',False)
L=Lens(cfg['name'],cfg['zl'],cfg['zs'],P,mesh=(30,30) if fast else None,os_=1 if fast else None); data0=L.ds; F=float(hst.inputs(cfg['name']).get('F',1.0))
ds.BASE=[tuple(c) for c in cfg.get('base_clumps',[])]
tp=cfg['name'][:12]+'_truthsrc.pkl'
if 'src_params' in cfg:
    from mock import src_basis
    sp=json.load(open(cfg['src_params'])); t=ds.tracer(L,Qt); g=list(t.galaxies); g[-1]=al.Galaxy(redshift=L.zs,bulge=src_basis(*sp))
    fixed=al.FitImaging(dataset=data0,tracer=al.Tracer(galaxies=g)).model_obj_linear_light_profiles_to_light_profiles
elif os.path.exists(tp): fixed=pickle.load(open(tp,'rb'))
else:
    p,fixed,c2=fit_truth(L,Qt); pickle.dump(fixed,open(tp,'wb')); print('truth source fitted chi2',round(c2,3),flush=True)
G=list(fixed.galaxies); lensG,srcG,baseG=G[0],G[-1],G[1:-1]
pts=json.load(open(cfg['arc_pts'])) if 'arc_pts' in cfg else S.arc_positions(L,Qt,spacing=0.15)
pts=[p for p in pts if all(np.hypot(p[0]-c[0],p[1]-c[1])>0.4 for c in ds.BASE)]; print('arc pts',len(pts),flush=True)
# source peak for knot amplitudes
sg=al.Grid2D.uniform(shape_native=(41,41),pixel_scales=0.02); S0=float(np.max(np.array(srcG.image_2d_from(grid=sg))))
def make(kind,rng):
    lg=lensG; extra={}; meta=dict(kind=kind)
    if kind in ('null_m1',) or (kind=='inj' and rng.random()<0.15):
        a=rng.uniform(0.01,0.02); th=rng.uniform(0,2*np.pi); meta['m1']=[a*np.cos(th),a*np.sin(th)]
        extra['mpx1']=al.mp.PowerLawMultipole(m=1,centre=(Qt['my'],Qt['mx']),einstein_radius=Qt['te'],slope=2.05,multipole_comps=tuple(meta['m1']))
    if kind=='null_m34':
        m=int(rng.choice([3,4])); a=rng.uniform(0.01,0.02); th=rng.uniform(0,2*np.pi); meta['m34']=[m,a*np.cos(th),a*np.sin(th)]
        extra['mpx34']=al.mp.PowerLawMultipole(m=m,centre=(Qt['my'],Qt['mx']),einstein_radius=Qt['te'],slope=2.0,multipole_comps=tuple(meta['m34'][1:]))
    if extra:
        kw={k:v for k,v in lg.__dict__.items() if k not in ('redshift',) and not k.startswith('_')}
        kw=dict(redshift=lg.redshift,**{k:v for k,v in kw.items() if hasattr(v,'deflections_yx_2d_from') or hasattr(v,'image_2d_from')}); kw.update(extra); lg=al.Galaxy(**kw)
    sgal=srcG
    if kind=='null_knots' or (kind=='inj' and rng.random()<0.15):
        kn={}
        for i in range(int(rng.integers(2,5))):
            kn['knot%d'%i]=al.lp.Gaussian(centre=(float(rng.uniform(-.15,.15)),float(rng.uniform(-.15,.15))),ell_comps=(0.,0.),intensity=float(rng.uniform(0.2,0.6)*S0),sigma=float(rng.uniform(0.02,0.05)))
        meta['knots']=len(kn); sgal=al.Galaxy(redshift=srcG.redshift,bulge=srcG.bulge,**kn)
    gal=[lg]+list(baseG)
    if kind=='inj':
        p=pts[rng.integers(len(pts))]; y=p[0]+rng.uniform(-.07,.07); x=p[1]+rng.uniform(-.07,.07); lm=float(rng.choice(cfg['masses']))
        meta['clump']=[float(y),float(x),lm]; gal.append(pj.pj_galaxy(y+0.0013,x+0.0017,10**lm,L.zl,L.zs,Qt['te']))
    gal.append(sgal)
    truth=np.array(al.FitImaging(dataset=data0,tracer=al.Tracer(galaxies=gal)).model_data.native)
    return truth,meta
rng=np.random.default_rng(7000+wid); kinds=cfg['kinds']; probs=cfg['probs']
done=sum(1 for l in open(out) if json.loads(l).get('w')==wid) if os.path.exists(out) else 0
import time as _t
_dl=float(os.environ.get('DS_DEADLINE','0') or 0); _last=float(os.environ.get('DS_EST','0') or 0)
for k in range(ntr):
    kind=str(rng.choice(kinds,p=probs)); truth,meta=make(kind,rng)
    n1=noise(L,F,10**6*wid+2*k); n2=noise(L,F,10**6*wid+2*k+1)
    if k<done: continue
    if _dl and _t.time()+1.15*_last>_dl: print('STOP: not enough time left for another fake',flush=True); break
    _t0=_t.time()
    halves=[truth+np.sqrt(2)*n1,truth+np.sqrt(2)*n2]; full=(halves[0]+halves[1])/2
    L.ds=data0; L.set_data(full)
    Q0=dict(Qt);
    for kk,sd in (('te',0.01),('e0',0.02),('e1',0.02),('g1',0.01),('g2',0.01)): Q0[kk]=Qt[kk]+rng.normal(0,sd)
    rec=dict(w=wid,k=k,**meta)
    if cfg.get('mode')=='targeted':
        # verify + validate at ONE position: the true clump (inj) or a random arc position (all other classes)
        try:
            from model import MACRO, MP
            Qb,_,vb=ds.nm_fit(lambda Q,c_:ds.post(L,Q),Q0,MACRO+MP,None,cfg['det']['nf_macro'])
            p=meta['clump'][:2] if kind=='inj' else list(pts[rng.integers(len(pts))])
            sc={m:ds.post(L,Qb,[(p[0],p[1],m)]) for m in (8.5,9.0,9.5)}; lm0=max(sc,key=sc.get)
            Qs,_,vs=ds.nm_fit(lambda Q,c_:ds.post(L,Q),Qb,MACRO+MP,None,cfg['det']['nf_verify']); v0=max(vb,vs); Q0b=Qs if vs>vb else Qb
            Qc,c,vc=ds.nm_fit(lambda Q,c_:ds.post(L,Q,[c_]),Q0b,MACRO+MP,(p[0],p[1],lm0),cfg['det']['nf_verify'])
            Qm,cm,vm=ds.nm_fit(lambda Q,c_:ds.post(L,Q,[c_]),Q0b,MACRO+MP,(-p[0],-p[1],lm0),cfg['det']['nf_verify'])
            cand=dict(start=[p[0],p[1],lm0],clump=list(c),d=vc-v0,mirror=vm-v0,off=float(np.hypot(c[0]-p[0],c[1]-p[1])),Qc={k:Qc[k] for k in MACRO+MP},Q0={k:Q0b[k] for k in MACRO+MP})
            cc={kk:cand[kk] for kk in ('start','clump','d','mirror','off')}
            cc['hp']=ds.hide_and_predict(L,P,cand,cfg['R'],cfg['nf_hp'])
            cc['split']=ds.split_test(L,P,cand,halves)
            cc['adv']=ds.lopsided_adversary(L,P,cand,cfg['nf_adv'])
            rec.update(mode='targeted',v0=v0,cands=[cc])
        except Exception as e:
            rec['error']=repr(e)[:300]
        open(out,'a').write(json.dumps(rec)+'\n'); c0=rec.get('cands',[{}])[0]
        print('FAKE',kind,meta.get('clump'),'-> d %.1f mir %.1f hp %.1f split %s adv %.1f'%(c0.get('d',np.nan),c0.get('mirror',np.nan),c0.get('hp',{}).get('score',np.nan),[round(v,1) for v in c0.get('split',{}).get('d_half',[])],c0.get('adv',{}).get('d_m1',np.nan)),rec.get('error',''),flush=True); _last=max(_last,_t.time()-_t0)
        continue
    try:
        det=ds.detect(L,Q0,cfg['det'])
        rec.update(v0=det['v0'],n_pts=det['n_pts'],screen_top=det['screen_top'][:3],cands=[])
        for c in det['cands']:
            cc={kk:c[kk] for kk in ('screen','start','clump','d','mirror','off')}
            if c['d']>cfg['val_dmin'] and c['mirror']<cfg['val_mirror'] and c['off']<0.3:
                cc['hp']=ds.hide_and_predict(L,P,c,cfg['R'],cfg['nf_hp'])
                cc['split']=ds.split_test(L,P,c,halves)
                cc['adv']=ds.lopsided_adversary(L,P,c,cfg['nf_adv'])
            rec['cands'].append(cc)
    except Exception as e:
        rec['error']=repr(e)[:300]
    open(out,'a').write(json.dumps(rec)+'\n'); print('FAKE',json.dumps({kk:rec.get(kk) for kk in ('k','kind','clump','m1','knots')}),'->',[(round(c['d'],1),round(c['mirror'],1),round(c.get('hp',{}).get('score',np.nan),1),[round(v,1) for v in c.get('split',{}).get('d_half',[])],round(c.get('adv',{}).get('d_m1',np.nan),1)) for c in rec.get('cands',[])],rec.get('error',''),flush=True)
    _last=max(_last,_t.time()-_t0)
