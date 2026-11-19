"""Build the smooth lens model for a new lens (no subhalo).  usage: prepare.py NAME zl zs thetaE outdir
Writes outdir/<NAME>_start.json, <NAME>_smooth.json and configs-style blind_<NAME>.json."""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
from scipy.ndimage import gaussian_filter
import isophote2 as I, hst
from model import Lens, MP, MACRO; import stages as S
name,zl,zs,teL=sys.argv[1],float(sys.argv[2]),float(sys.argv[3]),float(sys.argv[4]); od=sys.argv[5]; os.makedirs(od,exist_ok=True)
Dn=hst.inputs(name); pix=float(Dn['pix']) if 'pix' in Dn else 0.05; c=Dn['img'].shape[0]//2
img=gaussian_filter(Dn['img'],0.6); p=[0.85,0.3,float(c),float(c)]
for a in (0.3,0.4,0.5): p=I.fit(img,a/pix,p)
q,pa,cy,cx=p; ly=(c-cy)*pix; lx=(cx-c)*pix; phi=np.pi-pa; f=(1-q)/(1+q); ell=(f*np.sin(2*phi),f*np.cos(2*phi))
if np.hypot(ly,lx)>0.25 or np.hypot(*ell)>0.35: ly=lx=0.0; ell=(0.,0.)
P=dict(my=ly,mx=lx,e0=ell[0],e1=ell[1],te=teL,slope=2.0,g1=0.,g2=0.,ly=ly,lx=lx,a0=ell[0],a1=ell[1],b0=ell[0],b1=ell[1],reg=1.0,
       a3c=0.,a3s=0.,a4c=0.,a4s=0.,mask_in=max(0.4,teL-0.6),mask_out=min(2.45,teL+0.6),mesh=[30,30],os=1)
json.dump(P,open(f'{od}/{name}_start.json','w'),indent=1); print('start',name,'q %.2f te %.2f'%(q,teL),flush=True)
L=Lens(name,zl,zs,P); Q=dict(P)
for it in range(3):
    qq,v=S.smooth_fit(L,Q,0.01,250); Q.update(qq); fit=L.fit(Q)
    sm=dict(keys=MACRO+MP,x=[float(Q[k]) for k in MACRO+MP],logpost=float(v),nev=250*(it+1),chi2=float(fit.chi_squared/L.ds.mask.pixels_in_mask))
    json.dump(sm,open(f'{od}/{name}_smooth.json','w'),indent=1); print('smooth',sm['nev'],'logpost %.1f chi2 %.3f'%(v,sm['chi2']),flush=True)
cfg=dict(name=name,zl=zl,zs=zs,start=f'{od}/{name}_start.json',fit=f'{od}/{name}_smooth.json',fast=True,
         det=dict(nf_macro=150,spacing=0.3,pre_thr=2.0,max_peaks=3,nf_verify=150,nf_screen=80,presel=0),R=0.2,nf_hp=100,nf_adv=120,chi2=sm['chi2'])
json.dump(cfg,open(f'{od}/blind_{name}.json','w'),indent=1)
