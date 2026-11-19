import numpy as np, autolens as al
PIX=0.05; ZL=0.222; ZS=0.609; SKY_SHIFT=0.02
_CACHE={}
def inputs(name):
    if name not in _CACHE: _CACHE[name]=dict(np.load(__import__('os').path.join(__import__('os').environ.get('DS_DATA','data'),'%s_inputs.npz'%name)))
    return _CACHE[name]
def set_z(zl,zs):
    global ZL,ZS; ZL,ZS=zl,zs
def make_dataset(mask_in=0.8,mask_out=1.85,half=50,extra=None,psf=None,noise_scale=1.0,os_lp=2,os_pix=2,name='j0946'):
    global PIX
    D=inputs(name); PIX=float(D['pix']) if 'pix' in D else 0.05; h=D['img'].shape[0]//2
    if 'pix' in D: half=int(round(2.5/PIX))
    s=slice(h-half,h+half+1); sh=float(D['sky_shift']) if 'sky_shift' in D else SKY_SHIFT
    f=D['img'][s,s]+sh; r=D['noise'][s,s]*noise_scale; p=D['psf'] if psf is None else psf
    ds=al.Imaging(data=al.Array2D.no_mask(f,pixel_scales=PIX),noise_map=al.Array2D.no_mask(r,pixel_scales=PIX),
                  psf=al.Convolver(kernel=al.Array2D.no_mask(p/p.sum(),pixel_scales=PIX),normalize=True))
    n=f.shape[0]; yy,xx=np.indices(f.shape); rr=np.hypot(xx-(n-1)/2,yy-(n-1)/2)*PIX
    m=(rr>mask_out)|(rr<mask_in)|(~D['valid'][s,s])
    if extra is not None: m|=extra[s,s]
    return ds.apply_mask(mask=al.Mask2D(mask=m,pixel_scales=PIX)).apply_over_sampling(over_sample_size_lp=os_lp,over_sample_size_pixelization=os_pix)
def light_basis(centre,ell1,ell2=None,n=25,smin=0.03,smax=3.0):
    sig=np.logspace(np.log10(smin),np.log10(smax),n)
    prof=[al.lp_linear.Gaussian(centre=centre,ell_comps=ell1,sigma=float(x)) for x in sig]
    if ell2 is not None: prof+=[al.lp_linear.Gaussian(centre=centre,ell_comps=ell2,sigma=float(x)) for x in sig]
    prof.append(al.lp_linear.Gaussian(centre=(0.,0.),ell_comps=(0.,0.),sigma=500.))   # flat sky term
    return al.lp_basis.Basis(profile_list=prof)
def tracer(P,subhalo=None,mesh=(36,36)):
    """P dict: mass centre (my,mx), ell (e0,e1), te, slope, shear (g1,g2); light centre/ells; reg"""
    mass=al.mp.PowerLaw(centre=(P['my'],P['mx']),ell_comps=(P['e0'],P['e1']),einstein_radius=P['te'],slope=P['slope'])
    lens=al.Galaxy(redshift=ZL,bulge=light_basis((P['ly'],P['lx']),(P['a0'],P['a1']),(P['b0'],P['b1'])),mass=mass,shear=al.mp.ExternalShear(gamma_1=P['g1'],gamma_2=P['g2']))
    gals=[lens]
    if subhalo is not None: gals.append(subhalo)
    src=al.Galaxy(redshift=ZS,pixelization=al.Pixelization(mesh=al.mesh.RectangularUniform(shape=mesh),regularization=al.reg.Constant(coefficient=P['reg'])))
    return al.Tracer(galaxies=gals+[src])
P0=dict(my=0.,mx=0.,e0=0.05,e1=0.0,te=1.4,slope=2.0,g1=0.,g2=0.,ly=0.,lx=0.,a0=0.01,a1=-0.02,b0=0.1,b1=0.1,reg=300.)
if __name__=='__main__':
    import time
    ds=make_dataset(); print('mask pixels',ds.mask.pixels_in_mask)
    t=time.time(); fit=al.FitImaging(dataset=ds,tracer=tracer(P0)); print('eval time %.2fs'%(time.time()-t),'chi2_red %.2f'%(fit.chi_squared/ds.mask.pixels_in_mask),'logev %.0f'%fit.log_evidence)
    t=time.time(); fit=al.FitImaging(dataset=ds,tracer=tracer(P0)); print('eval time %.2fs'%(time.time()-t))
