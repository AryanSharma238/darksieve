"""Linear gravitational imaging (Koopmans 2005; Vegetti & Koopmans 2009) as a candidate PROPOSER.
Residual of the smooth model is modelled, to first order, as  r = -B [ grad_beta(s) . grad(dpsi) ]
where s is the lensed-source model, B the PSF blur, dpsi a potential correction on a coarse grid.
Solve regularised least squares for dpsi, then dkappa = 0.5 * laplacian(dpsi). Positive dkappa peaks near the arcs = candidates."""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.ndimage import gaussian_filter, maximum_filter
import autolens as al

def _lensed_source_and_jacobian(fit, tracer):
    # lensed source image (PSF-convolved model of the last plane) on the native grid
    src=np.array(fit.model_images_of_planes_list[-1].native)
    n=src.shape[0]
    # Jacobian of the lens mapping by finite differences of deflections on a uniform grid
    ps=fit.dataset.pixel_scales[0]
    g=al.Grid2D.uniform(shape_native=(n,n),pixel_scales=ps)
    a=np.array(tracer.deflections_yx_2d_from(grid=g)).reshape(n,n,2)   # (y,x) components; native rows go from +y (top) down
    # derivatives w.r.t. image coords: row index increases -> y decreases
    day_dy=-np.gradient(a[...,0],ps,axis=0); day_dx=np.gradient(a[...,0],ps,axis=1)
    dax_dy=-np.gradient(a[...,1],ps,axis=0); dax_dx=np.gradient(a[...,1],ps,axis=1)
    A11=1-day_dy; A12=-day_dx; A21=-dax_dy; A22=1-dax_dx       # A = d beta / d theta  (rows: beta_y, beta_x)
    return src,(A11,A12,A21,A22)

def propose(fit, tracer, mask, noise, step=3, lam=None, kmax=4, min_sep=0.4):
    """returns list of (score, y, x) candidate positions and the dkappa map"""
    src,(A11,A12,A21,A22)=_lensed_source_and_jacobian(fit,tracer)
    n=src.shape[0]; ps=fit.dataset.pixel_scales[0]
    m=~np.array(mask)
    resid=np.nan_to_num(np.array(fit.residual_map.native))
    # image-plane gradient of the lensed source, then source-plane gradient grad_beta s = A^{-T} grad_theta I
    Iy=-np.gradient(src,ps,axis=0); Ix=np.gradient(src,ps,axis=1)
    det=A11*A22-A12*A21; det=np.where(np.abs(det)<1e-3,np.sign(det)*1e-3+1e-12,det)
    # A^{-T} = (1/det) [[A22,-A21],[-A12,A11]]
    gy=( A22*Iy - A21*Ix)/det; gx=(-A12*Iy + A11*Ix)/det
    # unknown dpsi on a coarse grid (every `step` pixels), bilinear interpolation to pixels
    nc=(n-1)//step+1; cy=np.arange(nc)*step;
    idx=np.argwhere(m); N=len(idx);
    rows=[];cols=[];vals_y=[];vals_x=[]
    # build d(dpsi)/dy and d(dpsi)/dx at each masked pixel via bilinear gradient of the coarse grid
    for k,(i,j) in enumerate(idx):
        fi=i/step; fj=j/step; i0=min(int(fi),nc-2); j0=min(int(fj),nc-2); ti=fi-i0; tj=fj-j0
        h=step*ps
        # psi(i,j) = sum w * psi_c ; d/drow and d/dcol of the bilinear form
        for (di,dj) in ((0,0),(0,1),(1,0),(1,1)):
            wi=(1-ti) if di==0 else ti; wj=(1-tj) if dj==0 else tj
            dwi=(-1 if di==0 else 1)/h; dwj=(-1 if dj==0 else 1)/h
            c=(i0+di)*nc+(j0+dj)
            # derivative wrt row (pointing -y) and col (+x)
            d_row=dwi*wj; d_col=wi*dwj
            rows.append(k); cols.append(c); vals_y.append(-d_row); vals_x.append(d_col)
    Dy=sp.csr_matrix((vals_y,(rows,cols)),shape=(N,nc*nc)); Dx=sp.csr_matrix((vals_x,(rows,cols)),shape=(N,nc*nc))
    G=-(sp.diags(gy[m])@Dy+sp.diags(gx[m])@Dx)      # r_unblurred = G dpsi
    # PSF blur as a sparse operator on masked pixels (small kernel)
    psf=np.array(fit.dataset.psf.kernel.native); ks=psf.shape[0]//2
    pos={(i,j):k for k,(i,j) in enumerate(idx)}; br=[];bc=[];bv=[]
    thr=psf.max()*1e-3
    for k,(i,j) in enumerate(idx):
        for a_ in range(-ks,ks+1):
            for b_ in range(-ks,ks+1):
                w=psf[ks+a_,ks+b_]
                if w<thr: continue
                q=pos.get((i+a_,j+b_))
                if q is not None: br.append(k); bc.append(q); bv.append(w)
    B=sp.csr_matrix((bv,(br,bc)),shape=(N,N))
    M=B@G; W=sp.diags(1/np.array(noise)[m])
    # curvature (Laplacian) regularisation on the coarse grid
    I=sp.identity(nc); D1=sp.diags([1.,-2.,1.],[-1,0,1],shape=(nc,nc))
    H=sp.kron(I,D1)+sp.kron(D1,I)
    MW=W@M; rW=W@resid[m]
    AtA=(MW.T@MW).tocsr(); HtH=(H.T@H).tocsr()
    if lam is None: lam=1e-2*AtA.diagonal().mean()/max(HtH.diagonal().mean(),1e-12)
    ridge=1e-6*AtA.diagonal().mean()
    sol=spl.spsolve((AtA+lam*HtH+ridge*sp.identity(AtA.shape[0])).tocsc(),MW.T@rW)
    psi=sol.reshape(nc,nc)
    # dkappa = 0.5 laplacian on the coarse grid
    h=step*ps
    lap=(np.roll(psi,1,0)+np.roll(psi,-1,0)+np.roll(psi,1,1)+np.roll(psi,-1,1)-4*psi)/h**2
    lap[0,:]=lap[-1,:]=lap[:,0]=lap[:,-1]=0
    dk=0.5*lap
    # only near data: coarse cells with masked pixels nearby
    cover=np.zeros((nc,nc),bool)
    for (i,j) in idx: cover[min(int(round(i/step)),nc-1),min(int(round(j/step)),nc-1)]=True
    cover=maximum_filter(cover,size=3)
    dks=gaussian_filter(dk,1.0)*cover
    peaks=(dks==maximum_filter(dks,size=5))&(dks>0)
    cand=[]
    for (a_,b_) in np.argwhere(peaks):
        y=( (n-1)/2 - a_*step )*ps; x=( b_*step - (n-1)/2 )*ps
        cand.append((float(dks[a_,b_]),float(y),float(x)))
    cand.sort(reverse=True); out=[]
    for c in cand:
        if all(np.hypot(c[1]-o[1],c[2]-o[2])>min_sep for o in out): out.append(c)
        if len(out)>=kmax: break
    return out,dks
