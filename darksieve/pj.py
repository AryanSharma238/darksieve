import numpy as np, autolens as al
from astropy.cosmology import Planck15 as cosmo
import astropy.units as u, astropy.constants as c
def sigma_crit_arcsec2(zl,zs):
    Dl=cosmo.angular_diameter_distance(zl); Ds=cosmo.angular_diameter_distance(zs); Dls=cosmo.angular_diameter_distance_z1z2(zl,zs)
    sc=(c.c**2/(4*np.pi*c.G)*Ds/(Dl*Dls)).to(u.Msun/u.Mpc**2)          # Msun / Mpc^2
    pc_per_arcsec=(Dl*u.arcsec.to(u.rad)).to(u.Mpc)                      # Mpc per arcsec
    return float((sc*pc_per_arcsec**2).value)                              # Msun / arcsec^2
def pj_pars(M,zl,zs,alpha,ra=0.01):
    """Vegetti+10 family: M=pi*b0*rs*Sigma_c ; rt=sqrt(alpha_sub*alpha) -> rs^3 = M*alpha/(pi*Sigma_c). returns (ra, rs, b0) in arcsec"""
    Sc=sigma_crit_arcsec2(zl,zs); rs=(M*alpha/(np.pi*Sc))**(1/3); b0=rs**2/alpha
    return ra,rs,b0,Sc
def pj_galaxy(y,x,M,zl,zs,alpha=1.4):
    ra,rs,b0,_=pj_pars(M,zl,zs,alpha)
    return al.Galaxy(redshift=zl,mass=al.mp.dPIEMassB0Sph(centre=(y,x),ra=ra,rs=rs,b0=b0))
if __name__=='__main__':
    zl,zs=0.222,0.609; Sc=sigma_crit_arcsec2(zl,zs); print('Sigma_crit = %.3e Msun/arcsec^2 = %.0f Msun/pc^2'%(Sc,Sc/(cosmo.angular_diameter_distance(zl).to('pc').value*4.848e-6)**2))
    M=3.5e9; ra,rs,b0,_=pj_pars(M,zl,zs,1.43); print('M=3.5e9: ra %.3f rs %.3f arcsec (rt=%.2f kpc) b0 %.4f'%(ra,rs,rs*cosmo.kpc_proper_per_arcmin(zl).value/60,b0))
    # numeric mass check
    g=al.Grid2D.uniform(shape_native=(1201,1201),pixel_scales=0.005); m=al.mp.dPIEMassB0Sph(centre=(0,0),ra=ra,rs=rs,b0=b0)
    kap=np.array(m.convergence_2d_from(grid=g)); print('numerical M_tot within 3": %.3e Msun (target %.2e)'%(kap.sum()*0.005**2*Sc,M))
