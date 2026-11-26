import numpy as np
from scipy import integrate
from scipy.optimize import brentq
c=299792.458; H0=67.66; Om=0.3111; OL=1-Om; Oc=0.2621   # Planck18 (as used by colossus 'planck18'); Oc is their DM normalisation
G=6.674e-11; Ms=1.989e30; cm=299792458.
E=lambda z: np.sqrt(Om*(1+z)**3+OL)
_zg=np.linspace(0,6,6001); _Dc=np.concatenate([[0],integrate.cumulative_trapezoid(c/H0/E(_zg),_zg)])
Dc=lambda z: np.interp(z,_zg,_Dc)
def zm(dm,c2=200.):
    a=dm; b2=855.; d=(b2-a)**2-4*b2*(c2-a)
    if d<0: return np.nan
    z=((a-b2)+np.sqrt(d))/(2*b2); return z if z>0 else 1e-3
def dt_ms(y,M,zl): return 4*G*M*Ms/cm**3*(1+zl)*(y/2*np.sqrt(y*y+4)+np.log((np.sqrt(y*y+4)+y)/(np.sqrt(y*y+4)-y)))*1e3
def ymin(zl,M,w,ymax=6):
    f=lambda y: dt_ms(y,M,zl)-w
    return brentq(f,1e-12,ymax) if f(ymax)>0 else ymax
def ymax_R(R): return np.sqrt((1+R)/np.sqrt(R)-2)
def tau(zs,M,w,Rmax,eff=None,nz=60):
    if not np.isfinite(zs) or Rmax<1: return 0.
    zl=np.linspace(zs/nz/2,zs*(1-1/nz/2),nz); dz=zs/nz; tot=0
    ym=ymax_R(Rmax)
    for z in zl:
        Dl=Dc(z)/(1+z); Ds=Dc(zs)/(1+zs); Dls=(Dc(zs)-Dc(z))/(1+zs)
        y0=ymin(z,M,w)
        if ym<=y0: continue
        if eff is None: area=ym**2-y0**2
        else:
            ys=np.linspace(y0,ym,40); R=lambda y:(y*y+2+y*np.sqrt(y*y+4))/(y*y+2-y*np.sqrt(y*y+4))
            area=np.trapezoid([2*y*eff(dt_ms(y,M,z),R(y)) for y in ys],ys)
        tot+=1.5*Oc*H0**2/(c*H0*E(z))*Dls*Dl/Ds*area*(1+z)**2*dz*c/c  # H(z)=H0E; units: Mpc^-1 * Mpc
    return tot
def sample():
    B='/home/claude/pbh/data_zhou/Z/FRB_data/CHIME_cat2_frb/'
    a=np.concatenate([np.load(B+'chimefrbcat2_unique_first_repeaters.npy',allow_pickle=True),np.load(B+'chimefrbcat2_unique_non_repeater.npy',allow_pickle=True)])
    ok=np.isfinite(a['snr_fitb'])&(a['snr_fitb']>0)
    a=a[ok]; return a['dm_exc_ymw16'],a['width_fitb']*1e3,a['snr_fitb']
if __name__=='__main__':
    dm,w,snr=sample(); print('N sources',len(dm))
    for M in [100,300,467,609,1000,3000]:
        T=sum(tau(zm(d),M,wi,s/10) for d,wi,s in zip(dm,w,snr))
        print(f"M={M:5d}  Sum tau={T:7.2f}   1/Sum(tau)={1/T:.3f} (their 'f_up')   95%CL zero-detection 3.0/Sum={3/T:.3f}")
