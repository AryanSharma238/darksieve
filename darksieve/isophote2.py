import numpy as np, sys
from scipy.ndimage import map_coordinates, gaussian_filter
from scipy.optimize import minimize
t=np.linspace(0,2*np.pi,240,endpoint=False)
def samp(img,A,q,pa,cy,cx):
    x=A*np.cos(t); y=A*q*np.sin(t)
    return map_coordinates(img,[cy+x*np.sin(pa)+y*np.cos(pa),cx+x*np.cos(pa)-y*np.sin(pa)],order=3)
def harm(I,ms):
    M=np.c_[[np.ones_like(t)]+[f(m*t) for m in ms for f in (np.cos,np.sin)]].T
    return np.linalg.lstsq(M,I,rcond=None)[0]
def e2qpa(e1,e2):
    e=min(np.hypot(e1,e2),0.8); return (1-e)/(1+e),0.5*np.arctan2(e2,e1)
def fit(img,A,p0):
    q0,pa0,cy0,cx0=p0; e=(1-q0)/(1+q0); x0=[e*np.cos(2*pa0)+0.02,e*np.sin(2*pa0)+0.02,cy0,cx0]
    def cost(p):
        q,pa=e2qpa(p[0],p[1]); c=harm(samp(img,A,q,pa,p[2],p[3]),[1,2]); return (c[1:]**2).sum()/c[0]**2
    r=minimize(cost,x0,method='Nelder-Mead',options={'xatol':1e-5,'fatol':1e-14,'maxiter':3000}); q,pa=e2qpa(r.x[0],r.x[1]); return np.array([q,pa,r.x[2],r.x[3]])
def run(name,radii,ps=0.05):
    D=np.load('/home/claude/lens/hst/%s_inputs.npz'%name); img=gaussian_filter(D['img'],0.6); p=[0.85,0.3,200.,200.]; out=[]
    for a in radii:
        A=a/ps; p=fit(img,A,p); q,pa,cy,cx=p
        I=samp(img,A,q,pa,cy,cx); c=harm(I,[1,2,3,4])
        dI=(samp(img,A*1.05,q,pa,cy,cx).mean()-samp(img,A*0.95,q,pa,cy,cx).mean())/(0.1*A)
        s=-dI*A
        a3,b3,a4,b4=c[5]/s,c[6]/s,c[7]/s,c[8]/s
        out.append(dict(a=a,q=q,pa=np.degrees(pa)%180,cy=cy,cx=cx,A3=np.hypot(a3,b3),A4=np.hypot(a4,b4),a4=a4))
        print('%s a=%.2f q=%.3f pa=%5.1f ctr=(%.2f,%.2f)  |m3|=%.4f |m4|=%.4f  a4(cos,+disky)=%+.4f'%(name,a,q,np.degrees(pa)%180,cy,cx,np.hypot(a3,b3),np.hypot(a4,b4),a4))
    return out
if __name__=='__main__':
    for n in sys.argv[1:]: run(n,[0.2,0.3,0.4,0.5,0.6,0.7])
