"""Lens-specific detection threshold from the emulator.
P_FP(T) = E_prior[ sigmoid(a + b*log1p(score(error)) - g*T) ], prior over lens-model errors:
  null 0.15 / m=1 0.30 / m=3,4 0.15 / source knots 0.10 (relative weights of the benchmark), amplitudes U(0.01,0.02), random angle.
T_lens = smallest T with P_FP(T) <= alpha.   usage: decide.py emu.npz calib.json [alpha]"""
import sys, json, numpy as np
def load(emup,calp):
    Z=np.load(emup); w=json.load(open(calp))['w']; G=Z['G']; NN=len(G)
    noise=np.random.default_rng(1).standard_normal((1000,NN))@((G-G.mean(0))/np.sqrt(NN-1)); return Z,w,noise
def scores(Z,noise,n=400,seed=2):
    rng=np.random.default_rng(seed); half=Z['half']; kinds=rng.choice(['null','m1','m34','knots'],size=n,p=np.array([.15,.3,.15,.1])/.7); out=[]
    for k in kinds:
        v=np.zeros_like(half)
        if k=='m1': a=rng.uniform(.01,.02); t=rng.uniform(0,2*np.pi); v=Z['E_m1c']*a*np.cos(t)/.01+Z['E_m1s']*a*np.sin(t)/.01
        if k=='m34': m=int(rng.choice([3,4])); a=rng.uniform(.01,.02); t=rng.uniform(0,2*np.pi); v=Z['E_m%dc'%m]*a*np.cos(t)/.01+Z['E_m%ds'%m]*a*np.sin(t)/.01
        out.append(np.log1p(max(float(np.median((noise[:300]+v-half).max(1))),0)))
    return np.array(out)
def pfp(S,w,T): return float(np.mean(1/(1+np.exp(-(w[0]+w[1]*S-w[2]*T)))))
def threshold(S,w,alpha=0.05):
    for T in np.arange(3,60,0.5):
        if pfp(S,w,T)<=alpha: return float(T)
    return 60.0
if __name__=='__main__':
    Z,w,noise=load(sys.argv[1],sys.argv[2]); a=float(sys.argv[3]) if len(sys.argv)>3 else 0.05; S=scores(Z,noise)
    T=threshold(S,w,a); print(json.dumps(dict(emu=sys.argv[1],alpha=a,T_lens=T,pfp_at_5=pfp(S,w,5),pfp_at_12=pfp(S,w,12))))
