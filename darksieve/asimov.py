"""Asimov audit: noise-free image of the lens's smooth model (parametric source fitted to the real data), searched with the
pipeline's own pixelized-source model at the true lens parameters. Any arc position where a subhalo still raises the
evidence is a source-regularization hotspot that can produce false detections.
usage: asimov.py CFG out.json [spacing]"""
import sys, os, json, numpy as np
sys.path.insert(0,os.path.dirname(os.path.abspath(__file__)))
cfgp,outp=sys.argv[1],sys.argv[2]; sp=float(sys.argv[3]) if len(sys.argv)>3 else 0.3
sys.argv=['bench.py',cfgp,'0','0','/tmp/none.jsonl']
src=open(os.path.join(os.path.dirname(os.path.abspath(__file__)),'bench.py')).read(); g={'__file__':os.path.abspath(__file__)}
exec(src[:src.index('rng=np.random.default_rng(7000+wid)')],g)
L,Qt,data0,ds,al,srcG,S=(g[k] for k in ('L','Qt','data0','ds','al','srcG','S'))
t=ds.tracer(L,Qt); T0=np.array(al.FitImaging(dataset=data0,tracer=al.Tracer(galaxies=list(t.galaxies)[:-1]+[srcG])).model_data.native)
L.set_data(T0); v0=ds.post(L,Qt); pts=S.arc_positions(L,Qt,spacing=sp); res=[]
for i,p in enumerate(pts):
    gm=max((ds.post(L,Qt,[(p[0]+0.0013,p[1]+0.0017,m)])-v0,m) for m in (8.5,9.0,9.5,10.0))
    res.append([float(p[0]),float(p[1]),float(gm[0]),gm[1]])
    if i%20==0: print('pos',i,len(pts),flush=True)
res.sort(key=lambda r:-r[2]); n4=sum(r[2]>4 for r in res)
out=dict(name=g['cfg']['name'],n=len(res),max_gain=res[0][2],n_hot4=n4,top=res[:8],all=res); json.dump(out,open(outp,'w'),indent=1)
print('ASIMOV',json.dumps({k:out[k] for k in ('name','n','max_gain','n_hot4')}),'top',[(round(a,2),round(b,2),round(c,1)) for a,b,c,_ in res[:4]],flush=True)
