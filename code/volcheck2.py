import sys,json,time
from pathlib import Path
import numpy as np
from scipy.constants import c,h,e
from scipy.optimize import minimize_scalar
P=Path('.').resolve(); sys.path.insert(0,str(P/'model')); sys.path.insert(0,str(P))
from aperture_model import Aperture
from distributed_port import PlasmaPort,grid_shapes
import importlib.util
spec=importlib.util.spec_from_file_location('ar',P/'assess_resources.py')
src=open(P/'assess_resources.py').read().split("out={'halo_corrected'")[0]
ns={'__file__':str(P/'assess_resources.py')}; exec(src,ns)
port_halo,rad_halo,guide_loss=ns['port_halo'],ns['rad_halo'],ns['guide_loss']
m=float(sys.argv[1]); r=json.load(open(f'results/radiative_{m:g}.json')); v=r['volume_optimized']
f=m*1e-3*e/h;lam=c/f;ny=max(100,int(np.ceil(.1/lam*4)))
o=Aperture(.1,.1,lmax=23,nmax=4,ny=ny);Y=o.radiation_admittance(f);a,rr,N,d=[v[k] for k in ['a','r','N','d']]
resp=lambda z,w:o.fast_response(z*f,a,rr,N,d=d,Yrad=Y,port=o.gaussian_port(f,w))
zpk=minimize_scalar(lambda z:-resp(z,.05)[1],bounds=(.99,1.01),method='bounded',options={'xatol':1e-11}).x
w=minimize_scalar(lambda w:-resp(zpk,w)[1],bounds=(.03,.07),method='bounded',options={'xatol':1e-8}).x
port=o.gaussian_port(f,w); q=v['Q_response']
zs=np.unique(np.r_[np.linspace(.85,zpk-10/q,101),np.linspace(zpk-10/q,zpk+10/q,801),np.linspace(zpk+10/q,1.15,101)])
vals=np.array([o.fast_response(z*f,a,rr,N,d=d,Yrad=Y,port=port)[:2] for z in zs])
Ig=float(np.trapezoid(vals[:,1]**2,np.log(zs)))
rv={'best':v,'finite':dict(peak_over_nominal=float(zpk),waist_m=float(w),ny=ny)}
chi=rad_halo(m,rv)
obj=PlasmaPort(m,N*d); rows=[]
for K in [1,4,16,64]:
    cands=[obj.optimize(nx,ny_)[0] for nx,ny_ in grid_shapes(K)]
    for row in cands:
        hp=port_halo(obj,row); row['R_halo']=Ig/row['I']*(chi/hp['chi_peak'])**2
    rows.append(min(cands,key=lambda x:x['R_halo']))
cont=obj.optimize(uniform=True)[0]; hp=port_halo(obj,cont); cont['R_halo']=Ig/cont['I']*(chi/hp['chi_peak'])**2; rows.append(cont)
Rbest=min(x['R_halo'] for x in rows); gg=guide_loss(m); thr=-10*np.log10(.8*np.sqrt(Rbest))
res=dict(mass=m,N=N,L_mm=N*d*1e3,I_G=Ig,eta_I=float(np.sqrt(Ig/v['I_1D'])),chi_rad=chi,R_halo=[x['R_halo'] for x in rows],R_best=Rbest,parity_dB=thr,guide_cm=100*thr/gg['power_loss_dB_m'])
print(json.dumps(res),flush=True); json.dump(res,open(f'volcheck2_{m:g}.json','w'))
