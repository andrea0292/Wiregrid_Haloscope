from pathlib import Path
import sys,json,time
import numpy as np
from scipy.constants import c,h,e
from scipy.sparse import diags,kron,eye
from scipy.sparse.linalg import splu
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from distributed_port import PlasmaPort,single_mode_validation,spectrum_summary
from aperture_model import Aperture

def fd_check(obj,row,Nx,Ny):
    obj.pp2=row['plasma_ratio_squared']
    L,W=obj.L,obj.W;hx=L/(Nx+1);hy=W/(Ny+1)
    tx=diags([-np.ones(Nx-1),2*np.ones(Nx),-np.ones(Nx-1)],[-1,0,1])/hx**2
    ty=diags([-np.ones(Ny-1),2*np.ones(Ny),-np.ones(Ny-1)],[-1,0,1])/hy**2
    lap=kron(tx,eye(Ny))+kron(eye(Nx),ty)
    x=np.arange(1,Nx+1)*hx;y=np.arange(1,Ny+1)*hy
    gx=sum(np.exp(-((x-cc)/obj.w)**2) for cc in (np.arange(row['nx'])+.5)*L/row['nx'])
    gy=sum(np.exp(-((y-cc)/obj.w)**2) for cc in (np.arange(row['ny'])+.5)*W/row['ny'])
    u=np.outer(gx,gy).ravel();u/=np.linalg.norm(u)
    f=np.ones(Nx*Ny)*np.sqrt(hx*hy/(L*W))
    z=row['peak_over_nominal'];gamma=obj.Gamma/(2*np.pi*obj.nu)*z**(2/3)
    diag=obj.Q*(-z*z+obj.pp2*z/(z+1j*gamma))
    D=(obj.Q/obj.k**2*lap+diag*eye(Nx*Ny)).tocsc()
    sol=splu(D).solve(np.column_stack((z*z*f,u)).astype(complex))
    G=np.vdot(u,sol[:,1]);H=np.vdot(u,sol[:,0]);beta=z*(row['alpha']-1j*row['zeta'])
    t=H/(1+beta*G)
    power=obj.k*obj.L*obj.Q*row['zeta']*abs(t)**2
    u0=obj.profile(row['nx'],row['ny']);G0,H0=obj.greens(u0,np.array([z]))
    ref=obj.spectrum(G0,H0,row['alpha'],row['zeta'],np.array([z]))[0]
    return dict(Nx=Nx,Ny=Ny,power=float(power),modal_power=float(ref),relative_difference=float(power/ref-1))

out={'single_mode':single_mode_validation(),'port_refinement':[],'finite_difference':[],'aperture_refinement':[],'pickup_width_scan':[]}
L=json.loads((P/'results/radiative_1.json').read_text())['best']['L']
row=json.loads((P/'results/ports_1_same_active_volume.json').read_text())['rows'][1]
for scale,step in [(1.,.05),(1.5,.025),(2.,.025)]:
    obj=PlasmaPort(1.,L,basis_scale=scale,step=step)
    obj.pp2=row['plasma_ratio_squared']
    u=obj.profile(row['nx'],row['ny']);G,H=obj.greens(u)
    b=obj.spectrum(G,H,row['alpha'],row['zeta'])
    out['port_refinement'].append(dict(basis_scale=scale,step=step,modes=len(obj.f),I=float(np.trapezoid(b*b,np.log(obj.z))),**spectrum_summary(obj.z,b)))
obj=PlasmaPort(1.,L,basis_scale=2.)
for nx,ny in [(40,400),(80,800)]:
    res=fd_check(obj,row,nx,ny);out['finite_difference'].append(res);print('FD',res,flush=True)
for m in [.5,1.,2.]:
    r=json.loads((P/'results'/f'radiative_{m:g}.json').read_text());L=r['best']['L']
    for ww in [1.,2.,4.]:
        obj=PlasmaPort(m,L,pickup_width_pitches=ww)
        for K in [4,16]:
            rr,b,u=obj.optimize(1,K)
            rr['pickup_width_pitches']=ww;rr['R']=r['finite']['I_gaussian']/rr['I']
            out['pickup_width_scan'].append(rr)
            print('WIDTH',m,ww,K,rr['R'],flush=True)

r=json.loads((P/'results/radiative_1.json').read_text());g=r['best'];nu=1e-3*e/h
for lm,ny in [(23,323),(31,450)]:
    o=Aperture(.1,.1,lmax=lm,nmax=4,ny=ny);Y=o.radiation_admittance(nu);port=o.gaussian_port(nu,r['finite']['waist_m'])
    zpk=r['finite']['peak_over_nominal'];zs=np.linspace(zpk-8/g['Q_response'],zpk+8/g['Q_response'],401)
    bs=np.array([o.fast_response(z*nu,g['a'],g['r'],g['N'],d=g['d'],Yrad=Y,port=port)[1] for z in zs])
    I=float(np.trapezoid(bs*bs,np.log(zs)))
    out['aperture_refinement'].append(dict(modes=o.M,ny=ny,I_central=I,peak=float(bs.max())))
    print('APERTURE',out['aperture_refinement'][-1],flush=True)
# Add counting-filter integrals to the radiative result without changing it.
out['radiative_FWHM']={}
for m in [.5,1.,2.]:
    a=np.load(P/'results'/f'radiative_spectrum_{m:g}.npz')
    out['radiative_FWHM'][str(m)]=spectrum_summary(a['z'],a['gaussian'])
(P/'results/validation.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='pickup_width_scan'},indent=2),flush=True)
