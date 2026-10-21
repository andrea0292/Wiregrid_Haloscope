"""Optimize the area-limited grid and calculate its finite-aperture spectrum."""
from pathlib import Path
import sys,json,time
import numpy as np
from scipy.constants import c,h,e
from scipy.optimize import minimize_scalar,brentq
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from design_origin import geometry_peak
from sheet_model import sheets
from aperture_model import Aperture

def sheet_design(m,N):
    g=geometry_peak(m,N);a=g['a_um']*1e-6;r=g['radius_um']*1e-6;f=m*1e-3*e/h
    def p(d):return float(abs(sheets(f,a,r,N,d=d,external_drive=0)[0])**2)
    d=minimize_scalar(lambda d:-p(d),bounds=(a*.98,a*1.02),method='bounded',options={'xatol':1e-16}).x
    zs=np.linspace(.85,1.15,12001);b=abs(sheets(zs*f,a,r,N,d=d,external_drive=0)[0])**2
    I=np.trapezoid(b*b,np.log(zs))
    return dict(N=N,a=a,r=r,d=d,L=N*d,I_1D=float(I),Q_response=g['Q_response'])

def run(m,lmax=23,nmax=4,quadrature_per_lambda=4.,cache=True):
    name=P/'results'/f'radiative_{m:g}.json'
    if name.exists() and cache:return json.loads(name.read_text())
    start=time.time();cs=[sheet_design(m,N) for N in range(5,25)]
    best=max(cs,key=lambda d:d['I_1D']);vol=max(cs,key=lambda d:d['I_1D']/d['L']**2)
    f=m*1e-3*e/h;lam=c/f;ny=max(100,int(np.ceil(.1/lam*quadrature_per_lambda)))
    o=Aperture(.1,.1,lmax=lmax,nmax=nmax,ny=ny)
    Y=o.radiation_admittance(f);a,r,N,d=[best[k] for k in ['a','r','N','d']]
    def response(z,w):return o.fast_response(z*f,a,r,N,d=d,Yrad=Y,port=o.gaussian_port(f,w))
    zpk=minimize_scalar(lambda z:-response(z,.05)[1],bounds=(.99,1.01),method='bounded',options={'xatol':1e-11}).x
    waist=minimize_scalar(lambda w:-response(zpk,w)[1],bounds=(.03,.07),method='bounded',options={'xatol':1e-8}).x
    port=o.gaussian_port(f,waist)
    q=best['Q_response'];zs=np.unique(np.r_[np.linspace(.85,zpk-10/q,101),np.linspace(zpk-10/q,zpk+10/q,801),np.linspace(zpk+10/q,1.15,101)])
    vals=np.array([o.fast_response(z*f,a,r,N,d=d,Yrad=Y,port=port)[:2] for z in zs])
    Ig=np.trapezoid(vals[:,1]**2,np.log(zs));Ir=np.trapezoid(vals[:,0]**2,np.log(zs))
    # Direct finite-device source-work and passive-receiver checks at the peak.
    o.nu_rad=zpk*f
    chk=o.solve(zpk*f,a,r,N,d=d,receiver_w=waist,frozen_radiation=True)
    out=dict(mass_meV=m,candidates=cs,best=best,volume_optimized=vol,
        finite=dict(I_gaussian=float(Ig),I_radiated=float(Ir),eta_I=float(np.sqrt(Ig/best['I_1D'])),
            peak_gaussian=float(vals[:,1].max()),peak_radiated=float(vals[:,0].max()),
            waist_m=float(waist),peak_over_nominal=float(zpk),modes=o.M,ny=ny,
            energy_balance_relative=chk['balance_relative'],gaussian_absorptivity=chk['gaussian_absorptivity'],
            radiation_min_eigenvalue=float(np.linalg.eigvalsh(Y.real)[0]),reciprocity_residual=o.reciprocity_residual),
        seconds=time.time()-start)
    np.savez(P/'results'/f'radiative_spectrum_{m:g}.npz',z=zs,gaussian=vals[:,1],radiated=vals[:,0])
    name.write_text(json.dumps(out,indent=2));print('RADIATIVE',m,best,out['finite'],flush=True)
    return out

if __name__=='__main__':
    for m in [.5,1.,2.]:run(m)
