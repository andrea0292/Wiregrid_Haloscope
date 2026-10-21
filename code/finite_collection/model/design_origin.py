"""Reproduce the reference design from a target band and stated constraints.
For each integer N, optimize pitch at fixed target frequency, with d=a and
r=max(0.03*a,2.5 micrometers). Follow the first collective resonance.
The pitch is then adjusted so the physical frequency-spectrum maximum is at
the target, including Rs(nu). Stored data expose the layer-count dependence.
"""
from pathlib import Path
import json
import numpy as np
from scipy.constants import c,h,e,mu_0,epsilon_0
from scipy.optimize import brentq,minimize_scalar
from sheet_model import sheets
ROOT=Path(__file__).resolve().parent
FACTOR=1.8;XI=.03;RMIN=2.5e-6
def resistance(f):return FACTOR*2.9e-3*(f/3.9e9)**(2/3)
def geometry_peak(m,N):
 f0=m*1e-3*e/h;lam=c/f0
 def power(x):
  a=x*lam;r=max(XI*a,RMIN)
  return float(abs(sheets(f0,a,r,N,factor=FACTOR)[0])**2)
 # Broad grid resolves all candidate peaks in this subwavelength design range.
 x=np.linspace(.16,.49,1601);p=np.array([power(t) for t in x])
 peaks=np.flatnonzero((p[1:-1]>p[:-2])&(p[1:-1]>p[2:]))+1
 if not len(peaks):raise RuntimeError(('no peak',m,N))
 # The lowest-frequency collective mode is the first maximum as pitch grows.
 j=peaks[0]
 opt=minimize_scalar(lambda t:-power(t),bounds=(x[j-1],x[j+1]),method='bounded',options={'xatol':1e-13})
 a=opt.x*lam
 def frequency_peak(a):
  r=max(XI*a,RMIN)
  fn=lambda z:float(abs(sheets(z*f0,a,r,N,factor=FACTOR)[0])**2)
  # Restrict to the first resonance, using the target as the known local seed.
  op=minimize_scalar(lambda z:-fn(z),bounds=(.995,1.005),method='bounded',options={'xatol':1e-13})
  return op.x,-op.fun
 for _ in range(5):
  rel,_=frequency_peak(a)
  if abs(rel-1)<2e-9:break
  # The radius floor breaks exact scale invariance; Newton slope by differences.
  ep=2e-5;slope=(frequency_peak(a*(1+ep))[0]-frequency_peak(a*(1-ep))[0])/(2*ep)
  a*=np.exp(-(rel-1)/slope)
 r=max(XI*a,RMIN);rel,pk=frequency_peak(a)
 fn=lambda z:float(abs(sheets(z*f0,a,r,N,factor=FACTOR)[0])**2)-pk/2
 # Half-width roots are local; expand just until the first crossing.
 step=2e-5;lo=rel-step;hi=rel+step
 while fn(lo)>0:lo-=step
 while fn(hi)>0:hi+=step
 low=brentq(fn,lo,rel,xtol=2e-14);high=brentq(fn,rel,hi,xtol=2e-14)
 ell=np.log(a/(2*np.pi*r));Ls=mu_0*a/(2*np.pi)*ell
 wp=1/np.sqrt(epsilon_0*a*Ls);G=resistance(f0)/(mu_0*r*ell)
 return {'mass_meV':m,'N':N,'a_um':a*1e6,'radius_um':r*1e6,'length_mm':N*a*1e3,
         'boost':pk,'Q_material':wp/G,'Q_response':rel/(high-low),'resonance_over_target':rel,
         'a_over_lambda':a/lam,'surface_resistance_ohm':resistance(f0)}

