"""Stable sheet model: normalized axion drive, surface losses, and a PEC mirror.
Reused from the independently reconstructed model in the previous audit.
All lengths/frequencies SI; exp(-i omega t). beta is the outgoing field/E0.
"""
import numpy as np
from scipy.constants import c,mu_0,epsilon_0
Z0=np.sqrt(mu_0/epsilon_0)

def sheets(freq,a,r,N,d=None,factor=1.8,fixed_rs=None,gaps=None,Lfactors=None,external_drive=1.0):
    """Stable tridiagonal field solve; avoids evanescent transfer cancellation.
    Nodes are total E/E0 on each sheet. E=0 at mirror, outgoing at last sheet.
    A second right-hand side gives the reflection of unit incident radiation.
    external_drive is the axion-induced vacuum particular-field amplitude
    outside x=N*d relative to that inside: 1 for a uniform magnetic field,
    0 for a static field confined to the active depth. The inside drive is 1.
    """
    f=np.asarray(freq);shape=f.shape;f=f.ravel();w=2*np.pi*f;k=w/c
    if N==0:return external_drive*np.ones(shape,complex),np.zeros(shape)
    if d is None:d=a
    if gaps is None:gaps=np.repeat(d,N)
    if Lfactors is None:Lfactors=np.ones(N)
    Ls=mu_0*a/(2*np.pi)*np.log(a/(2*np.pi*r))
    rs=factor*2.9e-3*(f/3.9e9)**(2/3) if fixed_rs is None else fixed_rs*np.ones_like(f)
    rg=rs*a/(2*np.pi*r)
    phase=np.asarray(gaps)[:,None]*k[None,:]
    sn=np.sin(phase);cs=np.cos(phase);cot=cs/sn;isc=1/sn
    jump=-1j*Z0/(rg[None,:]-1j*w[None,:]*Ls*np.asarray(Lfactors)[:,None])
    diag=-cot-jump;rhs=((cs-1)/sn).astype(complex)
    diag[:-1]-=cot[1:];rhs[:-1]+=(cs[1:]-1)/sn[1:]
    diag[-1]+=1j;rhs[-1]-=1j*external_drive
    refl_rhs=np.zeros_like(diag);refl_rhs[-1]=2j
    # Thomas elimination. Only the last field is needed, so no back substitution.
    for j in range(1,N):
        fac=isc[j]/diag[j-1]
        diag[j]-=fac*isc[j]
        rhs[j]-=fac*rhs[j-1]
        refl_rhs[j]-=fac*refl_rhs[j-1]
    beta=rhs[-1]/diag[-1]+external_drive
    refl=refl_rhs[-1]/diag[-1]-1
    return beta.reshape(shape),(1-abs(refl)**2).reshape(shape)

