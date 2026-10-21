"""Same-frequency, field, volume, and surface-resistance comparison.

The ALPHA-type comparator is an ideal port-coupled plasma resonator, NOT a
forecast for ALPHA Phase I at THz frequencies. It uses the square-lattice
bulk wire formula (Millar et al., 2210.00017, Eqs. 53-59), retuned at each
target frequency, and a specified mode overlap. No GHz Q or amplifier
performance is carried over to THz frequencies.
"""
from pathlib import Path
import json
import numpy as np
from scipy.constants import c, h, e, mu_0, epsilon_0
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar
from sheet_model import sheets

ROOT = Path(__file__).resolve().parent
F_SQUARE = np.pi/6 + sum(2/np.expm1(2*np.pi*n)/n for n in range(1,65))
C_RECT = 64/np.pi**4
C_SLAB = 8/np.pi**2
XI = .03
RMIN = 2.5e-6

def resistance(nu):
    return 1.8*2.9e-3*(nu/3.9e9)**(2/3)

def bulk_design(mass_meV, xi=XI):
    nu = mass_meV*1e-3*e/h
    lam = c/nu
    ell = np.log(1/(2*np.pi*xi)) + F_SQUARE
    pitch = lam/np.sqrt(2*np.pi*ell)
    radius = xi*pitch
    if radius < RMIN:
        radius = RMIN
        lower = 2*np.pi*radius*np.exp(-F_SQUARE)*(1+1e-9)
        pitch = brentq(lambda a: a*a*(np.log(a/(2*np.pi*radius))+F_SQUARE)
                       -lam*lam/(2*np.pi), lower, 2*lam, xtol=1e-18)
        ell = np.log(pitch/(2*np.pi*radius))+F_SQUARE
    qu = 2*np.pi*nu*mu_0*radius*ell/resistance(nu)
    return {'frequency_GHz':nu/1e9,'wavelength_mm':lam*1e3,
            'pitch_um':pitch*1e6,'radius_um':radius*1e6,'Q_unloaded':qu,
            'surface_resistance_ohm':resistance(nu),
            'amplitude_transport_length_mm':lam*np.sqrt(2*qu)/(2*np.pi)*1e3}

