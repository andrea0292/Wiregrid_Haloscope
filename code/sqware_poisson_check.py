"""Low-count check of the SQWARE configuration-3 dark-photon comparison (Supplemental S5).

Signal rate s(t)=s0/(1+(2t/tl)^2) (Lorentzian power response, FWHM traversed in time tl),
known background b = 1e-4 /s, with b*tl = 0.062 counts per response width (common scan protocol).
Gaussian threshold: int s^2/b dt = 25.  Asimov: 2 int[(s+b)ln(1+s/b)-s] dt = 25.
The signal rate s0 is proportional to chi^2, so chi scales as s0^(1/2); the Gaussian SQWARE/grid ratio 1.40 becomes
1.40*(s0_A/s0_G)^(1/2). (Earlier versions used a fourth root, which understated the low-count penalty.)
Not used for the comparisons quoted in the current Letter (End Matter B uses sqware_check.py).
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import quad
b = 1e-4; tl = 0.062 / b
ZG = lambda s0: s0**2 * tl * np.pi / 4 / b
def ZA(s0):
    f = lambda t: (lambda s: (s + b) * np.log1p(s / b) - s)(s0 / (1 + (2 * t / tl)**2))
    return 4 * quad(f, 0, 200 * tl, limit=400)[0]
s0G = brentq(lambda s: ZG(s) - 25, 1e-8, 10); s0A = brentq(lambda s: ZA(s) - 25, 1e-8, 10)
print("chi factor %.3f -> ratio %.2f" % ((s0A / s0G)**0.5, 1.40 * (s0A / s0G)**0.5))
