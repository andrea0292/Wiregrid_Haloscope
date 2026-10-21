"""Penalty from random common-gap jitter during an integration.

The tabulated common-gap tolerance (SM, flatness table) is the gap error that shifts the resonance by the
detuning at which the Lorentzian power response drops by 10% (x = 1/3 in units of the half width).
A Gaussian jitter of rms sigma_d = m * tolerance smears the response; the scan figure of merit
I = int |beta|^4 dnu/nu drops by the ratio below and the minimum coupling rises as I^(-1/4).
"""
import numpy as np
x = np.linspace(-60, 60, 200001)          # detuning in units of the half width at half maximum
L = 1/(1 + x**2)
def ratio(s):
    k = np.exp(-0.5*(x/s)**2); k /= k.sum()
    Lm = np.convolve(L, k, mode='same')
    return (Lm**2).sum()/(L**2).sum()
if __name__ == '__main__':
    print('sigma_d/tolerance   I/I0    g penalty')
    for m in [0.5, 1, 2, 3, 5]:
        r = ratio(m/3); print(f'{m:>6}            {r:.3f}   {100*(r**-0.25-1):.1f}%')
