"""Independent FDFD cross-check of the multipole solver (lossless, 3 grids, r/a = 0.1).
Reproduces fdfd_crosscheck.csv (~15 min at dx = a/400)."""
import numpy as np, pandas as pd
from scipy.optimize import curve_fit
import os
HERE = os.path.dirname(os.path.abspath(__file__))
exec(open(os.path.join(HERE, '..', 'fullwave.py')).read()); exec(open(os.path.join(HERE, 'fdfd_check.py')).read())
a2 = 296e-6; r2 = 0.1*a2; zl3 = a2*np.arange(1, 4)
st = Stack(a2, r2, zl3, M=4)
ff = np.linspace(385e9, 415e9, 31)
lor = lambda f, A, f0, g, b: A/(1 + ((f - f0)/g)**2) + b
rows = []
curves = {'multipole M=4': np.array([abs(st.solve(f, 0.0, 'axion'))**2 for f in ff])}
for nd in [200, 400]:
    curves[f'FDFD dx=a/{nd}'] = np.array([fdfd_beta(f, a2, r2, zl3, a2/nd)[0]**2 for f in ff])
for lab, v in curves.items():
    p, _ = curve_fit(lor, ff/1e9, v, p0=[v.max(), ff[np.argmax(v)]/1e9, 2, 1])
    rows.append(dict(method=lab, f0_GHz=p[1], peak_beta2=p[0] + p[3], FWHM_GHz=2*abs(p[2])))
print(pd.DataFrame(rows)); pd.DataFrame(rows).to_csv('fdfd_crosscheck_repro.csv', index=False)
