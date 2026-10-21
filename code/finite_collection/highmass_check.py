"""Transverse-mode check at 4 meV with paraxial aperture admittance."""
from pathlib import Path
import sys,json,numpy as np
from scipy.constants import e,h
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from aperture_model import Aperture
r=json.loads((P/'results/radiative_4.json').read_text());g=r['best'];v=r['finite'];f=4e-3*e/h
rows=[]
for lm in [23,47,95]:
    o=Aperture(.1,.1,lmax=lm,nmax=4,ny=128);Y=np.eye(o.M,dtype=complex);port=o.gaussian_port(f,v['waist_m'])
    zs=np.linspace(v['peak_over_nominal']-4/g['Q_response'],v['peak_over_nominal']+4/g['Q_response'],121)
    b=np.array([o.fast_response(z*f,g['a'],g['r'],g['N'],d=g['d'],Yrad=Y,port=port)[1] for z in zs])
    row=dict(lmax=lm,modes=o.M,I_central=float(np.trapezoid(b*b,np.log(zs))),peak=float(b.max()),relative_peak_to_full60=float(b.max()/v['peak_gaussian']-1))
    rows.append(row);print(row,flush=True)
(P/'results/highmass_basis.json').write_text(json.dumps(rows,indent=2))
