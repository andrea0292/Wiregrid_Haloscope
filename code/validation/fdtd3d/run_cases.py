"""Finite-wire-length 3D FDTD runs (Supplemental check of wire ends / tile frames).

Usage (from this directory):
    python run_cases.py A      # 1-2 meV cartridge: periodic reference + bars, walls, open  (~35 min on 10 cores)
    python run_cases.py B      # 0.25-0.5 meV cartridge: periodic reference + bars, walls, open (~30 min)
    python run_cases.py smoke  # short run to check the installation (~1 min)
Raw field records go to raw/ (large, not version-controlled); reduce.py turns them into ../../../data/fdtd3d/*.npz.
"""
import os, sys, time, numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import fdtd3d as F
RAW = os.path.join(HERE, 'raw'); os.makedirs(RAW, exist_ok=True)

def Rs(nu):  # copper surface resistance used throughout the paper (83.6 mOhm at 250 GHz, anomalous-skin scaling)
    return 83.6e-3*(nu/250e9)**(2/3)

CART = {  # band: (a [m], N grids, centre mass [meV])  -- Table I designs
    '1-2': (237.7e-6, 11, 1.41421),
    '0.25-0.5': (945.6e-6, 15, 0.353553),
}

def case(band, frames, W=0.10, Nx=6, ymargin=60, nsteps=60000, nf=451, span=(0.985, 1.03), tag=None, yper=False):
    a, N, m = CART[band]; D = a/Nx; nu_c = m*241.799e9
    r_e = 0.18518*D                      # intrinsic FDTD wire radius (calibrated vs sheet model; r_e/a = 0.031)
    Rp = Rs(nu_c)/(2*np.pi*r_e)          # resistance per unit length
    p = dict(D=D, Nx=Nx, N=N, kg=Nx, npml=30, zair=60, yper=yper, Rp=Rp, frames=frames,
             Wc=int(round(W/D)), ymargin=ymargin)
    g = F.build(p); kobs = g['kN'] + 30
    freqs = np.linspace(span[0], span[1], nf)*nu_c
    t = time.time()
    o = F.run(g, nu_c, 0.04, freqs, nsteps, kobs)
    tag = tag or f"{band}_{frames}_W{int(W*1000)}mm" + ("_per" if yper else "")
    np.savez_compressed(os.path.join(RAW, f'res_{tag}.npz'), freqs=freqs, accE=o['accE'], accH=o['accH'], E0=o['E0'],
                        trace=o['trace'], D=D, j0=g['j0'], j1=g['j1'], dt=g['dt'], nu_c=nu_c, a=a, N=N)
    print(tag, 'done', round(time.time()-t), 's', 'ring-down tail', abs(o['trace'][-3000:]).max()/abs(o['trace']).max(), flush=True)

if __name__ == '__main__':
    which = sys.argv[1]
    if which == 'A':
        case('1-2', 'bars', yper=True, tag='1-2_per')
        for fr in ['bars', 'walls', 'open']:
            case('1-2', fr)
    elif which == 'B':
        case('0.25-0.5', 'bars', yper=True, nsteps=150000, tag='0.25-0.5_per')
        for fr in ['bars', 'walls', 'open']:
            case('0.25-0.5', fr, nsteps=150000)
    elif which == 'smoke':
        case('1-2', 'bars', yper=True, nsteps=20000, nf=91, tag='smoke_per')
