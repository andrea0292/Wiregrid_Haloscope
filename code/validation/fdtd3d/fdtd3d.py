"""
3D Yee FDTD for an axion-driven, mirror-backed stack of wire grids.
  x : periodic, Nx cells (one wire per grid per period, wire at i = 0)  -> infinite across the wires
  y : along the wires (and B0). Periodic (infinite wires, reference) or finite: wires occupy j0 <= y/D <= j1,
      terminated by frames ('bars', 'walls' or 'open'), then free space, then an absorbing layer.
  z : PEC mirror at k = 0, grids at k = kg*(l+1), free space above, absorbing layer on top.
Wires: lines of E_y nodes with conductivity 1/(R' D^2), R' = R_s/(2 pi r); inductance set by the intrinsic
       FDTD wire radius (calibrated against the sheet model with the infinite-wire run).
Drive: uniform J_y(t) in 0 < z < L (and over the wire span for finite wires). Discrete particular field
       E_p(w) = dt J(w) / (2 i eps0 sin(w dt/2)), E0 = -E_p.
"""
import numpy as np
import numba as nb

c0 = 299792458.0; mu0 = 4e-7*np.pi; eps0 = 1.0/(mu0*c0**2); eta0 = mu0*c0

@nb.njit(parallel=True, cache=True)
def update_H(Ex, Ey, Ez, Hx, Hy, Hz, daHx, dbHx, daHy, dbHy, daHz, dbHz, yper):
    Nx, Ny, Nz = Ex.shape
    for j in nb.prange(Ny):
        if yper:
            jp = (j+1) % Ny; okj = True
        else:
            jp = j+1; okj = jp < Ny
        for i in range(Nx):
            ip = (i+1) % Nx
            for k in range(Nz):
                kp = k+1; okk = kp < Nz
                ezjp = Ez[i, jp, k] if okj else 0.0
                eykp = Ey[i, j, kp] if okk else 0.0
                Hx[i, j, k] = daHx[i, j, k]*Hx[i, j, k] - dbHx[i, j, k]*((ezjp - Ez[i, j, k]) - (eykp - Ey[i, j, k]))
                exkp = Ex[i, j, kp] if okk else 0.0
                Hy[i, j, k] = daHy[i, j, k]*Hy[i, j, k] - dbHy[i, j, k]*((exkp - Ex[i, j, k]) - (Ez[ip, j, k] - Ez[i, j, k]))
                exjp = Ex[i, jp, k] if okj else 0.0
                Hz[i, j, k] = daHz[i, j, k]*Hz[i, j, k] - dbHz[i, j, k]*((Ey[ip, j, k] - Ey[i, j, k]) - (exjp - Ex[i, j, k]))

@nb.njit(parallel=True, cache=True)
def update_E(Ex, Ey, Ez, Hx, Hy, Hz, caEx, cbEx, caEy, cbEy, caEz, cbEz, src, Jt, yper, inv_d):
    Nx, Ny, Nz = Ex.shape
    for j in nb.prange(Ny):
        if yper:
            jm = (j-1) % Ny; okj = True
        else:
            jm = j-1; okj = jm >= 0
        for i in range(Nx):
            im = (i-1) % Nx
            for k in range(Nz):
                km = k-1; okk = km >= 0
                hzjm = Hz[i, jm, k] if okj else 0.0
                hykm = Hy[i, j, km] if okk else 0.0
                Ex[i, j, k] = caEx[i, j, k]*Ex[i, j, k] + cbEx[i, j, k]*((Hz[i, j, k] - hzjm) - (Hy[i, j, k] - hykm))*inv_d
                hxkm = Hx[i, j, km] if okk else 0.0
                Ey[i, j, k] = caEy[i, j, k]*Ey[i, j, k] + cbEy[i, j, k]*(((Hx[i, j, k] - hxkm) - (Hz[i, j, k] - Hz[im, j, k]))*inv_d - Jt*src[i, j, k])
                hxjm = Hx[i, jm, k] if okj else 0.0
                Ez[i, j, k] = caEz[i, j, k]*Ez[i, j, k] + cbEz[i, j, k]*((Hy[i, j, k] - Hy[im, j, k]) - (Hx[i, j, k] - hxjm))*inv_d

@nb.njit(parallel=True, cache=True)
def dft_accum(Ey, Hx, kobs, ph_e, ph_h, accE, accH):
    Nx, Ny, Nz = Ey.shape
    nf = ph_e.shape[0]
    for j in nb.prange(Ny):
        se = 0.0; sh = 0.0
        for i in range(Nx):
            se += Ey[i, j, kobs]; sh += 0.5*(Hx[i, j, kobs-1] + Hx[i, j, kobs])
        se /= Nx; sh /= Nx
        for f in range(nf):
            accE[f, j] += ph_e[f]*se; accH[f, j] += ph_h[f]*sh

def grade(n, npml, order=3):
    return (np.clip(n, 0, npml)/npml)**order

def build(p):
    D = p['D']; Nx = p['Nx']; N = p['N']; kg = p['kg']; npml = p['npml']
    kN = kg*N; Nz = kN + p['zair'] + npml; yper = p['yper']
    if yper:
        Ny = 1; j0 = 0; j1 = 1
    else:
        Wc = p['Wc']; Ny = Wc + 2*(p['ymargin'] + npml); j0 = p['ymargin'] + npml; j1 = j0 + Wc
    dt = 0.99*D/(c0*np.sqrt(3.0)); shape = (Nx, Ny, Nz)
    sigmax = p.get('pml_scale', 1.0)*0.8*4/(eta0*D)
    kk = np.arange(Nz); kz0 = Nz - npml
    pz_i = grade(kk - kz0, npml); pz_h = grade(kk + 0.5 - kz0, npml)
    jj = np.arange(Ny)
    if not yper:
        py_i = grade(npml - jj, npml) + grade(jj - (Ny-1-npml), npml)
        py_h = grade(npml - (jj+0.5), npml) + grade(jj + 0.5 - (Ny-1-npml), npml)
    else:
        py_i = np.zeros(Ny); py_h = np.zeros(Ny)
    one = np.ones(shape)
    sEx = sigmax*(py_i[None, :, None] + pz_i[None, None, :])*one
    sEy = sigmax*(py_h[None, :, None] + pz_i[None, None, :])*one
    sEz = sigmax*(py_i[None, :, None] + pz_h[None, None, :])*one
    sHx = sigmax*(py_h[None, :, None] + pz_h[None, None, :])*one
    sHy = sigmax*(py_i[None, :, None] + pz_h[None, None, :])*one
    sHz = sigmax*(py_h[None, :, None] + pz_i[None, None, :])*one
    kw = [kg*(l+1) for l in range(N)]
    wire = np.zeros(shape, bool)
    for kwl in kw:
        wire[0, j0:j1, kwl] = True
    if p['Rp'] > 0:
        sEy[wire] += 1.0/(p['Rp']*D*D)
    def coefE(s):
        a = s*dt/(2*eps0); return (1-a)/(1+a), (dt/eps0)/(1+a)
    def coefH(s):
        a = (s*mu0/eps0)*dt/(2*mu0); return (1-a)/(1+a), (dt/(mu0*D))/(1+a)
    caEx, cbEx = coefE(sEx); caEy, cbEy = coefE(sEy); caEz, cbEz = coefE(sEz)
    daHx, dbHx = coefH(sHx); daHy, dbHy = coefH(sHy); daHz, dbHz = coefH(sHz)
    pec = {c: np.zeros(shape, bool) for c in 'xyz'}
    if p['Rp'] <= 0:
        pec['y'] |= wire
    pec['x'][:, :, 0] = True; pec['y'][:, :, 0] = True
    pec['x'][:, :, -1] = True; pec['y'][:, :, -1] = True
    if not yper:
        for jb in (0, Ny-1):
            pec['x'][:, jb, :] = True; pec['z'][:, jb, :] = True
        fr = p.get('frames', 'bars')
        if fr == 'bars':
            for kwl in kw:
                pec['x'][:, j0, kwl] = True; pec['x'][:, j1, kwl] = True
        elif fr == 'walls':
            for jb in (j0, j1):
                pec['x'][:, jb, :kw[-1]+1] = True; pec['z'][:, jb, :kw[-1]] = True
        elif fr != 'open':
            raise ValueError(fr)
    for c, (ca, cb) in (('x', (caEx, cbEx)), ('y', (caEy, cbEy)), ('z', (caEz, cbEz))):
        ca[pec[c]] = 0.0; cb[pec[c]] = 0.0
    src = np.zeros(shape); src[:, :, 1:kN] = 1.0; src[:, :, kN] = 0.5
    if not yper and p.get('src_y', 'tile') == 'tile':
        m = np.zeros(Ny); m[j0:j1] = 1.0; src *= m[None, :, None]
    return dict(shape=shape, dt=dt, D=D, kN=kN, kw=kw, j0=j0, j1=j1, Ny=Ny, Nz=Nz,
                caEx=caEx, cbEx=cbEx, caEy=caEy, cbEy=cbEy, caEz=caEz, cbEz=cbEz,
                daHx=daHx, dbHx=dbHx, daHy=daHy, dbHy=dbHy, daHz=daHz, dbHz=dbHz, src=src, yper=yper)

def run(g, f0, rel_bw, freqs, nsteps, kobs, J0=1.0, progress=0):
    shape = g['shape']; dt = g['dt']
    Ex = np.zeros(shape); Ey = np.zeros(shape); Ez = np.zeros(shape)
    Hx = np.zeros(shape); Hy = np.zeros(shape); Hz = np.zeros(shape)
    w = 2*np.pi*np.asarray(freqs); tau = 1.0/(2*np.pi*f0*rel_bw); t0 = 5*tau
    Ny = shape[1]
    accE = np.zeros((len(freqs), Ny), complex); accH = np.zeros((len(freqs), Ny), complex)
    Jw = np.zeros(len(freqs), complex); inv_d = 1.0/g['D']; trace = np.zeros(nsteps)
    for n in range(nsteps):
        update_H(Ex, Ey, Ez, Hx, Hy, Hz, g['daHx'], g['dbHx'], g['daHy'], g['dbHy'], g['daHz'], g['dbHz'], g['yper'])
        th = (n+0.5)*dt
        Jt = J0*np.exp(-((th-t0)/tau)**2)*np.sin(2*np.pi*f0*th)
        Jw += Jt*np.exp(1j*w*th)
        update_E(Ex, Ey, Ez, Hx, Hy, Hz, g['caEx'], g['cbEx'], g['caEy'], g['cbEy'], g['caEz'], g['cbEz'],
                 g['src'], Jt, g['yper'], inv_d)
        te = (n+1)*dt
        dft_accum(Ey, Hx, kobs, np.exp(1j*w*te), np.exp(1j*w*th), accE, accH)
        trace[n] = Ey[0, Ny//2, kobs]
        if progress and n % progress == 0:
            print(n, float(np.abs(Ey[:, :, kobs]).max()), flush=True)
    E0 = -(dt*Jw/(2j*eps0*np.sin(w*dt/2)))
    return dict(freqs=np.asarray(freqs), accE=accE, accH=accH, E0=E0, trace=trace, dt=dt)
