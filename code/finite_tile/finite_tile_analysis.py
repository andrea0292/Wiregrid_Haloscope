"""
Helpers for the finite-tile analysis of the 13-grid, 1.131 meV wire-grid stack (requires finite_tile.py).

Conventions
* Wires along y (parallel to B0), shorted to conducting end walls.  E = E_y(x, z) is normal to those walls,
  so the 2D solution is exact *inside the parallel-plate region between the end walls*.  How that
  y-uniform aperture field couples to a Gaussian receiver after the walls end is NOT solved here: it is
  approximated by the paraxial top-hat/Gaussian overlap 0.890 (a 3D aperture calculation is needed).
* Mirror at z = 0: PEC strip of the tile width (boundary elements); grids at z = j d.
* Z_s = R_s (1 - i), R_s = 1.8 x Posen law 2.9 mOhm (nu/3.9 GHz)^(2/3).
* Boost = forward-hemisphere power / (E0^2 W / 2 Z0)  (bare mirror of the same width).
"""
import numpy as np
import scipy.special as sp
from scipy.optimize import minimize_scalar, brentq
from finite_tile import FiniteTile, c0, mu0, Z0

a, r, Nl = 296.116e-6, 8.8835e-6, 13
dR = 0.9814226085352032*a                    # spacing retuned to 1.13137 meV (periodic M=0 multipole solver)
eV, h = 1.602176634e-19, 6.62607015e-34
nu_of_m = lambda m_eV: m_eV*eV/h
nuT = nu_of_m(1.13137e-3); lam = c0/nuT
Rsf = lambda nu, sc=1.8: 2.9e-3*sc*(nu/3.9e9)**(2/3)
Zs = lambda nu: Rsf(nu)*(1 - 1j)
sig = 220e3/np.sqrt(2); vE = 232e3
tg, wg = np.polynomial.hermite.hermgauss(24)

def make_tile(dv, W=0.10, dz=None, **kw):
    Nw = int(round(W/a)); Wm = kw.pop('Wm', Nw*a)
    t = FiniteTile(a, r, dv, Nl, Nw, Wm=Wm, lam_ref=lam, **kw)
    if dz is not None:                      # per-wire z displacements (grid tilts)
        t.wz = t.wz + dz; t.px = np.r_[t.wx, t.sx]; t.pz = np.r_[t.wz, t.sz]
    return t

def fwd_power(t, q, nth=8192):
    th = np.linspace(-np.pi/2, np.pi/2, nth, endpoint=False) + np.pi/(2*nth)
    S = t.farfield(q, th); return (1/(2*Z0))*(2/(np.pi*t.k))*np.sum(np.abs(S)**2, axis=0)*(np.pi/nth)

def tile_metric(tile, f, kax=0.0):
    tile.factor(f, Zs(f)); q = tile.solve_axion(kax)[:, 0]; return fwd_power(tile, q)/(tile.W/(2*Z0)), q

def tile_peak(tile, f_lo, f_hi, n=13):
    fs = np.linspace(f_lo, f_hi, n); v = np.array([tile_metric(tile, f)[0] for f in fs]); i = np.argmax(v)
    res = minimize_scalar(lambda f: -tile_metric(tile, f)[0], bracket=(fs[max(i-1, 0)], fs[i], fs[min(i+1, n-1)]), tol=1e-7)
    f0, pk = res.x, -res.fun
    hh = lambda f: tile_metric(tile, f)[0] - pk/2
    step = fs[1] - fs[0]; xl = f0 - step
    while hh(xl) > 0: xl -= step
    xr = f0 + step
    while hh(xr) > 0: xr += step
    return f0, pk, f0/(brentq(hh, f0, xr, xtol=1e-6*f0) - brentq(hh, xl, f0, xtol=1e-6*f0))

def line_matrix(tile, xl, zl):
    D = np.hypot(xl[:, None] - tile.px[None, :], zl - tile.pz[None, :]); H = sp.hankel1(0, tile.k*D)
    H[:, tile.nw:] *= tile.slen[None, :]; return H

def gauss_u(x, w): return (2/(np.pi*w**2))**0.25*np.exp(-x**2/w**2)

def mode_boost(El, xl, W, ws):
    dx = xl[1] - xl[0]; return np.array([abs(np.sum(El*gauss_u(xl, w))*dx)**2/W for w in ws])

def coh_1d(L, k, wind=0.0):
    wy = L/(2*0.99); yy = np.linspace(-L/2, L/2, 4001); uy = gauss_u(yy, wy)
    ky = k*(wind + np.sqrt(2)*sig*tg)/c0
    a2 = np.array([abs(np.trapezoid(np.exp(1j*kk*yy)*uy, yy))**2 for kk in ky]); a0 = abs(np.trapezoid(uy, yy))**2
    return np.sum(wg*a2)/np.sqrt(np.pi)/a0

def collect_chain(tile, f0, wind_x=0.0):
    tile.factor(f0, Zs(f0)); k = 2*np.pi*f0/c0; W = tile.W
    z1 = tile.wz.max() + 2e-3; xl = np.arange(-W/2 - 0.03, W/2 + 0.03, (c0/f0)/5); HL = line_matrix(tile, xl, z1); dx = xl[1] - xl[0]
    q = tile.solve_axion(0.0)[:, 0]; El = HL @ q; flux = np.sum(abs(El)**2)*dx/W
    ws = np.linspace(0.3, 0.7, 41)*W; mb = mode_boost(El, xl, W, ws); i = np.argmax(mb); w = ws[i]
    Em = HL @ tile.solve_axion(k*(wind_x + np.sqrt(2)*sig*tg)/c0)
    ax_avg = np.sum(wg*abs((gauss_u(xl, w) @ Em)*dx)**2/W)/np.sqrt(np.pi)
    ay = coh_1d(W, k)*0.8901
    return dict(flux=flux, beta2_mode_x=mb[i], w_over_W=w/W, beta2_mode_x_halo=ax_avg, eta_y=ay, beta2_single_mode=ax_avg*ay)

def gauss_incident(w, z1, k, nk=1201):
    kx = np.linspace(-0.999*k, 0.999*k, nk); dk = kx[1] - kx[0]
    A = (w/(2*np.sqrt(np.pi)))*(2/(np.pi*w**2))**0.25*np.exp(-kx**2*w**2/4); kz = np.sqrt(k**2 - kx**2)
    m = A > 1e-12*A.max(); kx, A, kz = kx[m], A[m], kz[m]
    Pinc = (1/(2*Z0))*2*np.pi*np.sum(abs(A)**2*kz/k)*dk
    def E(x, z): return (np.exp(1j*(np.outer(x, kx) - np.outer(z - z1, kz)))*A[None, :]) @ np.ones_like(A)*dk
    return E, Pinc

def H_field(tile, q, x, z):
    x = np.asarray(x).ravel(); z = np.asarray(z).ravel()
    dx = x[:, None] - tile.px[None, :]; dz = z[:, None] - tile.pz[None, :]; D = np.maximum(np.hypot(dx, dz), 1e-12)
    dH = -tile.k*sp.hankel1(1, tile.k*D); wgt = np.r_[np.ones(tile.nw), tile.slen]
    Ex = (dH*dx/D*wgt[None, :]) @ q; Ez = (dH*dz/D*wgt[None, :]) @ q
    om = 2*np.pi*tile.nu
    return -Ez/(1j*om*mu0), Ex/(1j*om*mu0)
