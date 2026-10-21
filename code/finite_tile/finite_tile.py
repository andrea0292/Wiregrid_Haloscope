"""
Finite-width (non-periodic) two-dimensional solver for an axion-driven wire-grid stack
in front of a finite conducting mirror.

Geometry (x across the wires, z normal to the stack, wires along y || B0):
  - grids l = 1..Nl at z_l = l*d, each with Nw wires at x_j = (j - (Nw-1)/2)*a, radius r;
  - a perfectly conducting mirror strip of width Wm at z = 0 (boundary-element, pulse basis);
  - optional perfectly conducting side walls (vertical strips) at x = +-Xw.
Wires: monopole (thin-wire) model with Leontovich condition E = zeta dE/dn on the surface,
zeta = i Z_s/(k Z0), Z_s = R_s - i X_s.  (The periodic multipole solver shows the dipole
correction changes beta^2 by 0.35% and the frequency by 0.007%.)

Fields: E = E_y(x, z), e^{-i omega t}.  Axion drive with momentum k_a = (kx, ky=0, kz):
particular solution -E0 exp(i k_a.r) (k_a << k).  The scattered part E_s is radiated by line
sources on the conductors and satisfies, on every conductor, E_s - zeta dE_s/dn = E0 exp(i k_a.r).
For the passive problem the drive is an incident field instead.
"""
import numpy as np
import scipy.special as sp
import scipy.linalg as sla

c0 = 2.99792458e8; mu0 = 4e-7*np.pi; Z0 = mu0*c0; eps0 = 1/(mu0*c0**2)
GAMMA = 0.5772156649015329

class FiniteTile:
    def __init__(self, a, r, d, Nl, Nw, Wm=None, mirror_margin=0.0, dseg=None, side_walls=False,
                 wall_x=None, wall_top=None, lam_ref=None):
        self.a, self.r, self.d, self.Nl, self.Nw = a, r, d, Nl, Nw
        xs = (np.arange(Nw) - (Nw-1)/2)*a
        zs = d*np.arange(1, Nl+1)
        X, Zw = np.meshgrid(xs, zs)
        self.wx, self.wz = X.ravel(), Zw.ravel()
        self.nw = self.wx.size
        W = Nw*a
        self.W = W
        lam = lam_ref if lam_ref is not None else 4*a
        dseg = dseg or lam/12
        Wm = Wm if Wm is not None else W + 2*mirror_margin
        nseg = int(np.ceil(Wm/dseg)); self.dseg = Wm/nseg
        self.sx = (np.arange(nseg) + 0.5)*self.dseg - Wm/2
        self.sz = np.zeros(nseg)
        self.slen = np.full(nseg, self.dseg)
        self.sdir = np.zeros(nseg, int)            # 0 = horizontal (along x), 1 = vertical
        if side_walls:
            xw = wall_x if wall_x is not None else W/2 + a/2
            ztop = wall_top if wall_top is not None else zs[-1] + d/2
            nv = int(np.ceil(ztop/dseg)); dz = ztop/nv
            zc = (np.arange(nv) + 0.5)*dz
            for sgn in (-1, 1):
                self.sx = np.r_[self.sx, np.full(nv, sgn*xw)]
                self.sz = np.r_[self.sz, zc]
                self.slen = np.r_[self.slen, np.full(nv, dz)]
                self.sdir = np.r_[self.sdir, np.ones(nv, int)]
        self.ns = self.sx.size
        self.px = np.r_[self.wx, self.sx]; self.pz = np.r_[self.wz, self.sz]
        self.N = self.nw + self.ns

    # ------------------------------------------------------------------
    def _matrix(self, k, zeta):
        nw, ns = self.nw, self.ns
        px, pz = self.px, self.pz
        D = np.hypot(px[:, None] - px[None, :], pz[:, None] - pz[None, :])
        np.fill_diagonal(D, 1.0)
        H = sp.hankel1(0, k*D)
        # segment sources: integrate over segment length (midpoint, with 2-pt Gauss for neighbours)
        L = self.slen
        G = H.copy()
        G[:, nw:] *= L[None, :]
        # refine near-neighbour segment interactions with 3-point Gauss-Legendre
        gx, gw = np.polynomial.legendre.leggauss(3)
        for jj in range(ns):
            j = nw + jj
            near = np.where(D[:, j] < 3*L[jj])[0]
            for i in near:
                if i == j: continue
                if self.sdir[jj] == 0:
                    xx = self.sx[jj] + 0.5*L[jj]*gx; zz = np.full(3, self.sz[jj])
                else:
                    xx = np.full(3, self.sx[jj]); zz = self.sz[jj] + 0.5*L[jj]*gx
                rr = np.hypot(px[i] - xx, pz[i] - zz)
                G[i, j] = 0.5*L[jj]*np.sum(gw*sp.hankel1(0, k*rr))
        # self terms
        r = self.r
        # wires: total field on own surface: c H0(kr) + a J0(kr); BC: c(H0+zeta k H1) + a(J0+zeta k J1) = rhs
        H0r = sp.hankel1(0, k*r); H1r = sp.hankel1(1, k*r)
        J0r = sp.j0(k*r); J1r = sp.j1(k*r)
        Jt = J0r + zeta*k*J1r
        A = G.copy()
        A[:nw, :] *= Jt                      # regular part from all other sources, times Jt
        idx = np.arange(nw)
        A[idx, idx] = H0r + zeta*k*H1r
        # segments: self term of pulse integral of H0 over own segment
        Ls = self.slen
        selfseg = Ls*(1 + (2j/np.pi)*(np.log(k*Ls/4) + GAMMA - 1))
        js = np.arange(nw, nw+ns)
        A[js, js] = selfseg
        return A

    def factor(self, nu, Zs):
        k = 2*np.pi*nu/c0
        zeta = 1j*Zs/(k*Z0)
        self.k, self.zeta, self.nu = k, zeta, nu
        A = self._matrix(k, zeta)
        self.lu = sla.lu_factor(A, overwrite_a=True, check_finite=False)
        return self

    def solve_axion(self, kax=0.0, kaz=0.0, E0=1.0):
        """Return source strengths for drive E0 exp(i (kax x + kaz z)) (columns if arrays)."""
        kax = np.atleast_1d(kax); kaz = np.broadcast_to(np.atleast_1d(kaz), kax.shape)
        rhs = E0*np.exp(1j*(np.outer(self.px, kax) + np.outer(self.pz, kaz)))
        return sla.lu_solve(self.lu, rhs, check_finite=False)

    def solve_incident(self, Einc):
        """Passive problem: Einc(x, z) callable giving the incident field; BC E_s - zeta dE_s/dn = -Einc."""
        rhs = -Einc(self.px, self.pz)
        return sla.lu_solve(self.lu, rhs, check_finite=False)

    # ------------------------------------------------------------------
    def field(self, q, x, z):
        """Scattered field E_s at points (x, z) from source strengths q (wires + segments)."""
        x = np.asarray(x); z = np.asarray(z)
        D = np.hypot(x.ravel()[:, None] - self.px[None, :], z.ravel()[:, None] - self.pz[None, :])
        D = np.maximum(D, self.r)
        H = sp.hankel1(0, self.k*D)
        H[:, self.nw:] *= self.slen[None, :]
        return (H @ q).reshape(x.shape + q.shape[1:]) if q.ndim > 1 else (H @ q).reshape(x.shape)

    def farfield(self, q, theta):
        """Far-field pattern S(theta) (theta from +z towards +x): E_s ~ sqrt(2/(pi k rho)) e^{i(k rho - pi/4)} S."""
        th = np.asarray(theta)
        ph = np.exp(-1j*self.k*(np.outer(np.sin(th), self.px) + np.outer(np.cos(th), self.pz)))
        # segment sinc factor depends on direction; apply exactly
        fac = np.ones((th.size, self.N), complex)
        for jj in range(self.ns):
            j = self.nw + jj
            proj = np.sin(th) if self.sdir[jj] == 0 else np.cos(th)
            fac[:, j] = self.slen[jj]*np.sinc(self.k*self.slen[jj]*proj/(2*np.pi))
        return (ph*fac) @ q

    def radiated_power(self, q, nth=8192):
        """Radiated power per unit length (W/m) for E0 = 1 V/m, integrated over all angles."""
        th = np.linspace(-np.pi, np.pi, nth, endpoint=False)
        S = self.farfield(q, th)
        return (1/(2*Z0))*(2/(np.pi*self.k))*np.sum(np.abs(S)**2, axis=0)*(2*np.pi/nth), th, S

    def wire_loss(self, q, Rs):
        """Ohmic loss per unit length on wires: (Rs/2) 2 pi r |H_phi|^2, H_phi ~ c k H1(kr)/(i omega mu0)."""
        c = q[:self.nw]
        Hphi = c*self.k*sp.hankel1(1, self.k*self.r)/(2*np.pi*self.nu*mu0)
        return np.sum(0.5*Rs*2*np.pi*self.r*np.abs(Hphi)**2, axis=0)
