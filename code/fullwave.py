"""
Two-dimensional full-wave (multipole / lattice-sum) solver for an axion-driven,
mirror-backed stack of periodic wire grids.

Geometry: wires parallel to y (and to B0), infinitely long; grids are periodic along x
with pitch a; layer l (l = 1..N) lies at height z_l above a perfectly conducting mirror
at z = 0.  Field E = E_y(x, z) e^{-i omega t}.  Units: lengths in metres, E0 = 1.

Axion drive.  In vacuum the total field obeys (lap + k^2) E = -k^2 E0, with particular
solution -E0.  Writing E = -E0 + E0 e^{ikz} + E_w, the remainder E_w solves the homogeneous
Helmholtz equation, vanishes on the mirror (handled by image arrays), is outgoing, and
satisfies on each wire surface the Leontovich condition E = zeta dE/dn with
zeta = i Z_s/(k Z0), Z_s = R_s - i X_s (e^{-i omega t}; X_s > 0 is inductive).

E_w = sum_{l,n} c_{l,n} sum_p H_n^(1)(k|r - r_{l,p}|) e^{i n theta}  - (images).
Same-layer coupling: exact 1D lattice sums (Kummer-accelerated with polylogarithms).
Inter-layer and image coupling: Floquet (plane-wave) expansions incl. evanescent orders.
The emitted zero-order amplitude normalised to a bare mirror is beta.
"""
import numpy as np
import scipy.special as sp
import mpmath

c0 = 2.99792458e8
mu0 = 4e-7*np.pi
Z0 = mu0*c0

# ---------------- lattice sums ----------------
def _asym_coeffs(l, K):
    # Hankel asymptotic: H_l(z) ~ sqrt(2/(pi z)) e^{i(z - l pi/2 - pi/4)} sum_k i^k a_k / z^k
    mu = 4*l*l; out = []; ak = 1.0
    for k in range(K+1):
        if k > 0:
            ak *= (mu - (2*k-1)**2)/(k*8.0)
        out.append(ak*(1j)**k)
    return out

def half_sum_H(l, x, P=3000, K=8):
    """sum_{p=1}^infty H_l^(1)(p x) for 0 < x < 2 pi (one-sided)."""
    co = _asym_coeffs(l, K)
    pref = np.sqrt(2/np.pi)*np.exp(-1j*(l*np.pi/2 + np.pi/4))
    p = np.arange(1, P+1, dtype=float)
    exact = sp.hankel1(l, p*x)
    asym = np.zeros_like(exact)
    for k, ck in enumerate(co):
        asym += pref*ck*np.exp(1j*p*x)*(p*x)**(-0.5-k)
    tail = 0j
    for k, ck in enumerate(co):
        tail += complex(pref*ck*x**(-0.5-k)*mpmath.polylog(0.5+k, mpmath.exp(1j*x)))
    return np.sum(exact-asym) + tail

def lattice_sums(x, L, P=3000, K=8):
    """S_l = sum_{p != 0} H_l(k|p|a) e^{i l phi_p}, phi_p = angle of vector (source -> origin);
    returns dict for l = -L..L (odd l vanish at normal incidence)."""
    S = {}
    for l in range(0, L+1):
        if l % 2:
            S[l] = 0j; S[-l] = 0j
        else:
            v = 2*half_sum_H(l, x, P, K)
            S[l] = v; S[-l] = v
    return S

def twersky_S0(k, a, mmax=20000):
    """Independent spectral (Twersky-type) evaluation of S_0 for ka < 2 pi."""
    g = 0.5772156649015329
    m = np.arange(1, mmax+1)
    kx = 2*np.pi*m/a
    kap = np.sqrt(kx**2 - k**2 + 0j)
    corr = (2/a)*2*np.sum((1/(1j*kap) - a/(2*np.pi*1j*m)))
    return 2/(k*a) - 1 - (2j/np.pi)*(g + np.log(k*a/(4*np.pi))) + corr

# ---------------- Floquet translation factors ----------------
def _fplus(n, kx, kz, k):
    return ((kz - 1j*kx)/k)**n if n >= 0 else ((kz + 1j*kx)/k)**(-n)
def _fminus(n, kx, kz, k):
    return (-(kz + 1j*kx)/k)**n if n >= 0 else ((1j*kx - kz)/k)**(-n)

class Stack:
    def __init__(self, a, r, z_layers, M=3, Mfl=15):
        self.a, self.r = a, r
        self.z = np.asarray(z_layers, float)
        self.N = len(self.z); self.M = M; self.Mfl = Mfl
        self.q = np.arange(-M, M+1)

    def _setup(self, nu, Zs):
        a, r, M = self.a, self.r, self.M
        k = 2*np.pi*nu/c0
        zeta = 1j*Zs/(k*Z0)
        q = self.q
        H = sp.hankel1(q, k*r); Hp = sp.h1vp(q, k*r)
        J = sp.jv(q, k*r);     Jp = sp.jvp(q, k*r)
        Ht = H - zeta*k*Hp; Jt = J - zeta*k*Jp
        S = lattice_sums(k*a, 2*M)
        m = np.arange(-self.Mfl, self.Mfl+1)
        kx = 2*np.pi*m/a
        kz = np.sqrt(k**2 - kx**2 + 0j); kz = np.where(np.imag(kz) < 0, -kz, kz)
        return k, Ht, Jt, S, kx, kz

    def _matrix(self, k, Ht, Jt, S, kx, kz):
        """Linear map c -> regular coefficients a (at every wire centre)."""
        N, M, a = self.N, self.M, self.a
        nq = 2*M+1
        A = np.zeros((N*nq, N*nq), complex)
        qs = self.q
        # precompute plane-wave translation pieces
        fp = np.array([[_fplus(n, kx, kz, k) for n in qs]])[0]   # (nq, nm)
        fm = np.array([[_fminus(n, kx, kz, k) for n in qs]])[0]
        up = np.array([(1j)**qq*((kx - 1j*kz)/k)**qq for qq in qs])   # local exp of upgoing wave
        dn = np.array([(1j)**qq*((kx + 1j*kz)/k)**qq for qq in qs])   # downgoing
        w = 2/(a*kz)
        for l0 in range(N):
            for l in range(N):
                blk = np.zeros((nq, nq), complex)
                if l == l0:
                    for iq, qq in enumerate(qs):
                        for i_n, n in enumerate(qs):
                            blk[iq, i_n] += S.get(n-qq, 0j)
                else:
                    dz = self.z[l0]-self.z[l]
                    if dz > 0:
                        blk += np.einsum('qm,nm,m->qn', up, fp, w*np.exp(1j*kz*dz))
                    else:
                        blk += np.einsum('qm,nm,m->qn', dn, fm, w*np.exp(1j*kz*(-dz)))
                # image of layer l (at -z_l), coefficients c'_n = -(-1)^n c_{l,-n}
                dzi = self.z[l0]+self.z[l]
                T = np.einsum('qm,nm,m->qn', up, fp, w*np.exp(1j*kz*dzi))   # acts on c'
                # c'_n = -(-1)^n c_{-n}:  T @ c' = T @ P @ c with P[n, -n] = -(-1)^n
                P = np.zeros((nq, nq))
                for i_n, n in enumerate(qs):
                    P[i_n, list(qs).index(-n)] = -(-1.0)**n
                blk += T @ P
                A[l0*nq:(l0+1)*nq, l*nq:(l+1)*nq] = blk
        return A

    def solve(self, nu, Zs, mode='axion'):
        k, Ht, Jt, S, kx, kz = self._setup(nu, Zs)
        N, M = self.N, self.M; nq = 2*M+1; qs = self.q
        A = self._matrix(k, Ht, Jt, S, kx, kz)
        # equations:  c_q Ht_q + Jt_q (A c)_q = rhs_q
        D = np.diag(np.tile(Ht, N)); Jd = np.diag(np.tile(Jt, N))
        Mat = D + Jd @ A
        rhs = np.zeros(N*nq, complex); b = np.zeros(N*nq, complex)
        for l0 in range(N):
            z0 = self.z[l0]
            if mode == 'axion':
                b[l0*nq:(l0+1)*nq] = np.exp(1j*k*z0)          # E0 e^{ikz}: all J_q coeffs = e^{ikz0}
                rhs[l0*nq + M] += 1.0                          # uniform -E0 in q = 0
            else:   # passive: incident e^{-ikz} plus mirror reflection -e^{ikz}
                b[l0*nq:(l0+1)*nq] = np.exp(-1j*k*z0)*(-1.0)**qs - np.exp(1j*k*z0)
        rhs = rhs - np.tile(Jt, N)*b
        cvec = np.linalg.solve(Mat, rhs)
        C = cvec.reshape(N, nq)
        # zero-order outgoing amplitude
        amp = 0j
        for l in range(N):
            for i_n, n in enumerate(qs):
                cn = C[l, i_n]; cimg = -(-1.0)**n*C[l, list(qs).index(-n)]
                amp += (2/(self.a*k))*(cn*np.exp(-1j*k*self.z[l]) + cimg*np.exp(1j*k*self.z[l]))
        out = (1.0 + amp) if mode == 'axion' else (-1.0 + amp)
        self._last = dict(k=k, C=C, A=A, b=b.reshape(N, nq), mode=mode, kx=kx, kz=kz)
        return out

    def field(self, X, Z):
        """Total E_y on points (X, Z) inside the stack region (|z - z_l| <= spacing/2 cells),
        from local multipole expansions of the last solution (axion mode includes -E0)."""
        L = self._last; k = L['k']; C = L['C']; qs = self.q
        areg = (L['A'] @ C.reshape(-1)).reshape(C.shape) + L['b']
        E = np.full(X.shape, np.nan, complex)
        for idx in np.ndindex(X.shape):
            x, z = X[idx], Z[idx]
            x = (x + self.a/2) % self.a - self.a/2
            l0 = int(np.argmin(np.abs(self.z - z)))
            rho = np.hypot(x, z - self.z[l0]); th = np.arctan2(z - self.z[l0], x)
            if rho < self.r: continue
            val = np.sum((C[l0]*sp.hankel1(qs, k*rho) + areg[l0]*sp.jv(qs, k*rho))*np.exp(1j*qs*th))
            E[idx] = val + (-1.0 if L['mode'] == 'axion' else 0.0)
        return E
