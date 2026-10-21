"""
2D effective-medium model of a port-(antenna-)coupled wire-medium resonator at the plasma frequency.

Envelope equation (E || wires || B, uniform along wires between conducting end plates):
    lap_perp E + k0^2 eps(x) E = -k0^2 E0,   eps = delta + i[1/Qu + gamma_p * 1_port(x)]
delta = 1 - wp^2/w^2 is the tuning variable.  Near wp, a wire medium has w^2 = wp^2 + c^2 k^2 for
k perpendicular AND parallel to the wires (longitudinal plasmon), so the same envelope equation holds
for transport along the wires; here ports are "column" (line) sinks spanning the column length.
Lateral walls: Dirichlet (metal walls, E_z = 0) -> fundamental mode sin*sin with C = 64/pi^4.
Unit-cell mode: Neumann (mirror) boundaries = infinite periodic lattice of ports.

Ports are absorbing patches (side a_p) whose extra loss gamma_p is optimized; absorbed power in the
patches is the extracted signal.  Scan metric: F = int P_port(delta)^2 d delta, normalized to the
audit's ideal single-mode port model, Eq. (1) at zeta = 2 with C = 64/pi^4, same volume:
    F_ref = (2 pi / 27) C^2 Qu V^2      (per unit column length, V = area)
eta_ext = sqrt(F/F_ref)  (power-efficiency equivalent; can exceed 1 because distributed loading of a
large medium approaches C_eff -> 1 > 64/pi^4).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

C_SQ = 64/np.pi**4

def lap2d(n, h, bc):
    """5-point Laplacian on n x n cell-centred grid. bc='D' (Dirichlet at walls) or 'N' (Neumann)."""
    main = -2*np.ones(n); off = np.ones(n-1)
    if bc == 'D':
        pass                      # ghost = -interior (cell-centred Dirichlet): diag -3 at ends
        main[0] = main[-1] = -3
    else:
        main[0] = main[-1] = -1   # ghost = interior (Neumann)
    T = sp.diags([off, main, off], [-1, 0, 1], format='csr')/h**2
    I = sp.identity(n, format='csr')
    return (sp.kron(T, I) + sp.kron(I, T)).tocsc()

def port_mask(n, h, centers, a_p):
    x = (np.arange(n)+0.5)*h
    X, Y = np.meshgrid(x, x, indexing='ij')
    m = np.zeros((n, n), bool)
    for (cx, cy) in centers:
        m |= (np.abs(X-cx) <= a_p/2+1e-12) & (np.abs(Y-cy) <= a_p/2+1e-12)
    return m.ravel()

def lattice_centers(W, nside):
    s = W/nside
    xs = (np.arange(nside)+0.5)*s
    return [(x, y) for x in xs for y in xs]

def solve_F(W, h, lam, Qu, mask, gamma_p, bc='D', uniform=False, ndelta=161, span=None):
    """Return F = int P_port^2 d delta, and the spectrum."""
    n = int(round(W/h)); h = W/n
    k0 = 2*np.pi/lam
    L = lap2d(n, h, bc)
    loss = np.full(n*n, 1/Qu)
    if uniform:
        loss = loss + gamma_p
        pm = np.ones(n*n, bool)
    else:
        loss = loss + gamma_p*mask
        pm = mask
    A0 = L + sp.diags(1j*k0**2*loss)
    rhs = -k0**2*np.ones(n*n)
    # resonance centre estimate: fundamental mode shift
    if bc == 'D':
        d0 = 2*(np.pi/(k0*W))**2
    else:
        d0 = 0.0
    Gt = 1/Qu + gamma_p*(pm.mean())
    if span is None:
        span = 25*max(Gt, 1/Qu)
    # nonuniform grid denser near d0
    u = np.linspace(-1, 1, ndelta)
    deltas = d0 + span*np.sign(u)*np.abs(u)**2
    P = np.empty(ndelta)
    I = sp.identity(n*n, format='csc')
    for i, dl in enumerate(deltas):
        A = (A0 + k0**2*dl*I).tocsc()
        E = spla.spsolve(A, rhs)
        P[i] = np.sum((gamma_p*pm)*np.abs(E)**2)*h*h
    F = np.trapezoid(P**2, deltas)
    return F, deltas, P

def lap_quarter(n, h):
    """Quarter box [0,W/2]^2: Dirichlet at x=0 / y=0 (metal walls), Neumann at x=W/2 / y=W/2 (symmetry)."""
    main = -2*np.ones(n); off = np.ones(n-1)
    main[0] = -3; main[-1] = -1
    T = sp.diags([off, main, off], [-1, 0, 1], format='csr')/h**2
    I = sp.identity(n, format='csr')
    return (sp.kron(T, I) + sp.kron(I, T)).tocsc()

def solve_F_box_sym(W, h, lam, Qu, nside, a_p, gamma_p, ndelta=121, span=None, uniform=False):
    """Full Dirichlet box of side W with an nside x nside symmetric port lattice, solved on a quarter."""
    Wq = W/2; n = int(round(Wq/h)); h = Wq/n
    k0 = 2*np.pi/lam
    L = lap_quarter(n, h)
    x = (np.arange(n)+0.5)*h
    X, Y = np.meshgrid(x, x, indexing='ij')
    if uniform:
        pm = np.ones(n*n, bool)
    else:
        s = W/nside; cs = (np.arange(nside)+0.5)*s
        m = np.zeros((n, n), bool)
        for cx in cs:
            for cy in cs:
                m |= (np.abs(X-cx) <= a_p/2+1e-12) & (np.abs(Y-cy) <= a_p/2+1e-12)
        pm = m.ravel()
    loss = 1/Qu + gamma_p*pm
    A0 = L + sp.diags(1j*k0**2*loss)
    rhs = -k0**2*np.ones(n*n)
    d0 = 2*(np.pi/(k0*W))**2
    Gt = 1/Qu + gamma_p*pm.mean()
    if span is None:
        span = 25*max(Gt, 1/Qu)
    u = np.linspace(-1, 1, ndelta)
    deltas = d0 + span*np.sign(u)*np.abs(u)**2
    P = np.empty(ndelta)
    I = sp.identity(n*n, format='csc')
    for i, dl in enumerate(deltas):
        E = spla.spsolve((A0 + k0**2*dl*I).tocsc(), rhs)
        P[i] = 4*np.sum((gamma_p*pm)*np.abs(E)**2)*h*h
    return np.trapezoid(P**2, deltas), deltas, P

def F_ref(V, Qu, C=C_SQ):
    return 2*np.pi/27*C**2*Qu*V**2

def optimise_gamma(func, lo=-4, hi=1.5, n=12):
    """1D log-scan + golden refinement of F over gamma_p; func(g)->F."""
    gs = np.logspace(lo, hi, n)
    Fs = np.array([func(g) for g in gs])
    i = int(np.argmax(Fs))
    a = np.log10(gs[max(i-1, 0)]); b = np.log10(gs[min(i+1, n-1)])
    phi = (np.sqrt(5)-1)/2
    x1 = b-phi*(b-a); x2 = a+phi*(b-a); f1 = func(10**x1); f2 = func(10**x2)
    for _ in range(8):
        if f1 > f2:
            b, x2, f2 = x2, x1, f1; x1 = b-phi*(b-a); f1 = func(10**x1)
        else:
            a, x1, f1 = x1, x2, f2; x2 = a+phi*(b-a); f2 = func(10**x2)
    best = max([(Fs.max(), gs[i]), (f1, 10**x1), (f2, 10**x2)])
    return best  # (F, gamma)
