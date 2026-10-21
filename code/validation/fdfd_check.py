"""
Independent finite-difference frequency-domain (FDFD) check of the axion boost of a
mirror-backed stack of perfectly conducting wire grids (lossless), meshing the wires
explicitly.  Unknown: E_s = E + E0 (E0 = 1), which satisfies the homogeneous Helmholtz
equation in vacuum, E_s = 1 on every conductor (total E = 0), and is outgoing above the
stack (PML).  Periodic along x with period a.  |beta| = |<E_s>_x| above the stack.
"""
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spla

def fdfd_beta(nu, a, r, zl, dx, buf=1.0, pml=0.6, c0=2.99792458e8):
    k = 2*np.pi*nu/c0
    Nx = int(round(a/dx)); dx = a/Nx
    Ltop = zl[-1] + buf*a; Lz = Ltop + pml*a
    Nz = int(round(Lz/dx)); dz = Lz/Nz
    x = (np.arange(Nx) + 0.5)*dx - a/2
    z = np.arange(1, Nz)*dz           # interior nodes; z=0 (mirror) and z=Lz Dirichlet
    nzi = len(z)
    # PML stretching s(z) on nodes and half nodes
    def sig(zz):
        t = np.clip((zz - Ltop)/(pml*a), 0, None)
        return 40.0*t**3
    zh = (np.arange(Nz) + 0.5)*dz
    s_n = 1 + 1j*sig(z); s_h = 1 + 1j*sig(zh)
    # 1D z operator (1/s_n) d/dz (1/s_h) d/dz on interior nodes
    main = -(1/s_h[1:] + 1/s_h[:-1])/(s_n*dz**2)
    up = (1/s_h[1:-1])/(s_n[:-1]*dz**2)
    lo = (1/s_h[1:-1])/(s_n[1:]*dz**2)
    Dz = sps.diags([lo, main, up], [-1, 0, 1], shape=(nzi, nzi), format='csr')
    # periodic x second difference
    ex = np.ones(Nx)
    Dx = sps.diags([ex[:-1], -2*ex, ex[:-1]], [-1, 0, 1], shape=(Nx, Nx), format='lil')
    Dx[0, -1] = 1; Dx[-1, 0] = 1; Dx = Dx.tocsr()/dx**2
    Lap = sps.kron(sps.identity(nzi), Dx) + sps.kron(Dz, sps.identity(Nx))
    A = (Lap + k**2*sps.identity(nzi*Nx)).tolil()
    rhs = np.zeros(nzi*Nx, complex)
    # mirror Dirichlet E_s(0) = 1 enters the first interior row
    first = np.arange(Nx)
    rhs[first] -= (1/s_h[0])/(s_n[0]*dz**2)*1.0
    # wires: Dirichlet E_s = 1 at nodes inside conductors
    X, Zg = np.meshgrid(x, z)
    inside = np.zeros_like(X, bool)
    for zc in zl:
        inside |= (X**2 + (Zg - zc)**2) <= r**2
    idx = np.flatnonzero(inside.ravel())
    A = A.tocsr()
    # replace rows
    A = A.tolil()
    for i in idx:
        A.rows[i] = [i]; A.data[i] = [1.0]
    A = A.tocsr()
    rhs[idx] = 1.0
    E = spla.spsolve(A.tocsc(), rhs).reshape(nzi, Nx)
    iz = np.argmin(np.abs(z - (zl[-1] + 0.6*a)))
    return abs(E[iz].mean()), E, x, z
