import numpy as np
eta0 = 4e-7*np.pi*299792458.0
def peak_q(f, b):
    i = int(np.argmax(b)); P = b[i]
    jl = i
    while jl > 0 and b[jl] > P/2: jl -= 1
    jh = i
    while jh < len(b)-1 and b[jh] > P/2: jh += 1
    fl = np.interp(P/2, [b[jl], b[jl+1]], [f[jl], f[jl+1]])
    fh = np.interp(P/2, [b[jh], b[jh-1]], [f[jh], f[jh-1]])
    return f[i], P, f[i]/(fh-fl)
def scanI(f, b2):
    return np.trapezoid(b2**2, f)/np.mean(f)
def finite_metrics(R):
    f = R['freqs']; D = float(R['D']); j0 = int(R['j0']); j1 = int(R['j1'])
    E = R['accE']/R['E0'][:, None]; H = R['accH']/R['E0'][:, None]
    Ny = E.shape[1]; y = (np.arange(Ny) + 0.5)*D; W = (j1-j0)*D; yc = 0.5*(j0+j1)*D
    flux = np.sum(np.real(-E*np.conj(H)), axis=1)*D*eta0          # = int |beta|^2 dy for a plane wave
    b2_tot = flux/W
    w = W/2; u = np.exp(-((y-yc)/w)**2)
    cG = (E*u[None, :]).sum(1)*D
    b2_G = np.abs(cG)**2/((u**2).sum()*D*W)
    b2_mid = np.abs(E[:, (j0+j1)//2])**2
    return dict(f=f, b2_tot=b2_tot, b2_G=b2_G, b2_mid=b2_mid, y=y-yc, E=E, W=W)
