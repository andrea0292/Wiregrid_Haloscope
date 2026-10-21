"""Independent recomputation of SQWARE (Mehrani et al., arXiv:2509.14320v2, PRL 136, 241001) benchmarks.

Model (as stated in the paper, re-implemented from scratch):
  eps_CRI(QW) = eps_bg + i*sigma/(eps0*omega),  sigma = n_e e^2 tau /(m* d_QW) / (1 - i(omega+omega_c) tau)   [Eq. 4 / S16,S18]
  eps_eff     = (eps_B d_B + eps_CRI d_QW)/(d_B + d_QW)                                                    [Eq. 5 / S21]
  omega_c     = e B cos(theta)/m*
  beta (open slab)    : solved by our own 1D boundary-value problem and compared with Eq. (3)/(S23)
  beta (mirror-backed): solved by our own 1D BVP and compared with Eq. (S24)
Outputs: sqware_results.csv, sqware_scan_<cfg>.csv, printed summary.
"""
import numpy as np, pandas as pd
from scipy.optimize import minimize_scalar, brentq

# ---------------- constants (SI) ----------------
e = 1.602176634e-19; me = 9.1093837e-31; eps0 = 8.8541878128e-12; c = 2.99792458e8
hbar = 1.054571817e-34; hbar_eVs = 6.582119569e-16
mstar = 0.067 * me
EPS_BG, EPS_B = 11.56, 10.6929          # paper SM Sec. I C, I D (room-T values they adopt)
n_e, d_QW = 3e11 * 1e4, 30e-9           # m^-2, m
meV = 1e-3 / hbar_eVs                   # rad/s per meV

CFG = {  # Table I
    1: dict(D=2e-3, W=3e-2, dB=90e-9, B=36., tau=1.7e-9, Gd=1e-3, eta=0.07, tobs=84*60, rng=(2., 5.)),
    2: dict(D=10e-3, W=3e-2, dB=150e-9, B=36., tau=1.7e-9, Gd=1e-3, eta=0.20, tobs=34*60, rng=(0.6, 2.)),
    3: dict(D=20e-3, W=5e-2, dB=90e-9, B=50., tau=4.0e-9, Gd=1e-4, eta=0.35, tobs=7*60, rng=(0.6, 2.)),
}
QUOTED = {1: dict(beta=70, m=5.0, th=86, Q=6e3), 2: dict(beta=115, m=1.6, th=72, Q=1e4), 3: dict(beta=390, m=1.1, th=60, Q=5e4)}

def wc(B, th):  return e * B * np.cos(np.radians(th)) / mstar

def eps_qw(w, B, th, tau):
    sig = n_e * e**2 * tau / (mstar * d_QW) / (1 - 1j * (w + wc(B, th)) * tau)
    return EPS_BG + 1j * sig / (eps0 * w)

def eps_eff(w, p, th):
    return (EPS_B * p['dB'] + eps_qw(w, p['B'], th, p['tau']) * d_QW) / (p['dB'] + d_QW)

def Omega2(p):   # lossless ENZ condition omega(omega+omega_c) = Omega^2
    return n_e * e**2 / (eps0 * mstar * (EPS_B * p['dB'] + EPS_BG * d_QW))

def w_enz(p, th):
    a = wc(p['B'], th); return (-a + np.sqrt(a * a + 4 * Omega2(p))) / 2

def th_enz(p, w):  # tilt angle whose lossless ENZ sits at omega
    cth = (Omega2(p) / w - w) / (e * p['B'] / mstar)
    return np.degrees(np.arccos(cth)) if 0 <= cth <= 1 else np.nan

# ---------------- boost factors ----------------
def beta_paper(eps, w, D, th):            # Eq. (3) == Eq. (S23)
    n = np.sqrt(eps); Dl = w * n * D / c
    return np.abs(-(1 - eps) * np.sin(Dl / 2) / (eps * np.sin(Dl / 2) + 1j * n * np.cos(Dl / 2))) * np.sin(np.radians(th)) / np.sqrt(2)

def beta_paper_mirror(eps, w, D, th):     # Eq. (S24)
    n = np.sqrt(eps); Dl = w * n * D / c
    return np.abs(1 / eps - 1 - (1 / eps) / np.cos(Dl)) / np.abs(1 - 1j / n * np.tan(Dl)) * np.sin(np.radians(th)) / np.sqrt(2)

def beta_bvp(eps, w, D, th, mirror=False):
    """Own 1D solution: slab 0<z<D, source E_a = -E0/eps inside, -E0 outside (E0=1).
    Unknowns: inside A e^{ikz}+B e^{-ikz}; right R e^{ik0(z-D)}; left L e^{-ik0 z} (or PEC at z=0).
    Continuity of E and dE/dz.  Returns |R| * sin(th)/sqrt(2)."""
    k0 = w / c; k = k0 * np.sqrt(eps); ea = -1 / eps; eo = -1.0
    ep, em = np.exp(1j * k * D), np.exp(-1j * k * D)
    if not mirror:
        M = np.array([[1, 1, 0, -1], [1j*k, -1j*k, 0, 1j*k0], [ep, em, -1, 0], [1j*k*ep, -1j*k*em, -1j*k0, 0]], complex)
        rhs = np.array([eo - ea, 0, eo - ea, 0], complex)
        A, B, R, L = np.linalg.solve(M, rhs)
    else:
        M = np.array([[1, 1, 0], [ep, em, -1], [1j*k*ep, -1j*k*em, -1j*k0]], complex)
        rhs = np.array([-ea, eo - ea, 0], complex)
        A, B, R = np.linalg.solve(M, rhs)
    return np.abs(R) * np.sin(np.radians(th)) / np.sqrt(2)

# ---------------- resonance finder ----------------
def resonance(p, th, mirror=False):
    """Peak of beta(omega) at fixed theta in the fundamental (M=1) window above ENZ, its FWHM in beta."""
    w0 = w_enz(p, th)
    f = (lambda w: beta_paper_mirror(eps_eff(w, p, th), w, p['D'], th)) if mirror else (lambda w: beta_paper(eps_eff(w, p, th), w, p['D'], th))
    ws = w0 * (1 + np.linspace(-2e-3, 0.08, 40001)); bs = np.array([f(x) for x in ws])
    # restrict to first lobe: Re(eps) between 0 and (2*lambda/2D)^2 (below the M=3 resonance)
    re = np.real(eps_eff(ws, p, th)); lam = 2 * np.pi * c / ws
    ok = (re < (1.5 * lam / p['D'])**2) if not mirror else (re < (1.0 * lam / p['D'])**2)
    i = np.argmax(np.where(ok, bs, 0))
    r = minimize_scalar(lambda x: -f(x), bracket=(ws[max(i-1,0)], ws[i], ws[min(i+1,len(ws)-1)]))
    wp, bp = r.x, -r.fun
    try:
        lo = brentq(lambda x: f(x) - bp / 2, ws[max(i - 4000, 0)], wp)
        hi = brentq(lambda x: f(x) - bp / 2, wp, ws[min(i + 4000, len(ws) - 1)])
        Q = wp / (hi - lo)
    except ValueError:
        Q = np.nan
    eps_p = eps_eff(wp, p, th)
    return dict(theta=th, m_meV=wp / meV, beta=bp, Q=Q, eps_re=eps_p.real, eps_im=eps_p.imag,
                wtau=wp * p['tau'], wpwc_tau=(wp + wc(p['B'], th)) * p['tau'], m_enz_meV=w0 / meV)

# ---------------- signal / sensitivity (natural units, eV) ----------------
T_eV2 = 195.353            # 1 T in eV^2 (Heaviside-Lorentz)
m_eVinv = 5.0677307e6       # 1 m in eV^-1
eV2_W = e / hbar_eVs        # 1 eV^2 -> W   (eV * eV/hbar)
def rho_eV4(rho_GeVcm3): return rho_GeVcm3 * 1e9 * (1.973269804e-5) ** 3

def gKSVZ(m_meV, EN=0.0):   # GeV^-1 ; m_a = 5.69 ueV (1e12 GeV/f_a)
    fa = 5.691e9 / m_meV    # GeV
    return (1 / 137.036) / (2 * np.pi) / fa * abs(EN - 1.924)

def photon_rate(g_GeV, B, A_m2, beta, m_meV, eta, rho=0.3, a0sq_factor=2.0):
    """Gamma = eta * beta^2 g^2 B^2 a0^2 A /(2 omega), a0^2 = a0sq_factor*rho/m^2  (paper: factor 2; their quoted 1.52 eV <-> factor 1)."""
    g = g_GeV * 1e-9; m = m_meV * 1e-3
    a0sq = a0sq_factor * rho_eV4(rho) / m**2
    P = beta**2 * g**2 * (B * T_eV2)**2 * a0sq * (A_m2 * m_eVinv**2) / 2      # eV^2
    return eta * P / m / hbar_eVs                                             # s^-1

def g_limit_bkg(B, A, beta, m_meV, eta, Gd, tobs, rho=0.3, a0f=2.0):
    need = 2 * np.sqrt(Gd / tobs)
    return np.sqrt(need / photon_rate(1.0, B, A, beta, m_meV, eta, rho, a0f))

def g_limit_free(B, A, beta, m_meV, eta, tobs, rho=0.3, a0f=2.0, N=3):
    return np.sqrt(N / tobs / photon_rate(1.0, B, A, beta, m_meV, eta, rho, a0f))

if __name__ == '__main__':
    rows = []
    # --- 0. check our own BVP against Eq. (3)/(S23) and (S24)
    p = CFG[2]; th = 72.437; w = 1.5554 * meV
    for ww in w * (1 + np.array([-1e-3, 0, 1e-4, 1e-3, 1e-2])):
        ee = eps_eff(ww, p, th)
        rows.append(dict(quantity='BVP_vs_Eq3', config=2, m_meV=ww/meV, ours=beta_bvp(ee, ww, p['D'], th), paper=beta_paper(ee, ww, p['D'], th)))
        rows.append(dict(quantity='BVP_vs_S24', config=2, m_meV=ww/meV, ours=beta_bvp(ee, ww, p['D'], th, True), paper=beta_paper_mirror(ee, ww, p['D'], th)))
    pd.DataFrame(rows).to_csv('sqware_formula_check.csv', index=False)
    print(pd.DataFrame(rows).to_string())

    out = []
    for k, p in CFG.items():
        Om = np.sqrt(Omega2(p)); f = d_QW / (d_QW + p['dB'])
        wp2D = n_e * e**2 / (eps0 * mstar * d_QW)      # 3D-equivalent QW plasma freq^2
        rec = dict(config=k, N_wells=p['D'] / (d_QW + p['dB']), growth_days_1um_per_h=p['D'] / 1e-6 / 24,
                   fill=f, hw_p_QW_meV=np.sqrt(wp2D) * hbar_eVs * 1e3, hw_p_QW_screened_meV=np.sqrt(wp2D / EPS_BG) * hbar_eVs * 1e3,
                   hOmega_ENZmax_theta90_meV=Om * hbar_eVs * 1e3, m_ENZ_theta0_meV=w_enz(p, 0) / meV,
                   hwc_full_meV=e * p['B'] / mstar * hbar_eVs * 1e3,
                   theta_lo=th_enz(p, p['rng'][0] * meV), theta_hi=th_enz(p, p['rng'][1] * meV),
                   B2A_full=p['B']**2 * p['W']**2)
        rec['B2A_par_lo'] = (p['B'] * np.sin(np.radians(rec['theta_lo'])))**2 * p['W']**2
        rec['B2A_par_hi'] = (p['B'] * np.sin(np.radians(rec['theta_hi'])))**2 * p['W']**2
        # scan theta
        scan = [resonance(p, t) for t in np.arange(5, 89.9, 0.5)]
        scm = [resonance(p, t, mirror=True) for t in np.arange(5, 89.9, 0.5)]
        S = pd.DataFrame(scan); Sm = pd.DataFrame(scm).add_suffix('_mirror')
        S = pd.concat([S, Sm[['beta_mirror', 'Q_mirror', 'm_meV_mirror']]], axis=1)
        S['B2A_par'] = (p['B'] * np.sin(np.radians(S.theta)))**2 * p['W']**2
        S['g_lim_bkg'] = [g_limit_bkg(p['B'], p['W']**2, b, m, p['eta'], p['Gd'], p['tobs']) for b, m in zip(S.beta, S.m_meV)]
        S['g_lim_free3'] = [g_limit_free(p['B'], p['W']**2, b, m, p['eta'], p['tobs']) for b, m in zip(S.beta, S.m_meV)]
        S['g_KSVZ'] = gKSVZ(S.m_meV); S['g_DFSZ'] = gKSVZ(S.m_meV, 8/3)
        S.to_csv(f'sqware_scan_cfg{k}.csv', index=False)
        imax = S.beta.idxmax()
        inr = S[(S.m_meV >= p['rng'][0]) & (S.m_meV <= p['rng'][1])]
        rec.update(beta_max=S.beta[imax], m_at_betamax=S.m_meV[imax], th_at_betamax=S.theta[imax], Q_max=S.Q.max(),
                   beta_max_mirror=S.beta_mirror.max(), m_at_betamax_mirror=S.m_meV_mirror[S.beta_mirror.idxmax()],
                   beta_range_min=inr.beta.min(), beta_range_max=inr.beta.max(), Q_range_min=inr.Q.min(), Q_range_max=inr.Q.max(),
                   wtau_range_lo=p['rng'][0]*meV*p['tau'], wtau_range_hi=p['rng'][1]*meV*p['tau'],
                   power_gain_beta2=S.beta[imax]**2,
                   eps_im_at_peak=S.eps_im[imax], eps_re_at_peak=S.eps_re[imax])
        # quoted-point checks
        q = QUOTED[k]; r = resonance(p, q['th'])
        rec.update(quoted_beta=q['beta'], quoted_m=q['m'], quoted_th=q['th'], quoted_Q=q['Q'],
                   ours_at_quoted_th_m=r['m_meV'], ours_at_quoted_th_beta=r['beta'], ours_at_quoted_th_Q=r['Q'])
        # signal rates at beta_max point (their beta & ours), KSVZ, rho 0.3 & 0.4
        m0 = q['m']; A = p['W']**2
        for lab, b in [('theirbeta', q['beta']), ('ourbeta', S.beta[imax])]:
            for rho in (0.3, 0.4):
                G = photon_rate(gKSVZ(m0), p['B'], A, b, m0, p['eta'], rho)
                rec[f'Gsig_KSVZ_{lab}_rho{rho}_Hz'] = G
                rec[f'tneed_bkg_KSVZ_{lab}_rho{rho}_s'] = 4 * p['Gd'] / G**2     # S > 2 sqrt(B)
                rec[f'tneed_3ph_KSVZ_{lab}_rho{rho}_s'] = 3 / G                  # 3 signal photons
            rec[f'Gsig_DFSZ_{lab}_rho0.3_Hz'] = photon_rate(gKSVZ(m0, 8/3), p['B'], A, b, m0, p['eta'], 0.3)
            rec[f'g95_bkg_{lab}'] = g_limit_bkg(p['B'], A, b, m0, p['eta'], p['Gd'], p['tobs'])
            rec[f'g95_free3_{lab}'] = g_limit_free(p['B'], A, b, m0, p['eta'], p['tobs'])
            rec[f'g95_bkg_{lab}_a0eq1p52'] = g_limit_bkg(p['B'], A, b, m0, p['eta'], p['Gd'], p['tobs'], a0f=1.0)
        rec['g_KSVZ_at_m'] = gKSVZ(m0); rec['g_DFSZ_at_m'] = gKSVZ(m0, 8/3)
        rec['dark_counts_per_tobs'] = p['Gd'] * p['tobs']
        # Eq.(6) normalisation check (beta=100, B=36, eta=.2, A=9cm2, Gd=1mHz, m=1meV, tobs=30 d)
        rec['eq6_norm_a0sq_2rho'] = g_limit_bkg(36, 9e-4, 100, 1.0, 0.2, 1e-3, 30*86400)
        rec['eq6_norm_a0_1p52'] = g_limit_bkg(36, 9e-4, 100, 1.0, 0.2, 1e-3, 30*86400, a0f=1.0)
        # scan bookkeeping: steps of one FWHM at Q_max
        dm = q['m'] / q['Q']; nsteps = (p['rng'][1] - p['rng'][0]) / dm
        rec['n_steps_at_quotedQ'] = nsteps; rec['ttot_days_at_quotedQ'] = nsteps * p['tobs'] / 86400
        rec['n_steps_logQ_ours'] = np.trapezoid(np.interp(np.linspace(*p['rng'], 400), S.m_meV, S.Q) / np.linspace(*p['rng'], 400), np.linspace(*p['rng'], 400))
        rec['ttot_days_ours'] = rec['n_steps_logQ_ours'] * p['tobs'] / 86400
        out.append(rec)
    R = pd.DataFrame(out).set_index('config').T
    R.to_csv('sqware_results.csv')
    pd.set_option('display.width', 200); print(R.to_string())
    # a0 check
    print('a0 at 1 meV, rho=0.3: sqrt(2rho)/m =', np.sqrt(2 * rho_eV4(0.3)) / 1e-3, ' sqrt(rho)/m =', np.sqrt(rho_eV4(0.3)) / 1e-3)
