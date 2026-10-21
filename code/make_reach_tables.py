"""End-to-end driver for the axion reach tables behind Fig. 3 (and, through dark_photon/dp_tools.py, Fig. 4).

    python make_reach_tables.py            -> ../data/reach_0p25_8_fw_repro.csv, ../data/reach_0p25_8_fw_Rs2_repro.csv
    python make_reach_tables.py --compare  -> also prints the deviation from the stored ../data/reach_0p25_8_fw*.csv
    --overlap-at-4 : also let the 2-4 and 4-8 meV cartridges both contribute within +-5% of 4 meV (as at the other
                     junctions); this lowers g by up to 16% (counters) / 30% (KIDs) at the three grid points nearest 4 meV

Chain (all in this file or in reproduce_final.py):
  1. Designs: four cartridges 0.25-4 meV (r = 0.03 a, N = 15, 13, 11, 9) and the 4-8 meV cartridge (2r = 5 um, N = 8),
     each centred at the geometric band centre with d = a (impedance-sheet model, Rs = RS_SCALE * 83.6 mOhm (nu/250 GHz)^(2/3)).
  2. Tuning trajectories d(m) across each band (build_traj).
  3. Scan response per mass: J_l = int |beta|^4 dnu/nu over the band (heterodyne) and J_c = same within the half-maximum
     window (counters/KIDs), from the sheet model, multiplied by the circular-wire correction (sqrt(I_wire/I_sheet))^2
     interpolated quadratically in log m between the lower, centre and upper settings of data/fullwave_reconciliation*.csv
     (phase sqrt(3), i.e. Z_s = R_s (1 - i sqrt 3)).
  4. Halo factor xi = chi_their(m, w = W/2, L = W): Gaussian receiver mode truncated at the square tile of side W,
     Gaussian halo with 1D dispersion 220/sqrt(2) km/s and the 232 km/s mean wind IN the tile plane.
     Tiles: W = 10 cm up to 4 meV, 5 cm for the 4-8 meV cartridge.
  5. Single-tile signal P1 = 7.33e-26 W * eta_d * eta_ap * xi * (W/m)^2 (meV/m)^2 (g = 1e-12 GeV^-1, 10 T, rho = 0.4),
     eta_d = 0.5, eta_ap = 0.65; tiles add incoherently; SNR = 5.
  Columns (the Letter plots lin10, kid10, target, each multiplied by 10^(1/4) for the 1 m^2 array):
     target : counters, 10 m^2, 4 yr over the 0.25-4 meV design's four octaves (1 live yr per octave), total dark rate 4 /s
              (4 mHz per 10 cm tile, 1 mHz per 5 cm tile), converter at 0.1 K
     kid10  : KIDs, NEP 3e-19 W/rtHz per tile entering as a counter with D = NEP^2/(2 (h nu)^2), same exposure
     lin10  : heterodyne (noise 10 h nu, line width 1e-6 nu), same exposure; = lin_ref / 40^(1/4)
     lin_ref: heterodyne, 1 m^2, 1 yr in total (kept for reference)
"""
import os, sys, numpy as np, pandas as pd
from scipy.optimize import brentq
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
_argv = sys.argv; sys.argv = sys.argv[:1]
import reproduce_final as rf
sys.argv = _argv
DATA = os.path.join(HERE, '..', 'data')
hP, kB, yr = rf.hP, rf.kB, rf.yr
BANDS = [(0.25, 0.5), (0.5, 1), (1, 2), (2, 4)]; NS = [15, 13, 11, 9]
KEY = {(0.25, 0.5): '0.25-0.5', (0.5, 1): '0.5-1', (1, 2): '1-2', (2, 4): '2-4', (4, 8): '4-8'}
TILE = {k: 0.10 for k in KEY}; TILE[(4, 8)] = 0.05

def design_fixed_r(mc_meV, N, r):
    nuc = rf.nu_of_m(mc_meV*1e-3)
    g = lambda a: rf.peak_info(a, r, a, N, nuc, 0.3, 601)[0]/nuc - 1
    a0 = 60e-6*(5.657/mc_meV); a = brentq(g, 0.75*a0, 1.35*a0, xtol=1e-5*a0)
    f0, pk, Q, I = rf.peak_info(a, r, a, N, nuc); return a, f0, pk, Q, I

def designs():
    des = {}
    for (lo, hi), N in zip(BANDS, NS):
        a, f0, pk, Q, I = rf.design_band(np.sqrt(lo*hi), N); des[(lo, hi)] = dict(N=N, a=a, r=0.03*a, b=pk, Q=Q, I=I)
    traj = rf.build_traj(des, [0.25, 0.5, 1, 2, 4])
    a5, f05, pk5, Q5, I5 = design_fixed_r(np.sqrt(32), 8, 2.5e-6)
    des[(4, 8)] = dict(N=8, a=a5, r=2.5e-6, b=pk5, Q=Q5, I=I5)
    traj.update(rf.build_traj({(4, 8): des[(4, 8)]}, [4, 8]))
    return des, traj

def wire_correction():
    R = pd.concat([pd.read_csv(os.path.join(DATA, 'fullwave_reconciliation.csv')),
                   pd.read_csv(os.path.join(DATA, 'fullwave_reconciliation_4_8.csv'))])
    R = R[R.phase > 1.5]; tab = {}
    for band, g in R.groupby('band'):
        lo, hi = [float(x) for x in band.split('-')]
        tab[band] = np.polyfit(np.log(g.m/lo)/np.log(hi/lo), g.sqrtI_retuned**2, 2)
    def corr(key, m):
        lo, hi = key; x = np.clip(np.log(m/lo)/np.log(hi/lo), 0, 1); return float(np.polyval(tab[KEY[key]], x))
    return corr

def J_terms(m_a, des, traj, corr):
    nua = rf.nu_of_m(m_a*1e-3); out = []
    for key, tr in traj.items():
        if not (tr['lo']/1.05 <= m_a <= tr['hi']*1.05): continue
        D = des[key]; ms = np.clip(m_a, tr['lo'], tr['hi'])
        dv = np.exp(np.interp(np.log(ms), np.log(tr['m']), np.log(tr['d'])))
        fpk, pk = rf.peak_fast(D['a'], D['r'], dv, D['N'], rf.nu_of_m(ms*1e-3), 0.03, 121)
        w = fpk/D['Q']; y = np.linspace(-8*w, 8*w, 1601); b = rf.spec_conf(D['a'], D['r'], dv, D['N'], fpk + y)
        half = y[b >= pk/2]; Qr = fpk/(half.max() - half.min())
        nlo, nhi = rf.nu_of_m(tr['lo']*1e-3), rf.nu_of_m(tr['hi']*1e-3); mask = ((nua - y) >= nlo) & ((nua - y) <= nhi)
        Jl = np.trapezoid(np.where(mask, b**2, 0), y)/nua; mc = mask & (np.abs(y) < fpk/(2*Qr))
        Jc = np.trapezoid(np.where(mc, b**2, 0), y)/nua
        cf = corr(key, ms); out.append(dict(key=key, Jl=Jl*cf, Jc=Jc*cf, Qr=Qr, fpk=fpk, b=pk))
    return out

def snr2(m_a, terms, A, T, mode, Dtot=None, Dpix=None, Tc=0.1, eta=0.5, eta_ap=0.65):
    """sum over tiles and settings of SNR^2 at g = 1e-12 GeV^-1."""
    nua = rf.nu_of_m(m_a*1e-3); s = 0.0
    for t in terms:
        W = TILE[t['key']]; nt = A/W**2; ch = rf.chi_their(m_a, w=W/2, L=W); rho = T*yr/(4*np.log(2))
        P1 = 7.33e-26*eta*eta_ap*ch*W**2/m_a**2
        if mode == 'lin':
            s += nt*rho*t['Jl']*P1**2/((10*hP*nua)**2*1e-6*nua)
        else:
            Dd = Dtot/nt if Dpix is None else Dpix(nua)
            th = eta*(t['fpk']/t['Qr'])/np.expm1(min(hP*t['fpk']/(kB*Tc), 700))
            s += nt*rho*t['Jc']/(Dd + th)*(P1/(hP*nua))**2
    return s

g_of = lambda s: 1e-12*(25/s)**0.25
KID = lambda nu: (3e-19)**2/(2*(hP*nu)**2)

def table(ms, des, traj, corr, overlap_at_4=False):
    rows = []
    for m in ms:
        t = J_terms(m, des, traj, corr)
        if not overlap_at_4:
            # convention of the stored tables: the 0.25-4 and 4-8 meV parts were computed separately, so at the
            # 4 meV junction only one side contributes (internal junctions at 0.5, 1, 2 meV do include both neighbours)
            t = [x for x in t if (x['key'] == (4, 8)) == (m >= 4.0 - 1e-9)]
        lin_ref = g_of(snr2(m, t, 1.0, 1.0, 'lin'))
        rows.append(dict(m=m, lin_ref=lin_ref, target=g_of(snr2(m, t, 10.0, 4.0, 'cnt', Dtot=4.0, Tc=0.1)),
                         kid10=g_of(snr2(m, t, 10.0, 4.0, 'cnt', Dpix=KID, Tc=0.1)), lin10=lin_ref/40**0.25))
    return pd.DataFrame(rows)

if __name__ == '__main__':
    ms = pd.read_csv(os.path.join(DATA, 'reach_0p25_8_fw.csv')).m.values        # same 301-point mass grid
    corr = wire_correction(); out = {}
    rf.RS_SCALE = 1.0; des1, traj1 = designs()
    for tag, scale in (('', 1.0), ('_Rs2', 2.0)):
        rf.RS_SCALE = scale
        if scale == 1.0:
            des, traj = des1, traj1
        else:   # same hardware; resonance width and peak re-evaluated at the band centre with the doubled Rs
            des = {k: dict(v) for k, v in des1.items()}
            for k, D in des.items():
                if k == (4, 8): continue          # as in the stored tables: 4-8 meV keeps its Rs x1 window width
                mc = np.sqrt(k[0]*k[1]); f0, pk = rf.peak_fast(D['a'], D['r'], D['a'], D['N'], rf.nu_of_m(mc*1e-3), 0.05, 201)
                ff = np.linspace(f0*0.95, f0*1.05, 4001); bb = rf.spec_conf(D['a'], D['r'], D['a'], D['N'], ff)
                j = bb.argmax(); lo_ = j; hi_ = j
                while bb[lo_] > pk/2: lo_ -= 1
                while bb[hi_] > pk/2: hi_ += 1
                D['Q'] = f0/(ff[hi_] - ff[lo_]); D['b'] = pk
            traj = rf.build_traj({k: v for k, v in des.items() if k != (4, 8)}, [0.25, 0.5, 1, 2, 4])
            traj.update(rf.build_traj({(4, 8): des[(4, 8)]}, [4, 8]))
        T = table(ms, des, traj, corr, overlap_at_4='--overlap-at-4' in sys.argv)
        T.to_csv(os.path.join(DATA, f'reach_0p25_8_fw{tag}_repro.csv'), index=False); out[tag] = T
        print(tag or 'Rs x1', {k: (round(v['a']*1e6, 2), v['N']) for k, v in des.items()}, flush=True)
    rf.RS_SCALE = 1.0
    if '--compare' in sys.argv:
        for tag in ('', '_Rs2'):
            S = pd.read_csv(os.path.join(DATA, f'reach_0p25_8_fw{tag}.csv'))
            for c in ('target', 'kid10', 'lin10', 'lin_ref'):
                r = out[tag][c].values/S[c].values - 1
                print(f'{tag or "Rs x1":6s} {c:8s} max |dev| = {np.nanmax(np.abs(r)):.2e}  (at m = {S.m.values[np.nanargmax(np.abs(r))]:.3f} meV)')
