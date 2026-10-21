"""
Driver reproducing the finite-tile numbers (decision memo Sec. 1, Letter Appendix F).
  python run_finite_tile.py [--quick] [--with20]
--quick : 3 cm tile only (~2 min).  Default: 1, 3, 10 cm + endpoints + tilts (~45 min, ~4 GB RAM).
--with20: also the 20 cm tile (~11k unknowns, ~13 GB RAM, ~30 min).
Writes finite_tile_repro.json.
"""
import sys, os, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from finite_tile_analysis import *

quick = '--quick' in sys.argv; with20 = '--with20' in sys.argv
out = {}; T0 = time.time()
def log(*x): print(f'[{time.time()-T0:7.0f}s]', *x, flush=True)

# 1) tiles: forward boost, peak and Q_resp
widths = [0.03] if quick else [0.01, 0.03, 0.10] + ([0.20] if with20 else [])
tiles = {}
for W in widths:
    t = make_tile(dR, W=W); f0, pk, Q = tile_peak(t, 0.996*nuT, 1.004*nuT, 9)
    tiles[W] = (t, f0); out[f'W={W}'] = dict(f_over_nuT=f0/nuT, beta2_fwd=pk, Q_resp=Q)
    log('W', W, out[f'W={W}'])
# 2) edge fit beta2_inf (1 - delta/W): report the fit range explicitly
def edge_fit(Ws):
    b = np.array([out[f'W={w}']['beta2_fwd'] for w in Ws]); A = np.c_[np.ones(len(Ws)), -1/np.array(Ws)]
    c, *_ = np.linalg.lstsq(A, b, rcond=None); return dict(fit_widths=Ws, beta2_inf=c[0], delta_m=c[1]/c[0])
if not quick:
    out['edge_fit_3_10'] = edge_fit([0.03, 0.10]); out['edge_fit_1_3_10'] = edge_fit([0.01, 0.03, 0.10])
    if with20: out['edge_fit_3_10_20'] = edge_fit([0.03, 0.10, 0.20])
    log(out['edge_fit_3_10'], out['edge_fit_1_3_10'])
# 3) mirror margin and side walls (3 cm)
for lab, kw in [('margin_1cm', dict(Wm=None, mirror_margin=0.01)), ('side_walls', dict(side_walls=True))]:
    t = make_tile(dR, W=0.03, **kw); f0, pk, Q = tile_peak(t, 0.996*nuT, 1.004*nuT, 9)
    out[lab] = dict(f_over_nuT=f0/nuT, beta2_fwd=pk, Q_resp=Q); log(lab, out[lab])
# 4) end-wall and mirror losses relative to wire loss (3 cm, centre wire of layer 7)
t3, f3 = tiles[0.03]; t3.factor(f3, Zs(f3)); q3 = t3.solve_axion(0.0)[:, 0]; Nw = t3.Nw
jc = 6*Nw + Nw//2; xc, zc = t3.wx[jc], t3.wz[jc]
rr = np.geomspace(r, 0.5*min(a, dR), 200); tt = np.linspace(0, 2*np.pi, 128, endpoint=False)
RR, TT = np.meshgrid(rr, tt); Hx, Hz = H_field(t3, q3, xc + RR*np.cos(TT), zc + RR*np.sin(TT))
I_disk = np.trapezoid(np.trapezoid((abs(Hx)**2 + abs(Hz)**2).reshape(RR.shape)*RR, tt, axis=0), rr)
gx = np.linspace(xc - a/2, xc + a/2, 161); gz = np.linspace(zc - dR/2, zc + dR/2, 161); GX, GZ = np.meshgrid(gx, gz)
msk = np.hypot(GX - xc, GZ - zc) > 0.5*min(a, dR); Hx2, Hz2 = H_field(t3, q3, GX[msk], GZ[msk])
cell = I_disk + np.sum(abs(Hx2)**2 + abs(Hz2)**2)*(gx[1] - gx[0])*(gz[1] - gz[0])
Hxs, Hzs = H_field(t3, q3, xc + 1.0001*r*np.cos(tt), zc + 1.0001*r*np.sin(tt)); Hs2 = np.mean(abs(Hxs)**2 + abs(Hzs)**2)
Ks = -4*q3[t3.nw:]/(2*np.pi*f3*mu0)
out['losses'] = {f'endwall_over_wire_Ly={Ly}': 2*cell/(Ly*2*np.pi*r*Hs2) for Ly in (0.03, 0.10)}
out['losses']['mirror_over_wire'] = np.sum(0.5*Rsf(f3)*abs(Ks)**2*t3.slen)/t3.wire_loss(q3, Rsf(f3)); log(out['losses'])
# 5) receiver mode, halo, emissivity
for W in [w for w in widths if w >= 0.03]:
    t, f0 = tiles[W]
    out[f'collect_W={W}'] = dict(nowind=collect_chain(t, f0), wind_x=collect_chain(t, f0, vE))
    cc = out[f'collect_W={W}']['nowind']; t.factor(f0, Zs(f0)); k = 2*np.pi*f0/c0
    E, Pinc = gauss_incident(cc['w_over_W']*W, t.wz.max() + 2e-3, k); q = t.solve_incident(E)
    out[f'collect_W={W}']['absorptivity_x_mode'] = t.wire_loss(q, Rsf(f0))/Pinc; log('collect', W, out[f'collect_W={W}'])
if not quick:
    # 6) octave endpoints (10 cm), spacings from the periodic retune
    for m, dv in [(0.8, 1.6731459794546668*a), (1.6, 0.5226354153966059*a)]:
        t = make_tile(dv); tgt = nu_of_m(m*1e-3); f0, pk, Q = tile_peak(t, 0.996*tgt, 1.004*tgt, 9)
        out[f'endpoint_{m}meV'] = dict(f_meV=f0/nu_of_m(1e-3), beta2_fwd=pk, Q_resp=Q, **collect_chain(t, f0)); log(m, out[f'endpoint_{m}meV'])
    # 7) random rigid tilts of each grid (10 cm)
    rng = np.random.default_rng(7); tb = make_tile(dR); out['tilts'] = []
    for se in (1e-6, 3e-6):
        for rep in range(3):
            al = rng.normal(0, se, Nl)/(tb.W/2); dz = (al[:, None]*tb.wx.reshape(Nl, -1)).ravel()
            t = make_tile(dR, dz=dz); f0, pk, Q = tile_peak(t, 0.996*nuT, 1.004*nuT, 9)
            out['tilts'].append(dict(sigma_edge=se, rep=rep, beta2_fwd=pk, Q_resp=Q, **collect_chain(t, f0))); log(out['tilts'][-1])
json.dump(out, open('finite_tile_repro.json', 'w'), indent=1, default=float); log('done')
