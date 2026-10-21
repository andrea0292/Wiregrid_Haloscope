"""
comparator_model.py -- radiating mirror-backed wire-grid stack vs port-coupled wire-medium resonator
at m_a = 1.13137 meV (nu = 273.56 GHz), photon coupling.  Counting-limit scan-rate ratio
    R = R0 (V_s/V_p)^2 * F_s / F_p,     F = (eta * f_coh)^2 / N_rx      (audit Eqs. 5-6)
with equal background per receiving channel (D proportional to the number of channels N_rx).

Inputs fixed by the study (not re-derived): stack 13 grids, d = 3.85 mm, sheet beta^2 = 2066, Q_s = 313
(full-wave 2125/312; audit 2083/314); port model Qu = 473, C = 64/pi^4 (audit Eq. 1).
Extraction efficiency eta_ext(s, W) of a lattice of column couplers comes from the 2D effective-medium
solver wm_extraction.py (run_cells.py -> cells.json, run_box.py -> box.json).

Run:  python comparator_model.py      -> comparator_results.csv, equal_volume.json, R_vs_Nel.json
      python decision_table.py         -> decision_table.csv
      python make_figure.py            -> comparator_R_vs_envelope.png
"""
import json, os, sys, itertools, functools
import numpy as np
from scipy.special import erf

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.append(HERE)

# ---------------------------------------------------------------- constants & equal-volume benchmark
hbar = 1.054571817e-34; c = 2.99792458e8; e = 1.602176634e-19
eV4_to_Wm2 = e*(e/hbar)*(e/(hbar*c))**2
T_eV2 = 195.353
GeVcm3_eV4 = 1e9*(hbar*c/e*100)**3
m_a = 1.13137e-3
lam = 2*np.pi*hbar*c/(m_a*e)            # 1.0959 mm
nu = m_a*e/(2*np.pi*hbar)
Qu0 = 473.0; C_SQ = 64/np.pi**4; d_s0 = 3.85e-3

def P_dish_A(g_GeV=1e-12, B_T=10.0, rho=0.4):
    return (g_GeV*1e-9)**2*(B_T*T_eV2)**2*(rho*GeVcm3_eV4)/m_a**2*eV4_to_Wm2

def equal_volume(b, Qs, Qu=Qu0, d=d_s0, C=C_SQ):
    Pd = P_dish_A()
    Ps = b*Pd
    Pp = lambda z: C*z/(1+z)**2*Qu*(2*np.pi*d/lam)*Pd
    R0 = (Ps**2/Qs)/(Pp(2)**2/(Qu/3))
    thr = {R: float(np.sqrt(R/R0)) for R in (1, 2, 10)}
    return dict(P_dish_A=Pd, Ps_A=Ps, Pp_crit_A=Pp(1), Pp_scan_A=Pp(2), R0=R0, thresholds=thr,
                peak_ratio=Ps/Pp(1), peak_ratio_uniformC=Ps/Pp(1)*C)

# ---------------------------------------------------------------- halo coherence (audit Eq. 4, 10)
sig_v = 220e3/np.sqrt(2)/c
ell_a = lam/(2*np.pi*sig_v)              # 0.336 m

def F1(W):
    x = np.maximum(np.asarray(W, float)/ell_a, 1e-9)
    return 2/x**2*(x*np.sqrt(np.pi/2)*erf(x/np.sqrt(2)) + np.exp(-x**2/2) - 1)

def f_box(a, b, cdim=0.0):
    """fraction of the line-integrated signal in one fixed uniform-phase mode for a box a x b x cdim."""
    return float(F1(a)*F1(b)*(F1(cdim) if cdim > 0 else 1.0))

@functools.lru_cache(maxsize=None)
def f_tube_walls(w, L, nwalls, n=24):
    """uniform-mode fraction for the lined walls of one square tube (w x w cross-section, length L)."""
    t = (np.arange(n)+0.5)/n
    z = t*L; u = (t-0.5)*w
    pts = []
    for (ax, sgn) in ([(1, 1), (1, -1)] if nwalls == 2 else [(1, 1), (1, -1), (0, 1), (0, -1)]):
        U, Z = np.meshgrid(u, z)
        P = np.zeros((U.size, 3)); P[:, 2] = Z.ravel()
        P[:, ax] = sgn*w/2; P[:, 1-ax] = U.ravel()
        pts.append(P)
    P = np.vstack(pts)
    D2 = ((P[:, None, :]-P[None, :, :])**2).sum(-1)
    return float(np.exp(-D2/(2*ell_a**2)).mean())

@functools.lru_cache(maxsize=None)
def f_cyl(R, L, nphi=48, nz=40):
    ph = (np.arange(nphi)+0.5)/nphi*2*np.pi; z = (np.arange(nz)+0.5)/nz*L
    PH, Z = np.meshgrid(ph, z)
    P = np.c_[R*np.cos(PH.ravel()), R*np.sin(PH.ravel()), Z.ravel()]
    D2 = ((P[:, None, :]-P[None, :, :])**2).sum(-1)
    return float(np.exp(-D2/(2*ell_a**2)).mean())

# ---------------------------------------------------------------- port extraction tables
cells = [r for r in json.load(open(os.path.join(HERE, "cells.json"))) if r["a_p"] == 1.0 and r["s"] >= 3]
_s = np.array([r["s"] for r in cells]); _e = np.array([r["eta_ext"] for r in cells])
o = np.argsort(_s); _s, _e = _s[o], _e[o]
ETA_BULK_MAX = 1/C_SQ                     # uniform distributed loading of an infinite medium
ell_em0 = lam*np.sqrt(2*Qu0)/(2*np.pi)    # 5.36 mm

def eta_cell(s_m, fR=1.0):
    """bulk lattice extraction efficiency (relative to audit Eq. 1, zeta=2) for column spacing s [m];
    rescaled with the transport length ell_em ~ sqrt(Qu) when the surface resistance changes by fR."""
    s_mm = s_m*1e3*np.sqrt(fR)            # s/ell_em invariant mapping
    if s_mm <= _s[0]:
        return float(_e[0])
    if s_mm >= _s[-1]:
        return float(_e[-1]*(_s[-1]/s_mm)**2)
    return float(np.exp(np.interp(np.log(s_mm), np.log(_s), np.log(_e))))

box = json.load(open(os.path.join(HERE, "box.json")))
_bu = {r["W"]: r["eta_ext"] for r in box if r["kind"] == "uniform"}

def edge_factor(Wmin_m):
    """finite-box (metal side walls) penalty: uniform-loading box / bulk limit, -> 1 for large boxes."""
    Wmm = Wmin_m*1e3
    Ws = np.array(sorted(_bu)); ev = np.array([_bu[w] for w in Ws])/ETA_BULK_MAX
    if Wmm >= Ws[-1]:
        return float(1-(1-ev[-1])*Ws[-1]/Wmm)
    return float(np.interp(Wmm, Ws, ev))

# ---------------------------------------------------------------- loss scenarios (port side)
# alpha: cryogenic WR-3-class line loss at 274 GHz [dB/m]; Lj: excess loss per binary junction [dB];
# eta_coup: antenna/slot-to-line coupling efficiency; sig_phi: rms phase error [rad];
# exit: output run to leave the field region, as a fraction of the envelope axial length; fpack: active fill.
SCEN = {
    "port-optimistic":  dict(alpha=4.0,  Lj=0.10, eta_coup=0.95, sig_phi=0.10, exit=0.00, fpack=1.0,
                             alpha_os=0.1, L_trans=0.1, os_legs=True),
    "port-nominal":     dict(alpha=8.0,  Lj=0.25, eta_coup=0.90, sig_phi=0.30, exit=0.25, fpack=0.8,
                             alpha_os=0.5, L_trans=0.3, os_legs=False),
    "port-pessimistic": dict(alpha=20.0, Lj=0.50, eta_coup=0.80, sig_phi=0.50, exit=0.50, fpack=0.6,
                             alpha_os=2.0, L_trans=0.5, os_legs=False),
}
# oversized (overmoded) low-loss line: used for the exit run after a group combiner (all scenarios) and,
# in the optimistic scenario, for combiner legs longer than OS_LEG_MIN (with one extra transition).
OS_LEG_MIN = 0.05
# stack-side quasi-optical coherent combining of several tubes/faces into one receiver (beam-splitter tree)
# (one fixed set for all scenarios; 0.15 dB per binary beam-splitter stage, 0.3 rad rms phase error)
QO_LDB = 0.15; QO_SIG = 0.30

def series_eff(alpha_dBm, ell):
    """coherent series-fed (slotted-line) column coupler of length ell: (mean amplitude transmission)^2."""
    if ell <= 0:
        return 1.0
    k = alpha_dBm/20*np.log(10)            # amplitude Np/m
    return float(((1-np.exp(-k*ell))/(k*ell))**2)

def port_best(env, sc, fR=1.0, Nrx_max=None, Nel_max=None):
    """maximise the port figure V_p^2 * F_p over: column spacing s, used sub-volume (transverse fraction ft,
    wire-length fraction fl), axial segmentation n_z, transverse group size m x m columns, axial grouping k."""
    alpha = sc["alpha"]*fR
    best = None
    for ft in [1.0, 1/4, 1/16, 1/64, 1/256]:
        for fl in [1.0, 1/2, 1/4, 1/8]:
            A_perp = env["A_perp"]*ft*sc["fpack"]; Lw = env["L_wire"]*fl; Lax = env["L_axis"]
            side = np.sqrt(A_perp)
            Wmin = min(env["W_min"], side, Lw)
            V = A_perp*Lw
            for s in [3e-3, 4e-3, 5e-3, 6e-3, 7.5e-3, 10e-3, 15e-3, 20e-3, 30e-3, 40e-3, 60e-3, 80e-3, 0.12, 0.16, 0.25]:
                Ncol = A_perp/s**2
                if Ncol < 1:
                    continue
                ee = eta_cell(s, fR)*edge_factor(Wmin)
                for nz in [1, 2, 4, 8, 16, 32, 64]:
                    lseg = Lw/nz
                    if lseg < s:
                        continue
                    Nel = Ncol*nz
                    if Nel_max is not None and Nel > Nel_max*1.0001:
                        continue
                    ms = sorted(set([mm for mm in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024] if mm*mm <= Ncol]
                                    + [int(np.ceil(np.sqrt(Ncol)))]))
                    for k in [kk for kk in [1, 2, 4, 8, 16, 32, 64] if kk <= nz]:
                        for m in ms:
                            cols = min(m*m, Ncol)
                            g = min(m*s, side)
                            nports = cols*k
                            stages = int(np.ceil(np.log2(nports))) if nports > 1.0001 else 0
                            path = (g if cols > 1.0001 else 0) + (k*lseg/2 if k > 1 else 0)
                            if sc["os_legs"] and path > OS_LEG_MIN:
                                loss_dB = alpha*OS_LEG_MIN + sc["alpha_os"]*(path-OS_LEG_MIN) + sc["L_trans"]
                            else:
                                loss_dB = alpha*path
                            if sc["exit"] > 0:
                                loss_dB += sc["L_trans"] + sc["alpha_os"]*sc["exit"]*Lax
                            eta = (ee*sc["eta_coup"]*series_eff(alpha, lseg)
                                   * 10**(-(stages*sc["Lj"] + loss_dB)/10)*np.exp(-sc["sig_phi"]**2))
                            fc = f_box(g, g, k*lseg)
                            Nrx = float(np.ceil(Ncol*nz/(cols*k) - 1e-9))
                            if Nrx_max is not None and Nrx > Nrx_max:
                                continue
                            F = (eta*fc)**2/Nrx
                            fig = V**2*F
                            if best is None or fig > best["fig"]:
                                best = dict(fig=fig, F=F, V=V, s=s, m=m, nz=nz, k=k, eta=eta, eta_ext=ee, f_coh=fc,
                                            Nrx=Nrx, Nports=Nel, stages=stages, path=path, ft=ft, fl=fl)
    return best

# ---------------------------------------------------------------- radiating-stack configurations
def stack_configs(env, sc, fR=1.0):
    """explicit geometries; each returns effective area A (B^2-weighted), list of (f_coh, Nrx) options."""
    d = d_s0*fR**(-1/3)
    out = []
    kind = env["kind"]
    if env.get("face") is not None:            # end-face tile(s) (dipole-type: B in the face plane)
        a, b, shape = env["face"]
        A1 = a*b if shape == "rect" else np.pi*a**2/4
        for nfaces, lab in [(1, "single face"), (2, "back-to-back faces")]:
            opts = []
            for p in [a/q for q in (1, 2, 3, 4, 6, 8)]:
                if shape == "rect":
                    n_rx = (a/p)*(b/p); fc = f_box(p, p*b/a)
                else:
                    n_rx = max(1, round(A1/p**2)); fc = f_box(p, p) if p < a else f_box(0.886*a, 0.886*a)
                opts.append((fc, nfaces*n_rx, p, 1))
                if nfaces == 2:   # both faces combined coherently into one channel per patch pair
                    opts.append((fc, n_rx, p, 2))
            out.append(dict(config=lab, A=nfaces*A1, d=d, opts=opts, tube_w=np.nan))
    # tube bundles along the optical axis: lined walls whose plane contains B.  Each tube is divided into
    # BREAD-type focusing cells (coaxial reflector, Liu et al. 2111.12103: L_cell = 2 sqrt2 r_cell, here
    # r_cell = w_int/2).  Cells not at the tube exit must send their focus out through a line (exit loss of
    # the same scenario as the port design: transition + oversized line over sc["exit"]*L_axis).
    nwalls = 4 if kind == "solenoid" else 2
    for w in [0.02, 0.03, 0.05, 0.075, 0.10, 0.15, 0.25, 0.375, 0.5]:
        w_int = w-2*d
        if w_int < 0.012:
            continue
        Lc = min(env["L_axis"], np.sqrt(2)*w_int)
        nc = int(np.floor(env["L_axis"]/Lc+1e-9))
        if env["cross"] == "square":
            nt = int(np.floor(env["side"]/w+1e-9))**2
        else:
            R = env["R"]; xs = (np.arange(-int(R/w)-1, int(R/w)+1))*w
            nt = sum(1 for x in xs for y in xs if max(abs(x), abs(x+w))**2+max(abs(y), abs(y+w))**2 <= R**2)
        if nt < 1:
            continue
        A = nt*nc*nwalls*w_int*Lc
        if nc > 1 and sc["exit"] > 0:
            line = 10**(-(sc["L_trans"]+sc["alpha_os"]*sc["exit"]*env["L_axis"])/10)
            eta_line = (1+(nc-1)*line)/nc          # last cell exits in free space
        else:
            eta_line = 1.0
        fc = f_tube_walls(w_int, Lc, nwalls)
        opts = []
        for ng in [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]:
            if ng > nt:
                break
            G = np.sqrt(ng)*w
            fg = fc*(F1(G)/F1(w_int))**2
            opts.append((fg, max(1, nt/ng)*nc, Lc, ng))
        out.append(dict(config=f"{nwalls}-wall tube bundle w={w*100:.1f} cm (n={nt}x{nc} cells)", A=A, d=d, Nel=nt*nc,
                        opts=opts, tube_w=w, eta_line=eta_line,
                        comb=("line" if nc > 1 else "qo"), alpha=sc["alpha"]*fR, Lj=sc["Lj"], sig=sc["sig_phi"]))
    if kind == "solenoid" and env["cross"] == "circle":       # single BREAD barrel lining
        R, L = env["R"], env["L_axis"]
        opts = []
        for kphi, kz in [(1, 1), (2, 1), (4, 2), (8, 2), (8, 4), (16, 4), (16, 6), (24, 8)]:
            pa = 2*np.pi*R/kphi; pz = L/kz
            fc = f_cyl(R, L) if (kphi, kz) == (1, 1) else f_box(pa, pz)
            opts.append((fc, kphi*kz, (pa, pz), 1))
        out.append(dict(config="BREAD barrel lining", A=2*np.pi*R*L, d=d, opts=opts, tube_w=np.nan))
    if kind == "dipole" and env["cross"] == "circle":         # barrel lining with hoop wires (B_par ~ sin phi)
        R, L = env["R"], env["L_axis"]
        opts = []
        for kphi, kz in [(1, 1), (2, 1), (4, 2), (8, 2), (8, 4), (16, 4)]:
            pa = 2*np.pi*R/kphi; pz = L/kz
            fc = f_box(pa, pz)
            opts.append((fc, kphi*kz, (pa, pz), 1))
        out.append(dict(config="dipole barrel lining (hoop wires)", A=np.pi*R*L, d=d, opts=opts, tube_w=np.nan))
    return out

def stack_best(cfg, eta_s, Nrx_max=None, Nel_max=None):
    best = None
    if Nel_max is not None and cfg.get("Nel", 1) > Nel_max:
        return None
    for fc, nrx, p, ncomb in cfg["opts"]:
        if Nrx_max is not None and nrx > Nrx_max:
            continue
        if ncomb <= 1:
            eqo = 1.0
        elif cfg.get("comb", "qo") == "qo":     # free-space beam-splitter tree at the tube exits
            eqo = 10**(-np.ceil(np.log2(ncomb))*QO_LDB/10)*np.exp(-QO_SIG**2)
        else:                                   # in-field focal horns -> line combiner (port loss model)
            G = np.sqrt(ncomb)*cfg["tube_w"]
            eqo = 10**(-(np.ceil(np.log2(ncomb))*cfg["Lj"] + cfg["alpha"]*G)/10)*np.exp(-cfg["sig"]**2)
        F = (eta_s*cfg.get('eta_line', 1.0)*eqo*fc)**2/nrx
        if best is None or F > best[0]:
            best = (F, fc, nrx, ncomb)
    return best

# ---------------------------------------------------------------- envelopes
def envelopes():
    R_b = np.sqrt(10/(4*np.sqrt(2)*np.pi)); L_b = 2*np.sqrt(2)*R_b      # BREAD: A=2piRL=10 m^2, L=2sqrt2 R
    D_m = 1.25; L_m = 1.30                                                # MADMAX usable disc diameter, field length
    return [
        dict(name="10 cm cube (solenoid-type)", kind="solenoid", cross="square", side=0.10, V=1e-3,
             A_perp=1e-2, L_wire=0.10, L_axis=0.10, W_min=0.10, face=(0.10, 0.10, "rect"), B=10.0),
        dict(name="10 cm cube (dipole-type)", kind="dipole", cross="square", side=0.10, V=1e-3,
             A_perp=1e-2, L_wire=0.10, L_axis=0.10, W_min=0.10, face=(0.10, 0.10, "rect"), B=10.0),
        dict(name="MADMAX-like dipole (9 T, D=1.25 m, L=1.3 m)", kind="dipole", cross="circle", R=D_m/2,
             V=np.pi*D_m**2/4*L_m, A_perp=D_m*L_m*np.pi/4, L_wire=D_m, L_axis=L_m, W_min=min(D_m, L_m),
             face=(D_m, D_m, "disk"), B=9.0),
        dict(name="BREAD-like solenoid (10 T, R=0.75 m, L=2.12 m)", kind="solenoid", cross="circle", R=R_b,
             V=np.pi*R_b**2*L_b, A_perp=np.pi*R_b**2, L_wire=L_b, L_axis=L_b, W_min=2*R_b, face=None, B=10.0),
    ]

# receiver-channel and coupling-element budgets applied to BOTH architectures
# (elements = port couplers for the resonator; in-field focal cells/horns for the stack)
BUDGETS = [("unconstrained", None, None), ("Nrx<=100", 100, None), ("Nrx<=100,Nel<=1e4", 100, 1e4),
           ("Nrx<=100,Nel<=1e3", 100, 1e3), ("Nrx<=10,Nel<=1e2", 10, 1e2), ("Nrx=1", 1, None)]

# ---------------------------------------------------------------- main
if __name__ == "__main__":
    import pandas as pd
    ev_ours = equal_volume(2066, 313); ev_aud = equal_volume(2083, 314); ev_fw = equal_volume(2125, 312)
    R0 = ev_ours["R0"]
    json.dump(dict(audit=ev_aud, ours_sheet=ev_ours, ours_fullwave=ev_fw, ell_em_mm=ell_em0*1e3,
                   single_mode_size_mm=lam*np.sqrt(Qu0)*1e3, ell_a_m=ell_a),
              open(os.path.join(HERE, "equal_volume.json"), "w"), indent=1, default=float)
    rows = []
    eta_s_grid = [0.3, 0.5, 0.6, 0.68, 0.75, 0.9]   # 0.68 = lead finite-tile computed point (range 0.6-0.75)
    for env in envelopes():
        for fR in [0.5, 1.0, 2.0]:
            for blab, Nrx_max, Nel_max in BUDGETS:
                for scn, sc in SCEN.items():
                    cfgs = stack_configs(env, sc, fR)
                    pb = port_best(env, sc, fR, Nrx_max, Nel_max)
                    if pb is None:
                        continue
                    for cfg in cfgs:
                        Vs = cfg["A"]*cfg["d"]
                        for eta_s in eta_s_grid:
                            sb = stack_best(cfg, eta_s, Nrx_max, Nel_max)
                            if sb is None:
                                continue
                            R = R0*(Vs/pb["V"])**2*sb[0]/pb["F"]
                            rows.append(dict(envelope=env["name"], fR_surface_resistance=fR,
                                             budget=blab, Nrx_max=Nrx_max, Nelem_max=Nel_max,
                                             stack_Nelem=cfg.get("Nel", 1),
                                             port_scenario=scn, stack_config=cfg["config"], eta_s=eta_s,
                                             R=R, R0=R0, V_env_m3=env["V"], V_s_m3=Vs, V_p_m3=pb["V"],
                                             stack_fill=Vs/env["V"], stack_fcoh=sb[1], stack_Nrx=sb[2], stack_ncomb=sb[3],
                                             port_eta_p=pb["eta"], port_eta_ext=pb["eta_ext"],
                                             port_fcoh=pb["f_coh"], port_Nrx=pb["Nrx"], port_Nports=pb["Nports"],
                                             port_s_mm=pb["s"]*1e3, port_group_cols=pb["m"]**2, port_nz=pb["nz"],
                                             port_k=pb["k"], port_stages=pb["stages"], port_path_m=pb["path"],
                                             port_frac_transverse=pb["ft"], port_frac_wire=pb["fl"],
                                             alpha_dB_per_m=sc["alpha"]*fR, Lj_dB=sc["Lj"], eta_coup=sc["eta_coup"],
                                             sig_phi=sc["sig_phi"], fpack=sc["fpack"]))
    df = pd.DataFrame(rows)
    df.to_csv(os.path.join(HERE, "comparator_results.csv"), index=False)
    print(len(df), "rows written")
