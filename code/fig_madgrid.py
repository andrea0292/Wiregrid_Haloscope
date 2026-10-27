"""Fig. 3/4 variants with the conservative counter benchmark eta = 0.5, D = 1 mHz per channel (all tiles).
Counter curves scale exactly as (D/eta^2)^(1/4) (thermal term negligible): x0.707 below 4 meV (10 cm tiles, was 4 mHz) and x1 above 4 meV (unchanged).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
NAME_C = '0.25'
def on_curve(ax, m, y, x0, text, c, fac=1.0, size=14, dx=1.06):
    """label placed along a log-log curve at x0, rotated to its local slope, shifted by factor fac in y"""
    y0 = np.exp(np.interp(np.log(x0), np.log(m), np.log(y))); y1 = np.exp(np.interp(np.log(x0*dx), np.log(m), np.log(y)))
    p0, p1 = ax.transData.transform([(x0, y0), (x0*dx, y1)])
    ang = np.degrees(np.arctan2(p1[1]-p0[1], p1[0]-p0[0]))
    ax.text(x0, y0*fac, text, color=c, fontsize=size, rotation=ang, rotation_mode='anchor', ha='left', va='bottom')
def block(ax, x, y, lines, ha='left', name_size=17, size=11):
    t = ax.text(x, y, r'$\mathbf{MADGRID}$', transform=ax.transAxes, fontsize=name_size, color=NAME_C, ha=ha, va='top')
    for i, l in enumerate(lines):
        ax.text(x, y - 0.085 - 0.062*i, l, transform=ax.transAxes, fontsize=size, color=NAME_C, ha=ha, va='top')
from figures_reach import (style, lim, top_axis, xaxis, curve, D, L, OUT, C_HET, C_KID, C_CNT, GRAY)
F = 10 ** 0.25
KID_NEP = 6e-19                      # W/sqrt(Hz) per channel; the stored tables use 3e-19
KIDF = (KID_NEP / 3e-19) ** 0.5       # g, chi scale as NEP^(1/2) (KIDs enter as a dark rate NEP^2/2(h nu)^2)
def cfac(m): return np.where(m >= 4 - 1e-9, 1.0, 0.25**0.25)
def faint(ax, m, y, c, mask=None):
    if mask is None: mask = np.isfinite(y)
    ax.plot(m[mask], y[mask], color=c, lw=1.1, alpha=0.45)
def old(ax, m, y):
    ax.plot(m[m >= 0.8], y[m >= 0.8], color=C_CNT, lw=1.3, ls=(0, (1.2, 1.6)), alpha=0.9)

def fig_axion():
    A = pd.read_csv(os.path.join(D, "reach_0p25_8_fw.csv")); B = pd.read_csv(os.path.join(D, "reach_0p25_8_fw_Rs2.csv"))
    m = A.m.values; k = cfac(m)
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_a$ [meV]"); ax.set_yscale("log"); ax.set_ylim(5e-15, 1e-9)
    ax.set_ylabel(r"$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]")
    mc, gc = lim("CAST.txt"); mm = np.geomspace(0.25, 8, 400); gcast = np.interp(mm * 1e-3, mc, gc)
    ax.fill_between(mm, gcast, 1, color=GRAY, lw=0); ax.plot(mm, gcast, color='0.45', lw=1.2)
    ax.axhline(4.7e-11, color='0.3', ls='--', lw=1.8)
    mn, gn = lim("NuSTAR_Sun.txt"); s = (mn > 2.4e-4) & (mn < 8e-3); ax.plot(mn[s] * 1e3, gn[s], color='0.25', ls=':', lw=2.2)
    q = 2.04e-13 * mm; ax.fill_between(mm, 0.75 * q, 1.92 * q, color='#e9cf6a', alpha=0.75, lw=0)
    hm = m <= 4.0001
    for y, c, kk in [(A.lin10.values, C_HET, hm), (KIDF * A.kid10.values, C_KID, None), (k * A.target.values, C_CNT, None)]:
        faint(ax, m, y, c, kk)
    curve(ax, m, F * A.lin10.values, F * B.lin10.values, C_HET, mask=hm)
    curve(ax, m, F * KIDF * A.kid10.values, F * KIDF * B.kid10.values, C_KID)
    curve(ax, m, F * k * A.target.values, F * k * B.target.values, C_CNT, dashed_below=0.8)
    kw = dict(fontsize=15)
    ax.text(0.27, 2.0e-10, "CAST", color='0.3', fontsize=16)
    ax.text(0.62, 1.25e-10, "NuSTAR", color='0.25', fontsize=12)
    ax.text(0.27, 2.6e-11, "globular clusters", color='0.3', fontsize=12)
    on_curve(ax, m[hm], F * A.lin10.values[hm], 0.262, "heterodyne", C_HET, fac=0.27, size=15)
    ax.text(1.2, 6.0e-12, "KID arrays", color=C_KID, **kw)
    on_curve(ax, m, F * k * A.target.values, 2.35, "photon counters", C_CNT, fac=2.2, size=15)
    ax.text(0.31, 5.2e-14, "QCD axion", color='#a8860f', fontsize=15, rotation=18)
    block(ax, 0.975, 0.30, [r"$A=1\,{\rm m^2}$, $B_0=10\,{\rm T}$", r"1 live yr/octave, SNR$\,=5$", r"counters: $\eta_d=0.5$, $D=1\,{\rm mHz}$"], ha='right')
    hs = [Line2D([], [], color='0.3', lw=3), Line2D([], [], color='0.3', lw=1.1, alpha=0.45), Patch(color='0.3', alpha=0.22, lw=0)]
    ax.legend(hs, [r"$1\,{\rm m^2}$", r"$10\,{\rm m^2}$", r"$R_s\times2$"], loc='upper left', bbox_to_anchor=(0.27, 0.565), fontsize=10.5,
              frameon=False, handlelength=1.8, borderaxespad=0.4, labelspacing=0.3)
    top_axis(ax); return fig

def fig_dp():
    P = pd.read_csv(os.path.join(D, "dp_reach.csv")); m = P.m_meV.values; k = cfac(m)
    pj = pd.read_csv(os.path.join(L, "dp_projections_digitized.csv"))
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_X$ [meV]"); ax.set_yscale("log"); ax.set_ylim(5e-17, 1e-10)
    ax.set_ylabel(r"kinetic mixing $\chi$")
    b = pj[(pj.experiment == "BREAD (Liu et al.)") & (pj.scenario_phase == "KID/TES NEP 2e-19; 1000 days")].sort_values("mass_eV")
    ax.plot(b.mass_eV * 1e3, b.chi, color='#a1663e', ls=(0, (6, 3)), lw=2.0)
    ax.text(0.27, 1.35e-14, "BREAD (one KID/TES)", color='#a1663e', fontsize=12)
    for cfg, lab, xy in [("Config. 1 (HM), 300-d sweep, t_obs=84 min/step", "SQWARE 1", (3.0, 1.1e-14)),
                         ("Config. 2 (LM), 300-d sweep, t_obs=34 min/step", "SQWARE 2", (2.1, 6.0e-15)),
                         ("Config. 3 (LM-S), 300-d sweep, t_obs=7 min/step", "SQWARE 3", (2.1, 2.2e-15))]:
        s = pj[(pj.experiment == "SQWARE (Mehrani et al.)") & (pj.scenario_phase == cfg)].sort_values("mass_eV")
        ax.plot(s.mass_eV * 1e3, s.chi, color='0.35', ls=':', lw=2.0); ax.text(*xy, lab, color='0.35', fontsize=11)
    hm = np.isfinite(P.het.values)
    for y, c, kk in [(P.het.values, C_HET, hm), (KIDF * P.kid.values, C_KID, None), (k * P.cnt.values, C_CNT, None)]:
        faint(ax, m, y, c, kk)
    curve(ax, m, F * P.het.values, F * P.het_Rs2.values, C_HET, mask=hm)
    curve(ax, m, F * KIDF * P.kid.values, F * KIDF * P.kid_Rs2.values, C_KID)
    curve(ax, m, F * k * P.cnt.values, F * k * P.cnt_Rs2.values, C_CNT, dashed_below=0.8)
    ax.plot(m, KIDF * P.proto_kid.values, color=C_KID, lw=1.8, ls='-.')
    ax.text(0.52, 1.95e-13, "heterodyne", color=C_HET, fontsize=13)
    ax.text(4.15, 3.1e-13, "KID arrays", color=C_KID, fontsize=12)
    ax.text(0.3, 7.5e-13, "one tile, KID", color=C_KID, fontsize=12)
    on_curve(ax, m, F * k * P.cnt.values, 4.25, "photon counters", C_CNT, fac=1.6, size=12.5)
    block(ax, 0.025, 0.975, [r"$A=1\,{\rm m^2}$, no magnet", r"1 live yr/octave, SNR$\,=5$", r"counters: $\eta_d=0.5$, $D=1\,{\rm mHz}$"], ha='left')
    top_axis(ax); return fig

if __name__ == "__main__":
    style()
    for name, f in [("fig_axion_1m2", fig_axion), ("fig_dp_1m2", fig_dp)]:
        fig = f(); fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); plt.close(fig)
    print("written")
