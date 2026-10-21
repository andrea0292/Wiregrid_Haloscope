"""Exploratory versions of Fig. 1 and Fig. 2 for a 1 m^2 array (all curves).
Same detectors per tile, same per-tile background and exposure (1 live yr per octave): the number of tiles drops
by 10 and every curve rises by 10^(1/4) = 1.78 (g, chi ~ (N_tile T)^(-1/4)). Thin faint lines: 10 m^2 baseline.
Usage: python figures_reach_1m2.py -> ../paper/fig_axion_1m2.(pdf|png), ../paper/fig_dp_1m2.(pdf|png)
"""
import os, sys; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from figures_reach import (style, lim, top_axis, xaxis, curve, D, L, OUT, C_HET, C_KID, C_CNT, GRAY)
F = 10 ** 0.25

def faint(ax, m, y, c, mask=None):
    if mask is None: mask = np.isfinite(y)
    ax.plot(m[mask], y[mask], color=c, lw=1.1, alpha=0.45)

def fig_axion():
    A = pd.read_csv(os.path.join(D, "reach_0p25_8_fw.csv")); B = pd.read_csv(os.path.join(D, "reach_0p25_8_fw_Rs2.csv"))
    m = A.m.values
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_a$ [meV]"); ax.set_yscale("log"); ax.set_ylim(5e-15, 1e-9)
    ax.set_ylabel(r"$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]")
    mc, gc = lim("CAST.txt"); mm = np.geomspace(0.25, 8, 400); gcast = np.interp(mm * 1e-3, mc, gc)
    ax.fill_between(mm, gcast, 1, color=GRAY, lw=0); ax.plot(mm, gcast, color='0.45', lw=1.2)
    ax.axhline(4.7e-11, color='0.3', ls='--', lw=1.8)
    mn, gn = lim("NuSTAR_Sun.txt"); s = (mn > 2.4e-4) & (mn < 8e-3); ax.plot(mn[s] * 1e3, gn[s], color='0.25', ls=':', lw=2.2)
    q = 2.04e-13 * mm; ax.fill_between(mm, 0.75 * q, 1.92 * q, color='#e9cf6a', alpha=0.75, lw=0)
    hm = m <= 4.0001
    for y, c, k in [(A.lin10.values, C_HET, hm), (A.kid10.values, C_KID, None), (A.target.values, C_CNT, None)]:
        faint(ax, m, y, c, k)
    curve(ax, m, F * A.lin10.values, F * B.lin10.values, C_HET, mask=hm)
    curve(ax, m, F * A.kid10.values, F * B.kid10.values, C_KID)
    curve(ax, m, F * A.target.values, F * B.target.values, C_CNT, dashed_below=0.8)
    kw = dict(fontsize=15)
    ax.text(0.27, 2.0e-10, "CAST", color='0.3', fontsize=16)
    ax.text(0.62, 1.25e-10, "NuSTAR", color='0.25', fontsize=12)
    ax.text(0.27, 2.6e-11, "globular clusters", color='0.3', fontsize=12)
    ax.text(0.27, 6.0e-13, "heterodyne", color=C_HET, **kw)
    ax.text(1.2, 6.0e-12, "KID arrays", color=C_KID, **kw)
    ax.text(1.25, 1.6e-14, "photon counters", color=C_CNT, **kw)
    ax.text(0.31, 5.2e-14, "QCD axion", color='#a8860f', fontsize=15, rotation=18)
    ax.text(3.2, 1.15e-11, r"$1\,{\rm m^2}$, 10 T", fontsize=14)
    top_axis(ax); return fig

def fig_dp():
    P = pd.read_csv(os.path.join(D, "dp_reach.csv")); m = P.m_meV.values
    pj = pd.read_csv(os.path.join(L, "dp_projections_digitized.csv"))
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_X$ [meV]"); ax.set_yscale("log"); ax.set_ylim(1e-16, 1e-11)
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
    for y, c, k in [(P.het.values, C_HET, hm), (P.kid.values, C_KID, None), (P.cnt.values, C_CNT, None)]:
        faint(ax, m, y, c, k)
    curve(ax, m, F * P.het.values, F * P.het_Rs2.values, C_HET, mask=hm)
    curve(ax, m, F * P.kid.values, F * P.kid_Rs2.values, C_KID)
    curve(ax, m, F * P.cnt.values, F * P.cnt_Rs2.values, C_CNT, dashed_below=0.8)
    ax.plot(m, P.proto_kid.values, color=C_KID, lw=1.8, ls='-.')
    ax.text(0.52, 1.95e-13, "heterodyne", color=C_HET, fontsize=13)
    ax.text(4.15, 3.1e-13, "KID arrays", color=C_KID, fontsize=12)
    ax.text(0.3, 7.5e-13, "one tile, KID", color=C_KID, fontsize=12)
    ax.text(2.1, 1.9e-16, "photon counters", color=C_CNT, fontsize=14)
    ax.text(2.3, 3.5e-12, r"$1\,{\rm m^2}$, no magnet", fontsize=13)
    top_axis(ax); return fig

if __name__ == "__main__":
    style()
    for name, f in [("fig_axion_1m2", fig_axion), ("fig_dp_1m2", fig_dp)]:
        fig = f(); fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); plt.close(fig)
    print("written")
