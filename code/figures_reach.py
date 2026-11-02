"""Reach figures of the Letter (v14): Fig. 1 (axion) and Fig. 2 (dark photon).

Usage (from package/code):  python figures_reach.py  -> ../paper/fig_axion.pdf, ../paper/fig_dp.pdf (+ .png)
Inputs: ../data/reach_0p25_8_fw.csv, ../data/reach_0p25_8_fw_Rs2.csv, ../data/dp_reach.csv,
        ../data/limits/*.txt (AxionLimits files, github.com/cajohare/AxionLimits),
        ../data/limits/dp_projections_digitized.csv (digitized published projections).
"""
import os
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, "..", "data")
L = os.path.join(D, "limits")
OUT = os.path.join(HERE, "..", "paper")
SET1 = ['#d34a4c', '#4f86b5', '#5ba453', '#9b5fa6', '#f08a33', '#f2f25a', '#a1663e', '#ee8fc0', '#999999', '#d34a4c']
C_HET, C_KID, C_CNT = SET1[1], SET1[3], SET1[2]
GRAY = '#c9c9c9'


def style():
    import glob
    from matplotlib import font_manager
    for p in glob.glob("/usr/local/texlive/*/texmf-dist/fonts/opentype/public/cm-unicode/cmun*.otf"):
        try: font_manager.fontManager.addfont(p)
        except Exception: pass
    plt.rcParams.update({
        "font.family": "serif", "font.serif": ["CMU Serif", "Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset": "cm", "axes.linewidth": 1.8, "axes.labelsize": 19, "xtick.labelsize": 17.5,
        "ytick.labelsize": 17.5, "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True,
        "ytick.right": True, "xtick.minor.visible": True, "ytick.minor.visible": True,
        "xtick.major.size": 7, "ytick.major.size": 7, "xtick.minor.size": 3, "ytick.minor.size": 3,
        "lines.linewidth": 2.5, "legend.frameon": False, "savefig.bbox": "tight"})


def lim(fn):
    a = np.loadtxt(os.path.join(L, fn), comments="#")
    return a[:, 0], a[:, 1]


def yticks(ax):
    ax.yaxis.set_major_locator(matplotlib.ticker.LogLocator(base=10, numticks=30))
    ax.yaxis.set_minor_locator(matplotlib.ticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=30))
    ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())


def top_axis(ax):
    yticks(ax)
    t = ax.secondary_xaxis("top", functions=(lambda m: m * 0.2418, lambda f: f / 0.2418))
    t.set_xticks([0.1, 0.2, 0.5, 1, 1.9]); t.set_xticklabels(["0.1", "0.2", "0.5", "1", "1.9"])
    t.minorticks_off(); t.set_xlabel(r"$\nu$ [THz]", labelpad=8)


def xaxis(ax, lab):
    ax.set_xscale("log"); ax.set_xlim(0.25, 8)
    ax.set_xticks([0.25, 0.5, 1, 2, 4, 8]); ax.set_xticklabels(["0.25", "0.5", "1", "2", "4", "8"])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(lab)


def curve(ax, m, y, y2, c, mask=None, dashed_below=None):
    if mask is None: mask = np.isfinite(y)
    m, y, y2 = m[mask], y[mask], y2[mask]
    ax.fill_between(m, y, y2, color=c, alpha=0.22, lw=0)
    if dashed_below is None:
        ax.plot(m, y, color=c, lw=3)
    else:
        s = m >= dashed_below
        ax.plot(m[s], y[s], color=c, lw=3)
        ax.plot(m[m <= dashed_below * 1.001], y[m <= dashed_below * 1.001], color=c, lw=3, ls=(0, (5, 2.5)))


def fig_axion():
    A = pd.read_csv(os.path.join(D, "reach_0p25_8_fw.csv"))
    B = pd.read_csv(os.path.join(D, "reach_0p25_8_fw_Rs2.csv"))
    m = A.m.values
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_a$ [meV]"); ax.set_yscale("log"); ax.set_ylim(5e-15, 1e-9)
    ax.set_ylabel(r"$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]")
    mc, gc = lim("CAST.txt"); mm = np.geomspace(0.25, 8, 400)
    gcast = np.interp(mm * 1e-3, mc, gc)
    ax.fill_between(mm, gcast, 1, color=GRAY, lw=0); ax.plot(mm, gcast, color='0.45', lw=1.2)
    ax.axhline(4.7e-11, color='0.3', ls='--', lw=1.8)
    mn, gn = lim("NuSTAR_Sun.txt"); s = (mn > 2.4e-4) & (mn < 8e-3)
    ax.plot(mn[s] * 1e3, gn[s], color='0.25', ls=':', lw=2.2)
    q = 2.04e-13 * mm
    ax.fill_between(mm, 0.75 * q, 1.92 * q, color='#e9cf6a', alpha=0.75, lw=0)
    curve(ax, m, A.lin10.values, B.lin10.values, C_HET, mask=m <= 4.0001)
    curve(ax, m, A.kid10.values, B.kid10.values, C_KID)
    curve(ax, m, A.target.values, B.target.values, C_CNT)
    # same counters (50%, 4 mHz per channel), 1 m^2, same exposure: g scales as A^(-1/4)
    f1 = 10 ** 0.25
    g1, g1b = A.target.values * f1, B.target.values * f1
    ax.fill_between(m, g1, g1b, color=C_CNT, alpha=0.12, lw=0)
    s1 = m >= 0.8
    ax.plot(m[s1], g1[s1], color=C_CNT, lw=1.8, alpha=0.85)
    ax.plot(m[m <= 0.8008], g1[m <= 0.8008], color=C_CNT, lw=1.8, alpha=0.85, ls=(0, (4, 2.5)))
    kw = dict(fontsize=15)
    ax.text(0.27, 2.0e-10, "CAST", color='0.3', fontsize=16)
    ax.text(0.62, 1.25e-10, "NuSTAR", color='0.25', fontsize=12)
    ax.text(0.27, 2.6e-11, "globular clusters", color='0.3', fontsize=12)
    ax.text(0.27, 1.0e-12, "heterodyne", color=C_HET, **kw)
    ax.text(1.2, 6.5e-12, "KID arrays", color=C_KID, **kw)
    ax.text(1.25, 2.6e-14, r"photon counters, $10\,{\rm m^2}$", color=C_CNT, fontsize=14)
    ax.text(2.45, 1.45e-12, r"$1\,{\rm m^2}$", color=C_CNT, fontsize=13)
    ax.text(0.31, 5.2e-14, "QCD axion", color='#a8860f', fontsize=15, rotation=18)
    ax.text(3.0, 1.15e-11, r"$10\,{\rm m^2}$, 10 T", fontsize=14)
    top_axis(ax)
    return fig


def fig_dp():
    P = pd.read_csv(os.path.join(D, "dp_reach.csv"))
    m = P.m_meV.values
    pj = pd.read_csv(os.path.join(L, "dp_projections_digitized.csv"))
    fig, ax = plt.subplots(figsize=(6.15, 3.8))
    xaxis(ax, r"$m_X$ [meV]"); ax.set_yscale("log"); ax.set_ylim(1e-16, 1e-11)
    ax.set_ylabel(r"kinetic mixing $\chi$")
    # one published BREAD projection (KID/TES, 1000 days) and SQWARE configurations
    b = pj[(pj.experiment == "BREAD (Liu et al.)") & (pj.scenario_phase == "KID/TES NEP 2e-19; 1000 days")].sort_values("mass_eV")
    ax.plot(b.mass_eV * 1e3, b.chi, color='#a1663e', ls=(0, (6, 3)), lw=2.0)
    ax.text(0.27, 1.35e-14, "BREAD (one KID/TES)", color='#a1663e', fontsize=12)
    for cfg, lab, xy in [("Config. 1 (HM), 300-d sweep, t_obs=84 min/step", "SQWARE 1", (3.0, 1.1e-14)),
                         ("Config. 2 (LM), 300-d sweep, t_obs=34 min/step", "SQWARE 2", (2.1, 6.0e-15)),
                         ("Config. 3 (LM-S), 300-d sweep, t_obs=7 min/step", "SQWARE 3", (2.1, 2.2e-15))]:
        s = pj[(pj.experiment == "SQWARE (Mehrani et al.)") & (pj.scenario_phase == cfg)].sort_values("mass_eV")
        ax.plot(s.mass_eV * 1e3, s.chi, color='0.35', ls=':', lw=2.0)
        ax.text(*xy, lab, color='0.35', fontsize=11)
    curve(ax, m, P.het.values, P.het_Rs2.values, C_HET, mask=np.isfinite(P.het.values))
    curve(ax, m, P.kid.values, P.kid_Rs2.values, C_KID)
    curve(ax, m, P.cnt.values, P.cnt_Rs2.values, C_CNT)
    ax.plot(m, P.proto_kid.values, color=C_KID, lw=1.8, ls='-.')
    kw = dict(fontsize=15)
    ax.text(0.6, 1.65e-13, "heterodyne", color=C_HET, fontsize=14)
    ax.annotate("KID arrays", xy=(0.36, 4.6e-14), xytext=(0.27, 1.25e-13), color=C_KID, fontsize=14,
                arrowprops=dict(arrowstyle="-", color=C_KID, lw=1.2, shrinkA=0, shrinkB=2))
    ax.text(0.3, 7.5e-13, "one tile, KID", color=C_KID, fontsize=12)
    ax.text(2.1, 1.9e-16, "photon counters", color=C_CNT, fontsize=14)
    ax.text(2.2, 3.5e-12, r"$10\,{\rm m^2}$, no magnet", fontsize=13)
    top_axis(ax)
    return fig


if __name__ == "__main__":
    style()
    for name, f in [("fig_axion", fig_axion), ("fig_dp", fig_dp)]:
        fig = f()
        fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=150)
        plt.close(fig)
    print("written", OUT)
