"""Fig. 4 in wide context (same frame as the axion figure): MADGRID dark-photon reach over existing bounds, 0.01-20 meV.
Bounds/projections from AxionLimits (data/limits/axionlimits/dp) plus three recent measured limits not in the compilation
(USTC mm-wave dielectric haloscope, Dandelion first results, one-electron quantum cyclotron) from data/limits/dp_projections_digitized.csv."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_madgrid_context import *     # style, curve, on_curve, faint, F, cfac, colors, lower_edge, block2, AL, D, L
from figures_reach import yticks
ALD = os.path.join(AL, 'dp')
def ldd(f): a = np.loadtxt(os.path.join(ALD, f)); return a[:, 0] * 1e3, a[:, 1]
DM_FILES = ['DP_Combined_AxionSearchesRescaled.txt', 'DP_Combined_DarkMatterSearches.txt', 'BRASS-p.txt', 'DOSUE-RR.txt', 'DOSUE-RR-2.txt',
            'GigaBREAD.txt', 'MADMAX.txt', 'ORPHEUS.txt', 'QUALIPHIDE.txt', 'SHUKET.txt', 'SQuAD.txt', 'SUPAX.txt',
            'Tokyo-Knirck.txt', 'Tokyo-Tomita.txt']
C_PROJ = '0.35'

def bar(ax, m, c, **kw):
    """fill from c up to 1; single points widened to a visible sliver"""
    m, c = np.atleast_1d(m), np.atleast_1d(c)
    if len(m) == 1 or np.ptp(m) < 1e-3 * m.mean(): m = np.array([m.min() * 0.997, m.max() * 1.003]); c = np.full(2, c.min())
    ax.fill_between(m, c, 1, **kw)

def block_dp(ax, x, y, lines, ha='left', name_size=16, size=9, gap=0.068, step=0.044):
    """compact MADGRID block (tighter line spacing so the last line clears the x axis)"""
    ax.text(x, y, r'$\mathbf{MADGRID}$', transform=ax.transAxes, fontsize=name_size, color=NAME_C, ha=ha, va='top')
    for i, l in enumerate(lines):
        ax.text(x, y - gap - step * i, l, transform=ax.transAxes, fontsize=size, color=NAME_C, ha=ha, va='top')

def fig_dp_context(proj=False):
    P = pd.read_csv(os.path.join(D, "dp_reach.csv")); m = P.m_meV.values; k = cfac(m)
    pj = pd.read_csv(os.path.join(L, "dp_projections_digitized.csv"))
    fig, ax = plt.subplots(figsize=(6.15, 4.6))
    x0, x1 = 1e-2, 20.
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(x0, x1); ax.set_ylim(1e-17, 1e-9)
    ax.set_xticks([1e-2, 1e-1, 1, 10]); ax.set_xticklabels(['0.01', '0.1', '1', '10'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$m_X$ [meV]"); ax.set_ylabel(r"kinetic mixing $\chi$")
    mm = np.geomspace(x0, x1, 600)
    # solar dark photons in XENON1T (An et al. 2020) = DP_Combined_Stellar over 1-20 meV (light gray region)
    ms, cs = ldd('DP_Combined_Stellar.txt'); o = np.argsort(ms); cst = np.exp(np.interp(np.log(mm), np.log(ms[o]), np.log(cs[o])))
    ax.fill_between(mm, cst, 1, color=C_BND, lw=0, zorder=1); ax.plot(mm, cst, color='0.5', lw=1.0, zorder=1.1)
    # dark-photon dark-matter cosmology (dashed)
    mw, cw = ldd('Cosmology_Witte_inhomogeneous.txt'); s = (mw >= x0) & (mw <= x1) & (cw < 1e-3); o = np.argsort(mw[s])
    ax.plot(mw[s][o], cw[s][o], color='0.1', ls='--', lw=1.5, zorder=1.6)
    # direct dark-matter searches (dark gray)
    hk = dict(color=C_HALO, lw=0, zorder=1.2)
    for f in DM_FILES:
        a = np.loadtxt(os.path.join(ALD, f)); s = (a[:, 0] * 1e3 <= x1 * 1.5)
        if s.sum() >= 2: ax.fill_between(a[s, 0] * 1e3, a[s, 1], 1, **hk)
    for ex, ph in [("mm-wave dielectric haloscope (Wei et al., USTC)", "measured limit, 90% CL, rho=0.45, <cos^2>=1/3"),
                   ("Dandelion first results (Ourahou et al.)", "measured limit, 95% CL, rho=0.44, factor 2/3"),
                   ("One-electron quantum cyclotron (Fan et al.)", "measured limit, 90% CL, rho=0.3, <sin^2>=2/3")]:
        r = pj[(pj.experiment == ex) & (pj.scenario_phase == ph)].sort_values('mass_eV'); assert len(r), ex
        bar(ax, r.mass_eV.values * 1e3, r.chi.values, **hk)
    if proj:
        pk = dict(color=C_PROJ, lw=1.2, ls=(0, (4, 2)), zorder=1.5)
        for f, lab, xy in [('MADMAX.txt', 'MADMAX', (0.105, 3.5e-16)), ('ALPHA.txt', 'ALPHA', (0.0125, 2.2e-15))]:
            mx, cx = lower_edge(np.loadtxt(os.path.join(ALD, 'proj', f))); ax.plot(mx, cx, **pk)
            ax.text(*xy, lab, color='0.3', fontsize=10, zorder=5)
        for ph, lab, xy in [("KID/TES NEP 2e-19; 1000 days", "BREAD (KID)", (0.19, 3.1e-14))]:
            b = pj[(pj.experiment == "BREAD (Liu et al.)") & (pj.scenario_phase == ph)].sort_values("mass_eV")
            ax.plot(b.mass_eV * 1e3, b.chi, **pk); ax.text(*xy, lab, color='0.3', fontsize=10, zorder=5, ha='right', va='center')
        s = pj[(pj.experiment == "SQWARE (Mehrani et al.)") & (pj.scenario_phase == "Config. 3 (LM-S), 300-d sweep, t_obs=7 min/step")].sort_values("mass_eV")
        ax.plot(s.mass_eV * 1e3, s.chi, color=C_PROJ, lw=1.4, ls=':', zorder=1.5); ax.text(0.42, 1.6e-15, "SQWARE", color='0.3', fontsize=10, zorder=5)
    # MADGRID
    hm = np.isfinite(P.het.values) & (m < 4 - 1e-9)
    for y, c, kk in [(P.het.values, C_HET, hm), (KIDF * P.kid.values, C_KID, None), (k * P.cnt.values, C_CNT, None)]:
        faint(ax, m, y, c, kk)
    curve(ax, m, F * P.het.values, F * P.het_Rs2.values, C_HET, mask=hm)
    curve(ax, m, F * KIDF * P.kid.values, F * KIDF * P.kid_Rs2.values, C_KID)
    curve(ax, m, F * k * P.cnt.values, F * k * P.cnt_Rs2.values, C_CNT, dashed_below=0.8)
    # labels
    ax.text(18., 4.3e-10, "XENON1T", color='0.35', fontsize=13, ha='right')
    ax.text(0.13, 5.5e-10, "cosmology", color='0.25', fontsize=10, ha='left', va='center', zorder=6)
    ax.text(0.0365, 1.1e-13, "DM searches", color='0.3', fontsize=11, rotation=90, va='bottom', zorder=6)
    on_curve(ax, m[hm], F * P.het.values[hm], 0.62, "heterodyne", C_HET, fac=2.2, size=12.5)
    on_curve(ax, m, F * KIDF * P.kid.values, 4.6, "KIDs", C_KID, fac=1.5, size=12.5)
    ax.text(0.975, 0.122, "photon counters", color=C_CNT, fontsize=12.5, ha='right', va='bottom', transform=ax.transAxes)
    block_dp(ax, 0.025, 0.218, [r"$A=1\,{\rm m^2}$, no magnet", r"1 live yr/octave, SNR$\,=5$",
                              r"counters: $\eta_d=0.5$, $D=1\,{\rm mHz}$"], ha='left')
    hs = [Line2D([], [], color='0.3', lw=3), Line2D([], [], color='0.3', lw=1.1, alpha=0.45)]
    ax.legend(hs, [r"$1\,{\rm m^2}$", r"$10\,{\rm m^2}$"], loc='lower left', bbox_to_anchor=(0.40, 0.0),
              fontsize=10, frameon=False, handlelength=1.8, borderaxespad=0.2, labelspacing=0.25)
    yticks(ax)
    t = ax.secondary_xaxis("top", functions=(lambda x: x * 0.2418, lambda f: f / 0.2418))
    t.set_xticks([0.01, 0.1, 1]); t.set_xticklabels(["0.01", "0.1", "1"]); t.minorticks_off()
    t.set_xlabel(r"$\nu$ [THz]", labelpad=8)
    return fig

if __name__ == "__main__":
    style()
    for name, p in [("fig_dp_1m2_context", False), ("fig_dp_1m2_context_proj", True)]:
        fig = fig_dp_context(p); fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); plt.close(fig)
    print("ok")
