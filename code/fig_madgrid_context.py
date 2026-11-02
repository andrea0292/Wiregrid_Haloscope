"""Fig. 3 in wide context (BREAD-style): MADGRID axion reach on top of existing bounds from 1 ueV to 20 meV.
Bounds/projections from AxionLimits (C. O'Hare, github.com/cajohare/AxionLimits), files in data/limits/axionlimits/ (commit 7d375f4 of the compilation)."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fig_madgrid import *          # style, curve, on_curve, block, cfac, faint, F, colors, D
import matplotlib
AL = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 'limits', 'axionlimits')
def ld(f): a = np.loadtxt(os.path.join(AL, f)); return a[:, 0] * 1e3, a[:, 1]      # eV -> meV
HALO_FILES = ('ADMX', 'CAPP', 'HAYSTAC', 'ORGAN', 'QUAX', 'RADES', 'MADMAX', 'GigaBREAD', 'TASEH', 'CAST-CAPP', 'GrAHal', 'UF', 'RBF', 'SHAFT')
C_BND, C_HALO, C_QCD = '#d9d9d9', '#8c8c8c', '#e9cf6a'

def lower_edge(a, n=300):
    """lower envelope (in log g) of a closed projection polygon, meV"""
    x, y = np.log(a[:, 0] * 1e3), np.log(a[:, 1]); xs = np.linspace(x.min(), x.max(), n); lo = np.full(n, np.inf)
    xx = np.r_[x, x[0]]; yy = np.r_[y, y[0]]
    for i in range(len(x)):
        xa, xb, ya, yb = xx[i], xx[i+1], yy[i], yy[i+1]
        if xa == xb: continue
        sel = (xs >= min(xa, xb)) & (xs <= max(xa, xb))
        lo[sel] = np.minimum(lo[sel], ya + (yb - ya) * (xs[sel] - xa) / (xb - xa))
    ok = np.isfinite(lo); return np.exp(xs[ok]), np.exp(lo[ok])

def block2(ax, x, y, lines, ha='right', name_size=16, size=9):
    ax.text(x, y, r'$\mathbf{MADGRID}$', transform=ax.transAxes, fontsize=name_size, color=NAME_C, ha=ha, va='top')
    for i, l in enumerate(lines):
        ax.text(x, y - 0.075 - 0.05 * i, l, transform=ax.transAxes, fontsize=size, color=NAME_C, ha=ha, va='top')

def fig_axion_context(proj=False):
    A = pd.read_csv(os.path.join(D, "reach_0p25_8_fw.csv")); B = pd.read_csv(os.path.join(D, "reach_0p25_8_fw_Rs2.csv"))
    m = A.m.values; k = cfac(m)
    fig, ax = plt.subplots(figsize=(6.15, 4.6))
    x0, x1 = 1e-2, 20.
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(x0, x1); ax.set_ylim(1e-16, 1e-8)
    ax.set_xticks([1e-2, 1e-1, 1, 10]); ax.set_xticklabels(['0.01', '0.1', '1', '10'])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_xlabel(r"$m_a$ [meV]"); ax.set_ylabel(r"$g_{a\gamma\gamma}$ [GeV$^{-1}$]")
    mm = np.geomspace(x0, x1, 600)
    # QCD band and model lines
    q = 2.04e-13 * mm
    ax.fill_between(mm, 0.2 * q, 6 * q, color=C_QCD, alpha=0.35, lw=0, zorder=0)
    ax.plot(mm, 1.92 * q, color='#a8860f', lw=1.2, zorder=0.5); ax.plot(mm, 0.75 * q, color='#a8860f', lw=1.2, zorder=0.5)
    # existing bounds
    mc, gc = ld('CAST.txt'); gcast = np.interp(mm, mc, gc)
    ax.fill_between(mm, gcast, 1, color=C_BND, lw=0, zorder=1); ax.plot(mm, gcast, color='0.5', lw=1.0, zorder=1.1)
    mh, gh = ld('Haloscopes_Combined_microeV-meV.txt'); s = mh <= x1
    ax.fill_between(mh[s], gh[s], 1, color=C_HALO, lw=0, zorder=1.2)
    # individual published haloscope results (the combined file misses narrow and recent ones, e.g. MADMAX 2024, GigaBREAD 2025)
    import glob
    for f in sorted(glob.glob(os.path.join(AL, '*.txt'))):
        n = os.path.basename(f)
        if n.split('.')[0].endswith('_raw') or not n.startswith(HALO_FILES): continue
        d = np.loadtxt(f); ax.fill_between(d[:, 0] * 1e3, d[:, 1], 1, color=C_HALO, lw=0, zorder=1.2)
    ax.plot(mm, np.full_like(mm, 4.7e-11), color='0.3', ls='--', lw=1.4, zorder=1.3)
    mn, gn = ld('NuSTAR_Sun.txt'); s = (mn >= x0) & (gn <= 4.7e-11); ax.plot(mn[s], gn[s], color='0.2', ls=':', lw=2.0, zorder=1.35)
    if proj:
        for f, lab, xy, rot in [('MADMAX.txt', 'MADMAX', (0.082, 3.0e-15), 0), ('ALPHA.txt', 'ALPHA', (0.0175, 8.5e-15), 0),
                                ('DALI.txt', 'DALI', (0.028, 1.2e-15), 0), ('BREAD.txt', 'BREAD', (0.18, 3.6e-13), 0)]:
            mx, gx = lower_edge(np.loadtxt(os.path.join(AL, 'proj', f))); ax.plot(mx, gx, color='0.35', lw=1.2, ls=(0, (4, 2)), zorder=1.5)
            ax.text(*xy, lab, color='0.3', fontsize=10, rotation=rot, zorder=5, **({'ha': 'right', 'va': 'center'} if lab == 'BREAD' else {}))
    # MADGRID
    hm = m < 4 - 1e-9          # heterodyne: 2-4 meV cartridge and below (the m = 4.000 row belongs to 4-8 meV)
    for y, c, kk in [(A.lin10.values, C_HET, hm), (KIDF * A.kid10.values, C_KID, None), (k * A.target.values, C_CNT, None)]:
        faint(ax, m, y, c, kk)
    curve(ax, m, F * A.lin10.values, F * B.lin10.values, C_HET, mask=hm)
    curve(ax, m, F * KIDF * A.kid10.values, F * KIDF * B.kid10.values, C_KID)
    curve(ax, m, F * k * A.target.values, F * k * B.target.values, C_CNT)
    # labels
    ax.text(0.0125, 1.5e-9, "CAST", color='0.35', fontsize=17)
    ax.text(0.08, 5.1e-11, "globular clusters", color="0.25", fontsize=10, ha='left', va='bottom', zorder=6)
    ax.text(0.031, 1.45e-11, "NuSTAR", color="0.2", fontsize=10, ha='left', va='top')
    ax.text(0.0262, 9.5e-14, "haloscopes", color='0.3', fontsize=11, rotation=90, va='bottom', zorder=6)
    on_curve(ax, mm, 1.92 * q, 9.0, "KSVZ", '#8a6d0b', fac=1.25, size=11)
    on_curve(ax, mm, 0.75 * q, 9.0, "DFSZ", '#8a6d0b', fac=0.42, size=11)
    on_curve(ax, mm, 1.92 * q, 0.028, "QCD axion", '#a8860f', fac=1.9, size=14)
    on_curve(ax, m[hm], F * A.lin10.values[hm], 0.62, "heterodyne", C_HET, fac=2.5, size=12.5)
    on_curve(ax, m, F * KIDF * A.kid10.values, 4.6, "KIDs", C_KID, fac=1.45, size=12.5)
    ax.text(0.975, 0.340, "photon counters", color=C_CNT, fontsize=12.5, ha='right', va='top', transform=ax.transAxes)
    block2(ax, 0.975, 0.255, [r"$A=1\,{\rm m^2}$, $B_0=10\,{\rm T}$", r"1 live yr/octave, SNR$\,=5$",
                            r"counters: $\eta_d=0.5$, $D=1\,{\rm mHz}$"], ha='right')
    hs = [Line2D([], [], color='0.3', lw=3), Line2D([], [], color='0.3', lw=1.1, alpha=0.45)]
    ax.legend(hs, [r"$1\,{\rm m^2}$", r"$10\,{\rm m^2}$"], loc='lower right', bbox_to_anchor=(0.62, 0.012),
              fontsize=10, frameon=False, handlelength=1.8, borderaxespad=0.2, labelspacing=0.25)
    from figures_reach import yticks
    yticks(ax)
    t = ax.secondary_xaxis("top", functions=(lambda x: x * 0.2418, lambda f: f / 0.2418))
    t.set_xticks([0.01, 0.1, 1]); t.set_xticklabels(["0.01", "0.1", "1"]); t.minorticks_off()
    t.set_xlabel(r"$\nu$ [THz]", labelpad=8)
    return fig

if __name__ == "__main__":
    style()
    for name, p in [("fig_axion_1m2_context", False), ("fig_axion_1m2_context_proj", True)]:
        fig = fig_axion_context(p); fig.savefig(os.path.join(OUT, name + ".pdf")); fig.savefig(os.path.join(OUT, name + ".png"), dpi=150); plt.close(fig)
    print("ok")
