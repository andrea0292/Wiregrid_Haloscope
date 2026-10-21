"""Fig. S1 of the Letter (v14): boost of a uniformly driven slab on a mirror, Eq. (A1).

Usage (from package/code):  python figures_slab.py  -> ../paper/fig5.pdf (+ .png)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figures_reach import style, SET1

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "paper")


def beta2(L_over_lambda, eps):
    n = np.sqrt(eps + 0j)
    phi = n * 2 * np.pi * L_over_lambda
    b = (1 - np.cos(phi)) / (eps * np.cos(phi) - 1j * n * np.sin(phi))
    return np.abs(b) ** 2


def main():
    style()
    Qm = 300.0
    uopt = (np.pi / (4 * Qm)) ** (2 / 3)
    L = np.linspace(1e-3, 8, 4000)
    cases = [(1.0, 0.0, r"$\varepsilon=1$ (vacuum gap)", SET1[1], '-'),
             (0.050, 1 / Qm, r"$u=0.050$ (thin)", SET1[4], '-'),
             (uopt, 1 / Qm, rf"$u={uopt:.3f}$ (optimum)", SET1[0], '-'),
             (0.006, 1 / Qm, r"$u=0.006$ (thick)", SET1[3], '--')]
    fig, ax = plt.subplots(figsize=(6.15, 4.6))
    for u, im, lab, c, ls in cases:
        ax.plot(L, beta2(L, u + 1j * im), color=c, ls=ls, lw=2.5, label=lab)
    ax.set_yscale('log'); ax.set_xlim(0, 8); ax.set_ylim(0.5, 2e3)
    ax.set_xlabel(r'$L/\lambda$'); ax.set_ylabel(r'$\beta^2$')
    ax.legend(loc='lower center', bbox_to_anchor=(0.5, 1.01), ncol=2, fontsize=13, handlelength=1.8,
              columnspacing=1.2, title=rf'$Q_{{\rm m}}={Qm:.0f}$', title_fontsize=13)
    for ext in ('pdf', 'png'):
        fig.savefig(os.path.join(OUT, 'fig5.' + ext), dpi=150)
    i = np.argmax(beta2(L, uopt + 1j / Qm))
    print('optimum: beta2_max = %.0f at L/lambda = %.2f' % (beta2(L, uopt + 1j / Qm)[i], L[i]))


if __name__ == '__main__':
    main()
