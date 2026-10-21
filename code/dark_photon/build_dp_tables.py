"""Reproduce dp_projections.csv (dark-photon DM literature inputs, 0.1-20 meV window).

Inputs (fetched beforehand, see fetch commands in dp_literature.md, Sec. 0):
  dp/pdf/<arXiv>.pdf          primary papers (arxiv.org/pdf/<id>)
  dp/al/lim/*.txt             AxionLimits limit_data/DarkPhoton/*.txt   (cajohare/AxionLimits@master)
  dp/al/proj/*.txt            AxionLimits limit_data/DarkPhoton/Projections/*.txt
Digitisations are done from the *vector* drawing paths of the PDF figures (PyMuPDF),
calibrated on the axis tick-line coordinates; Bezier control points are included, so
expect ~0.03-0.05 dex scatter.
"""
import glob, os, re
import numpy as np, pandas as pd
import pymupdf as fitz

MLO, MHI = 1e-4, 2e-2           # eV window for the CSV
rows = []
def add(exp, scen, m, chi, src, notes):
    rows.append(dict(experiment=exp, scenario_phase=scen, mass_eV=m, chi=chi, source=src, notes=notes))

def pts(g):
    P = []
    for it in g['items']:
        for q in it[1:]:
            if hasattr(q, 'x'): P.append((q.x, q.y))
    return np.array(P)

# ---------------------------------------------------------------- 1. AxionLimits files
REPO = 'https://raw.githubusercontent.com/cajohare/AxionLimits/master/limit_data/DarkPhoton/'
for d, sub in [('lim', ''), ('proj', 'Projections/')]:
    for f in sorted(glob.glob(f'dp/al/{d}/*.txt')):
        name = os.path.basename(f)
        raw = open(f, errors='replace').read().splitlines()
        hdr = ' | '.join(l.strip() for l in raw if l.strip().startswith('#')) or '(no header lines in file)'
        dat = []
        for l in raw:
            s = l.strip()
            if not s or s[0] in '#%': continue
            try:
                v = [float(x) for x in s.replace(',', ' ').split()]
                if len(v) >= 2: dat.append(v[:2])
            except ValueError: pass
        if not dat: continue
        a = np.array(dat); sel = a[(a[:, 0] >= MLO) & (a[:, 0] <= MHI)]
        if not len(sel): continue
        note = f'AxionLimits {sub}{name}; raw file rows restricted to {MLO:g}-{MHI:g} eV; chi=1 rows are plotting-box edges'
        if name == 'Karshenboim.txt':
            # NOT (eV, chi): digitised Karshenboim 2010 atomic-spectroscopy bound, columns (lambda [MeV], alpha'),
            # see dp_literature.md Sec. 5. Recorded once as a flag, no numeric rows.
            add('Karshenboim (AxionLimits file)', 'DO NOT USE as DP-DM limit', np.nan, np.nan,
                REPO + sub + name + ' :: ' + hdr,
                'File has no header; columns are (lambda [MeV], alpha_prime) digitised from Karshenboim arXiv:1005.4859/1005.4872 Fig.1 '
                '(hydrogen spectroscopy, fifth-force, 1 sigma); unused in current AxionLimits plotting code; O\'Hare notebook (commit 60c13cf, 2021) '
                'converted it as m=x*1e6/(2pi) eV, chi=10.905*sqrt(y) -> m~0.17-230 eV, chi~6e-6-1e-4. Not a dark-matter bound.')
            continue
        kind = 'projection' if d == 'proj' else 'limit'
        for m, c in sel:
            add(name[:-4], f'AxionLimits {kind}', m, c, REPO + sub + name + ' :: ' + hdr, note)

# ---------------------------------------------------------------- 2. SQWARE Fig. 6 (arXiv:2509.14320v2, p.9)
p = fitz.open('dp/pdf/2509.14320v2.pdf')[8]; dr = p.get_drawings()
X = lambda x: 1e-3 * 10 ** ((x - 401.2) / ((465.8 - 401.2) / np.log10(2)))      # ticks: 1 meV @401.2, 2 meV @465.8
Y = lambda y: 10 ** (-16 + (266.8 - y) / ((266.8 - 144.3) / 4))                 # 1e-16 @266.8, 1e-12 @144.3
sq_col = {(0.25, 0.32, 0.45): 'Config. 3 (LM-S), 300-d sweep, t_obs=7 min/step',
          (0.36, 0.52, 0.71): 'Config. 2 (LM), 300-d sweep, t_obs=34 min/step',
          (0.95, 0.95, 0.35): 'Config. 1 (HM), 300-d sweep, t_obs=84 min/step'}
for g in dr:
    c = g.get('color')
    if g['type'] != 's' or not c or (g.get('width') or 0) < 1.4: continue
    key = tuple(round(v, 2) for v in c)
    if key not in sq_col: continue
    P = pts(g); P = P[np.argsort(P[:, 0])]
    for x, y in P:
        add('SQWARE (Mehrani et al.)', sq_col[key], X(x), Y(y), 'arXiv:2509.14320v2 Fig. 6 (p.9), Eq. (B2)',
            'digitised from vector path (incl. Bezier control points); random polarisation a0->a0/sqrt3; beta includes sin(theta)/sqrt2 (Eq. 3)')
add('SQWARE (Mehrani et al.)', 'Config. 2, 30-d single resonance', X(401.9), Y(241.4), 'arXiv:2509.14320v2 Fig. 6 (p.9)', 'digitised marker centre')
add('SQWARE (Mehrani et al.)', 'Config. 1, 30-d single resonance', X(517.8), Y(217.4), 'arXiv:2509.14320v2 Fig. 6 (p.9)', 'digitised marker centre')

# ---------------------------------------------------------------- 3. BREAD Fig. 2 left (arXiv:2111.12103v2, p.5)
p = fitz.open('dp/pdf/2111.12103.pdf')[4]; d5 = p.get_drawings()
X5 = lambda x: 1e-6 * 10 ** ((x - 91.4) / 31.4)      # eV ; 0.001 meV @91.4, 31.4 pt/decade
Y5 = lambda y: 10 ** (-6 - (y - 82.6) / 12.02)       # 1e-6 @82.6, 12.02 pt/decade
lab = {(0.32, 0.32, 0.32): 'IR Labs bolometer NEP 5e-14', (0.33, 0.15, 0.56): 'KID/TES NEP 2e-19',
       (0.99, 0.31, 0.16): 'QCDet DCR 4 Hz', (0.81, 0.07, 0.34): 'SNSPD DCR 1e-4 Hz', (0.0, 0.43, 0.17): 'Gentec NEP 1e-8'}
for g in d5:
    c = g.get('color')
    if g['type'] != 's' or not c or g['rect'].x0 > 300 or g['rect'].y0 > 205: continue
    key = tuple(round(v, 2) for v in c)
    if key not in lab or len(g['items']) != 1: continue
    w = round(g.get('width') or 0, 2); ds = str(g.get('dashes'))
    sty = '1 day, A=0.7 m2' if w < 0.5 else ('1000 days, NEP/100' if '3.489' in ds else ('1000 days' if '6.03' in ds else '10 days'))
    P = pts(g)
    m = np.logspace(np.log10(X5(P[:, 0].min())), np.log10(X5(P[:, 0].max())), 25)
    k = 10 ** np.interp(np.log10(m), np.log10(X5(P[:, 0])), np.log10(Y5(P[:, 1])))
    for mm, kk in zip(m, k):
        if MLO <= mm <= MHI:
            add('BREAD (Liu et al.)', f'{lab[key]}; {sty}', mm, kk, 'arXiv:2111.12103 Fig. 2 left (p.5), Eq. (2), Table I',
                'digitised; A=10 m2 (unless 0.7 m2), eps_s=0.5, SNR(Z)=5, rho=0.45, alpha_pol=sqrt(2/3) nominal')

# ---------------------------------------------------------------- 4. ALPHA Fig. 17 (arXiv:2210.00017, p.21)
p = fitz.open('dp/pdf/2210.00017.pdf')[20]
g = [g for g in p.get_drawings() if g.get('fill') and tuple(round(v, 2) for v in g['fill']) == (0.75, 0.37, 0.68)][0]
P = pts(g); mX = 10 ** (-6 + (P[:, 0] - 169.9) / 105.4); ch = 10 ** (-16 + (260.7 - P[:, 1]) / 43.93)
for lo, hi in zip(np.logspace(np.log10(2.07e-5), np.log10(1.86e-4), 11)[:-1], np.logspace(np.log10(2.07e-5), np.log10(1.86e-4), 11)[1:]):
    s = (mX >= lo) & (mX < hi)
    if s.any() and np.sqrt(lo * hi) >= MLO:
        add('ALPHA (Millar et al.)', 'Stage I+II discovery potential, fixed polarisation', np.sqrt(lo * hi), ch[s].min(),
            'arXiv:2210.00017 Fig. 17 (p.21)', 'digitised lower edge of red region; 5-45 GHz; outside 0.25-8 meV')

# ---------------------------------------------------------------- 5. published point/band limits read from text
lit = [
 ('One-electron quantum cyclotron (Fan et al.)', 'measured limit, 90% CL, rho=0.3, <sin^2>=2/3', 0.612276e-3, 3.2e-11, 'arXiv:2208.06519v3 Eq. (8), p.4; PRL 129, 261801 (2022)', 'narrow line (33 kHz detector bandwidth); AxionLimits file quotes 4.0e-11'),
 ('Tokyo-Knirck plane-parabolic mirror', 'measured limit, 90% CL, rho=0.45, alpha=sqrt(2/3)', 0.67e-3, 1e-8, 'arXiv:1806.05120v3 abstract & Sec. 4 (p.7); JCAP11(2018)031', 'band 0.67-0.92 meV, "chi >~ 1e-8 excluded"'),
 ('Tokyo-Knirck plane-parabolic mirror', 'measured limit, 90% CL, rho=0.45, alpha=sqrt(2/3)', 0.92e-3, 1e-8, 'arXiv:1806.05120v3 abstract & Sec. 4 (p.7)', 'band edge'),
 ('Dandelion first results (Ourahou et al.)', 'measured limit, 95% CL, rho=0.44, factor 2/3', 0.6e-3, 8.7e-10, 'arXiv:2602.18218v2 Sec. 4, Eq. (8), Fig. 8 (p.10)', '150-350 GHz, 24.7 h; flat value quoted'),
 ('Dandelion first results (Ourahou et al.)', 'measured limit, 95% CL, rho=0.44, factor 2/3', 1.4e-3, 8.7e-10, 'arXiv:2602.18218v2 Sec. 4 (p.10)', 'band edge'),
 ('mm-wave dielectric haloscope (Wei et al., USTC)', 'measured limit, 90% CL, rho=0.45, <cos^2>=1/3', 387.72e-6, 2.2e-11, 'arXiv:2503.23354v2 p.5, Fig. 4', 'band 387.72-391.03 ueV (93.7-94.5 GHz)'),
 ('mm-wave dielectric haloscope (Wei et al., USTC)', 'measured limit, 90% CL, best point', 387.79e-6, 6.68e-12, 'arXiv:2503.23354v2 p.5', 'best point'),
 ('mm-wave dielectric haloscope (Wei et al., USTC)', 'measured limit, 90% CL, rho=0.45, <cos^2>=1/3', 391.03e-6, 2.2e-11, 'arXiv:2503.23354v2 p.5', 'band edge'),
 ('Dark E-field Radio (Godfrey et al.)', 'Phase II projection (schematic Fig. 13)', 3e-3, 7.5e-10, 'AxionLimits Projections/DarkEfield.txt (cites 2101.02805)', 'text of 2101.02805 says Phase II extends to 20 GHz (83 ueV); THz extension only in schematic Fig. 13 dashed line'),
 ('Highly excited electron cyclotron (Fan et al.)', 'projection', np.nan, 2e-16, 'arXiv:2410.05549v3 abstract; PRD 111, 075022 (2025)', '0.1-2.3 meV (25-560 GHz); "down to eps~2e-16"; no curve digitised'),
 ('Plasma haloscope DP (Gelmini et al.)', 'projection, quantum-limited, 5 yr, S/N=3, <cos^2>=1/3', np.nan, 7e-16, 'arXiv:2006.06836 Sec. VI, Eq. (36)-(37), Fig. 2', '6-400 ueV; 0.8 m3, Q=100; AxionLimits ALPHA.txt is a box 20-200 ueV at 1.31e-15'),
]
for r in lit: add(*r)

df = pd.DataFrame(rows)
df.to_csv('dp_projections.csv', index=False)
print(df.groupby('experiment').size().to_string())
