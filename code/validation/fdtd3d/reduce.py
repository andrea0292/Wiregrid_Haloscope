"""Reduce raw FDTD records (raw/res_*.npz) to spectra + aperture profile at the peak (../../../data/fdtd3d/*.npz)
and write the summary table ../../../data/fdtd3d/fdtd3d_finite_wire_summary.csv."""
import os, sys, numpy as np, pandas as pd
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze as A
RAW = os.path.join(HERE, 'raw'); OUT = os.path.join(HERE, '..', '..', '..', 'data', 'fdtd3d'); os.makedirs(OUT, exist_ok=True)

def reduce(tag):
    R = dict(np.load(os.path.join(RAW, f'res_{tag}.npz'))); M = A.finite_metrics(R)
    ip = int(np.argmax(M['b2_tot']))
    d = dict(f=M['f'], b2_tot=M['b2_tot'], b2_G=M['b2_G'], b2_mid=M['b2_mid'], y=M['y'],
             prof=np.abs(M['E'][ip])**2, W=M['W'], nu_c=float(R['nu_c']), a=float(R['a']), N=int(R['N']))
    np.savez_compressed(os.path.join(OUT, f'{tag}.npz'), **d); return d

if __name__ == '__main__':
    rows = []
    for band in ['1-2', '0.25-0.5']:
        P = reduce(f'{band}_per'); nuc = P['nu_c']
        fp, bp_mid, qp = A.peak_q(P['f'], P['b2_mid']); _, bp_tot, _ = A.peak_q(P['f'], P['b2_tot'])
        Ip_mid = A.scanI(P['f'], P['b2_mid']); Ip_tot = A.scanI(P['f'], P['b2_tot'])
        rows.append(dict(band=band, ends='infinite (periodic)', W_over_lambda=np.nan, f_peak_over_nuc=fp/nuc,
                         beta2_peak=bp_mid, Q_r=qp, tot_ratio=1, mid_ratio=1, G_ratio=np.nan, sqrtI_tot=1, sqrtI_G=np.nan))
        for fr in ['bars', 'walls', 'open']:
            M = reduce(f'{band}_{fr}_W100mm')
            ft, bt, qt = A.peak_q(M['f'], M['b2_tot']); fm, bm, qm = A.peak_q(M['f'], M['b2_mid']); fg, bg, qg = A.peak_q(M['f'], M['b2_G'])
            rows.append(dict(band=band, ends=fr, W_over_lambda=0.10*nuc/2.99792458e8, f_peak_over_nuc=fm/nuc, beta2_peak=bm,
                             Q_r=qt, tot_ratio=bt/bp_tot, mid_ratio=bm/bp_mid, G_ratio=bg/bp_mid,
                             sqrtI_tot=np.sqrt(A.scanI(M['f'], M['b2_tot'])/Ip_tot),
                             sqrtI_G=np.sqrt(A.scanI(M['f'], M['b2_G'])/Ip_mid), df_over_nuc=(fm-fp)/nuc, Q_ratio=qt/qp))
    T = pd.DataFrame(rows); T.to_csv(os.path.join(OUT, 'fdtd3d_finite_wire_summary.csv'), index=False)
    print(T.round(4).to_string())
