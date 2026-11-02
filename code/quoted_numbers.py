"""Reach numbers quoted in the Letter (abstract, Sensitivity, Eq. (5), End Matter B), read off the Fig. 3/4 curves.
Uses the same benchmark scaling as fig_madgrid*.py: 1 m^2 (x 10^(1/4)), counters D = 1 mHz per channel
(x (1/4)^(1/4) below 4 meV), KIDs at NEP = KID_NEP (x KIDF).  Usage (from code/): python quoted_numbers.py
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, pandas as pd
from fig_madgrid import F, cfac, KIDF, KID_NEP, D, L
A = pd.read_csv(os.path.join(D, "reach_0p25_8_fw.csv")); P = pd.read_csv(os.path.join(D, "dp_reach.csv"))
m = A.m.values; het = F * A.lin10.values; kid = F * KIDF * A.kid10.values; cnt = F * cfac(m) * A.target.values
c = np.loadtxt(os.path.join(L, "CAST.txt")); gcast = np.interp(m * 1e-3, c[:, 0], c[:, 1]); ggc = 4.7e-11
hm = m < 4 - 1e-9          # heterodyne: 2-4 meV cartridge and below (the m = 4.000 row belongs to 4-8 meV)
def upto(y, ref, mask=True): return m[(y < ref) & mask].max()
ks, df = 2.04e-13 * 1.92 * m, 2.04e-13 * 0.75 * m
print(f"KID NEP = {KID_NEP:.1e} W/sqrt(Hz) per channel")
print(f"heterodyne below GC / CAST up to {upto(het, ggc, hm):.2f} / {upto(het, gcast, hm):.2f} meV")
print(f"KIDs       below GC / CAST up to {upto(kid, ggc):.2f} / {upto(kid, gcast):.2f} meV")
print(f"counters below DFSZ-I over {m[cnt < df].min():.2f}-{m[cnt < df].max():.2f} meV, below KSVZ over {m[cnt < ks].min():.2f}-{m[cnt < ks].max():.2f} meV")
print(f"counters / DFSZ-I and / KSVZ at {m[0]:.2f} meV: {cnt[0]/df[0]:.2f}, {cnt[0]/ks[0]:.2f}")
print(f"Eq. (5): counter reach at 1 meV = {np.interp(1.0, m, cnt):.3e} GeV^-1")
mp = P.m_meV.values; cntp = F * cfac(mp) * P.cnt.values; kidp = F * KIDF * P.kid.values
print(f"dark photons: counters min chi = {cntp.min():.2e} at {mp[np.argmin(cntp)]:.2f} meV (at 1 meV {np.interp(1.0, mp, cntp):.2e}); KIDs at 1 meV = {np.interp(1.0, mp, kidp):.2e}"
      f" ({np.interp(1.0, mp, kidp)/3e-14:.1f} x BREAD's 3e-14)")
