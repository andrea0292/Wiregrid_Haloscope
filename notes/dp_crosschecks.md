# Dark-photon cross-checks (Letter v13)

Each entry: quantity — method 1 — method 2 — agreement.

1. Conversion constant chi = 3.38e-15 (g/1e-12)(B/10 T)(meV/m)(1/(3 kappa))^1/2
   - analytic: chi*cos(theta) = g B0/m from E0 = g a0 B0 <-> chi sqrt(2 rho)  -> 1.95e-15*sqrt(3) = 3.38e-15
   - numeric: P/A = kappa chi^2 rho c = 1.92e4 W/m^2 kappa chi^2 reproduces Eq. (1) 7.33e-26 W/m^2 at chi=1.95e-15
   - independent assessment (user PDF, 19 Oct 2026): 3.3835e-15. Agreement < 0.1%.
2. All-space vs confined drive (15 settings): scan response change -0.5% to -6.0%, chi change <= 1.57% (sheet model, this work);
   independent assessment <= 1.53%. Agreement.
3. Polarization, fixed scenario: handbook (Caputo et al. 2021, Fig. 4, integer sidereal days, 46 deg): North 0.316, Zenith 0.306, West 0.186;
   Monte Carlo (this work, 2e5 polarizations, random window start, windows 0.9-4.3 d): split parallel/perpendicular-to-axis array, 5th percentile 0.299-0.314;
   single perpendicular orientation 5th percentile 0.049; north-pointing wires only: 5th percentile 0.237/0.236/0.250/0.253/0.257 for the five exposures (0.261 for >=1 d), i.e. 13-19% in chi;
   independent referee MC: 0.23-0.26 (north), 0.315-0.317 (split). The handbook factor (0.316) is a less strict criterion than the 5th percentile; the Letter quotes both.
   Data column cnt_fixed in dp_reach.csv uses kappa=0.299 (split array, 5th percentile).
4. Prototype, one tile 1-2 meV cartridge, KID NEP 3e-19, 1 day fixed: P_min = 3.61e-21 W, chi = 3.3e-13 (this work and assessment).
5. Dish formula vs BREAD Eq. (2)/Fig. 2: our expression with BREAD inputs (one sensor, NEP 2e-19, 1000 d, eta 0.5, kappa 2/3, BREAD's rho=0.45 GeV/cm^3) 3.25e-14 vs BREAD 3.07e-14 (digitized) / 3.1e-14 (Eq. 2). Agreement 5%.
6. Stack vs BREAD QCDet scenario (equal total dark rate 4/s, eta 0.5 detector, 10 m^2): published BREAD 1.72e-15, 2.4e-15, 3.42e-15 at 2,4,8 meV vs ours 6.1e-16, 1.08e-15, 1.68e-15 -> 2.8, 2.2, 2.0.
   Analytic ratio^4 = (kappa_s eta_c)^2 I (T_oct/ln2)/(kappa_d^2 T_d) with kappa_d=1 (BREAD normalization), I~6e3 at 2 meV -> 3.3. Consistent.
7. SQWARE identical-detector ratio: scaling (axion ratio x B_SQ/10T x sin(theta)/sqrt2) gives 10.2 (36 T) and 2.9 (50 T);
   independent common-scan reconstruction with SQWARE detectors (assessment) 10 and 1.4, which converted to our detectors (eta 0.5, 4 mHz) with eta^-1/2 D^1/4 gives 9.2 and 2.95: agreement 10% and 2% (v13 referee caught an inverted conversion, 11.5, in the earlier note).
8. SQWARE normal-field requirement at 1 meV: from scan of resonance vs tilt (this work) B_perp = 17.7 T (config 2), 26.7 T (config 3); assessment 17.7 / 26.7 T. Exact agreement.
9. Literature: all limits/projections verified in primary papers (dp_literature.md); AxionLimits 'Karshenboim.txt' is a mislabelled fifth-force bound and is not used;
   quantum-cyclotron limit taken from the paper (3.2e-11, 90% CL, rho=0.3) rather than AxionLimits (4.0e-11).
10. Independent referee pass (v13, dp_referee_report.md): mapping 3.3836e-15, prototype 3.305e-13, BREAD Eq.(2) 3.13e-14, dp_reach.csv = mapped axion curves to machine precision,
    handbook factors 0.316/0.306/0.186/0.024, N_modes and dish ratios, SQWARE 10.17/2.91, BREAD QCDet ratios 2.81/2.22/2.03 -- all reproduced with independent code.
11. Competitor projections added after verification in the primary papers: highly excited cyclotron eps ~ 2e-16 at 0.1 meV (2410.05549, text);
    Rydberg-BREAD (2511.00145, Fig. 3, DCR 0.01 Hz): bare dish (1-3)e-15, with dielectric 1e-19 - 1e-17 (read from figure); Rydberg transducer 2603.23337 (g ~ 1e-13 at 0.4 meV, 10 d);
    2607.15612 (5e-5 - 7e-4 eV); near-cutoff waveguide 2605.15820 (2.1e-15 near 0.1 meV).
12. v14 additions: halo factor wind dependence (own MC): normal wind 0.967/0.783, in-plane 0.932/0.615 at 1.41/4 meV (referee 0.967/0.932, 0.78/0.62).
    SQWARE cfg 3 Poisson (Asimov, Lorentzian response, 0.062 background counts per FWHM): threshold chi x1.27 -> ratio 1.78 (technical review: 2.26 with full response); quoted 1.8-2.3.
    One-tile KID curve: chi_kid(array) x (N_tile x 12)^(1/4), N_tile = 1000 (10 cm) below 4 meV, 4000 (5 cm) above -- v13 used 1000 above 4 meV (curve 1.41x too low there); fixed in v14.
13. v14 review: SQWARE axion ratio direction verified in sqware_headtohead.csv (g_ratio_corrected = g_SQ/g_ours: 4.20, 0.95; half collection 5.94, 1.34). Poisson diagnostic script sqware_poisson_check.py: 1.78 (Lorentzian); review's integrated-response value 2.26.
