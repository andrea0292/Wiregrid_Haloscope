# Comparator: radiating wire-grid stack vs port-coupled wire medium at m_a = 1.131 meV (273.6 GHz)

**Method.** Counting-limit scan-rate ratio R = R0 (V_s/V_p)^2 F_s/F_p, F = (eta f_coh)^2/N_rx (audit Eqs. 5–6). It assumes equal background per receiving channel, the same detector type for both designs (common QE cancels), and B^2 matched between the two. f_coh is the fraction of the halo-averaged signal in one fixed, wind-compensated uniform-phase mode (Gaussian kernel, l_a = 0.336 m, audit Eq. 4/10). Both designs may split their volume into independent channels or combine coherently, and both obey the same budgets: N_rx receivers and N_el coupling elements (port column couplers, or in-field stack focal cells). Code: `comparator_model.py`, `wm_extraction.py`; all numbers below are taken from `comparator_results.csv`, `decision_table.csv`, `equal_volume.json`, `cells.json`, `box.json` and `R_vs_Nel.json`.

**1. Equal volume (reproduced).** P_dish/A = 5.73e-26 W/m^2 (exact; the 1.1e-26 scaling gives 5.65e-26, a rounding difference). With the audit's b = 2083, Q_s = 314: P_s/A = 1.193e-22, P_p,crit/A = 9.82e-23 W/m^2, R0 = 0.9375, and eta_s/eta_p thresholds of 1.033 / 1.461 / 3.266 for R = 1 / 2 / 10. These match the audit. The audit's b, Q_s differ slightly from ours. Our sheet values (2066/313) give R0 = 0.925 (thresholds 1.040/1.470/3.288); our full-wave values (2125/312) give 0.982. This R0 spread of ±6% is the only discrepancy. Uniform port mode (C = 1): the peak ratio falls from 1.21 to 0.80, as the audit states.

**2. Port-coupled model.** Near ω_p a wire medium has ω² = ω_p² + c²k² for k both along and across the wires. The field envelope therefore obeys a lossy Helmholtz equation with field attenuation length l_em = λ√(2Q_u)/2π = 5.36 mm. The 2D finite-difference solver uses metal side walls, absorbing "column" couplers (1-mm patches) with optimized loading, and the scan metric ∫P²dδ normalized to audit Eq. (1) (ζ = 2, C = 64/π⁴).
- **One port:** η_ext = 1.02 at W = 10 mm. It falls to 0.34 at 20 mm, 0.097 at 30 mm, 0.022 at 50 mm and 0.0056 at 100 mm, i.e. ∝ W⁻² beyond about 2 l_em. This agrees with the ALPHA single-antenna limit d ≲ λ√Q = 24 mm (arXiv:2210.00017, Eq. 32 and Sec. VII).
- **Bulk lattice of spacing s:** η_ext = 1.52 / 1.46 / 1.02 / 0.63 / 0.27 / 0.035 for s = 3 / 5 / 7.5 / 10 / 15 / 40 mm. The value 1.52 = 1/C is the limit for uniform distributed loading. Each coupler drains at most about 40 mm² ≈ 1.4 l_em² of cross-section, independent of W.
- **Ports needed:** η_ext ≳ 0.9 requires N_port ≈ (W/7.5 mm)², e.g. 169–256 ports for W = 10 cm (box runs: 0.88 and 1.13).
- **Couplers along the wires:** each column coupler is a series-fed slotted line running along the wires.
- **Efficiency chain:** η_p = η_ext × η_coup × η_series × 10^(−[n_stages L_j + α·path + exit]/10) × e^(−σ_φ²). Coherent passive combining avoids the √N penalty but pays line loss, junction loss and phase errors. Independent receivers pay 1/N_rx in F. The optimizer picks the hybrid grouping and the used sub-volume.

**3. Loss inputs (274 GHz).**
- Room-temperature WR-3 waveguide: 0.015 dB/mm for CNC copper (WR-3 CNC/SU-8 filter paper, IEEE THz Trans. ~2017; DOI not verified); 0.02–0.07 dB/mm for silicon-micromachined guide (Beuerle et al. 2018, IEEE Xplore 8269833; DOI not verified).
- WR-3 theory with our cryogenic R_s = 88.8 mΩ: 8.6 dB/m. **No measured 4 K WR-3 attenuation was found.** The cryogenic values are scaled from room temperature using the anomalous-skin-effect law (Kurpiers et al., EPJ QT 4, 8 (2017), DOI 10.1140/epjqt/s40507-017-0059-7, covers only 3.5–12.8 GHz).
- Superconducting NbTiN/a-Si microstrip: tanδ = 2.1e-4 at 350 GHz (Hähnle et al., PRApplied 16, 014019 (2021), arXiv:2108.11095). My conversion gives ≈15 dB/m, so this is no better than waveguide.
- Combiners: a 4-way combiner shows about 1 dB, i.e. ≈0.5 dB per binary stage including lines (Zhang et al., JIMTW 35, 451 (2014), DOI 10.1007/s10762-014-0063-2). A WR-3 H-plane combiner shows 1.2 dB including 30 mm of guide (Wang, Zhao, Ding, IEEE Access 9, 22232 (2021); DOI not verified). A 2×2 WR-3 combiner shows 0.5 dB (MDPI Sensors 26, 1965 (2026); simulation and measurement, DOI not verified).
- Corrugated horn coupling: ≈98% Gaussian (Wylde, IEE Proc. H 131, 258 (1984), DOI 10.1049/ip-h-1.1984.0053).
- Port scenarios (optimistic / nominal / pessimistic): α = 4 / 8 / 20 dB/m; L_j = 0.1 / 0.25 / 0.5 dB; η_coup = 0.95 / 0.9 / 0.8; σ_φ = 0.1 / 0.3 / 0.5 rad; exit run in oversized line of 0 / 0.25 / 0.5 × L_axis; active fill 1.0 / 0.8 / 0.6.
- Stack-side quasi-optical combining: 0.15 dB per stage, σ_φ = 0.3 rad. Stack cells that are not at the tube exit pay the same exit loss as the port design.

**4. Envelopes and geometry.**
- **10 cm cube**, 10 T: solenoid-type and dipole-type.
- **MADMAX-like dipole**, 9 T: 1.35 m warm bore (CEA/Bilfinger design; Gardikiotis et al., Ann. Phys. 2024, DOI 10.1002/andp.202300046). We use a usable D = 1.25 m and L = 1.3 m. The sources give a field length of 1.3–2 m (B²A = 100 T²m², arXiv:1901.07401), not the ~1 m in the brief.
- **BREAD-like solenoid**, 10 T: A = 10 m², with L = 2√2 R taken from arXiv:2111.12103, giving R = 0.75 m and L = 2.12 m.
- **Stack configurations:**
  - single face, or back-to-back faces;
  - BREAD barrel lining;
  - dipole hoop-wire barrel, weighted by sin²φ;
  - tube bundles: 4 lined walls in a solenoid, 2 in a dipole. Each tube is divided into BREAD-type focusing cells (L_cell = 2√2 r_cell).
- **Stack fill:** in the cube the optimum fills 14% (one 4-wall cell). Finer cell bundles reach up to 42% but lose on receiver count. In the big magnets the fill is 0.6–1% with few elements; the geometric maximum is 22% (MADMAX-like) and 45% (BREAD-like), reached only with 10³–10⁵ in-field cells.
- **Computed stack point:** η_s = 0.68 (lead's finite-tile value: 0.756 × 0.9 optics), with range 0.6–0.75 and scan 0.3–0.9. The ≤1.5% halo factor of a 10 cm tile is double-counted, which is conservative for the stack.

**5. Results at η_s = 0.68 and R_s × 1.** Values are R_nominal, with the range [port-optimistic, port-pessimistic].

| Envelope | unconstrained | N_rx≤100, N_el≤10⁴ | N_rx≤100, N_el≤10³ | N_rx≤10, N_el≤10² | Gate (R≥2 robust) |
|---|---|---|---|---|---|
| 10 cm cube, solenoid-type | 0.050 [0.009, 0.83] | 0.050 | 0.050 | 0.23 [0.045, 2.8] | fails; undecided only at ≤10² elements |
| 10 cm cube, dipole-type | 0.013 [0.002, 0.21] | 0.013 | 0.013 | 0.057 [0.011, 0.70] | fails (robust) |
| MADMAX-like dipole | 0.004 [0.000, 0.27] | 0.027 [0.005, 3.7] | 0.53 [0.088, 33] | 17 [1.8, 723] | undecided at ≤10⁴; stack favored at ≤10² but not robust |
| BREAD-like solenoid | 0.016 [0.002, 0.69] | 0.27 [0.036, 16] | 3.3 [0.78, 213] | 22.5 [9.9, 1987] | **passes robustly only at ≤10² elements** |

R = 2 is crossed at the following element budgets (N_rx ≤ 100), quoted as nominal [optimistic–pessimistic]: cube 35 [19–135] (solenoid-type) and 14 [11–39] (dipole-type); MADMAX 273 [95–1.9e4]; BREAD 1.6e3 [295–2.8e5].

Sensitivity at nominal settings:
- R_s × 2 raises R; R_s × 0.5 lowers it. For MADMAX at N_el ≤ 10³: 3.0 and 0.18. For BREAD at N_el ≤ 10³: 23.6 and 0.70.
- R scales as η_s². For BREAD at N_el ≤ 10³: 0.64 at η_s = 0.3 and 5.7 at η_s = 0.9.

**Verdict.** The port-coupled resonator wins wherever it can be filled with couplers at spacing ≲ 7.5 mm: at every budget in a 10 cm envelope, and in any envelope if ≳10⁴–10⁵ column couplers are accepted. The unconstrained MADMAX optimum, for example, uses 1.6e5 couplers into 40 receivers. The radiating stack's advantage is resource efficiency: one large face or barrel needs only a handful of receivers and no in-volume couplers. It therefore wins in meter-scale magnets when the coupler budget is ≲10²–10³. The audit's factor-two gate is:
- **not met** in a 10 cm envelope;
- **undecided** in a MADMAX-size dipole: stack favored for N_el ≲ 270 nominal, but not robustly;
- **met robustly** only for a BREAD-size solenoid with ≤ 100 elements and ≤ 10 receivers, where R ≥ 7.7 even with optimistic port losses and η_s = 0.6.

The decisive quantity is the realizable number of low-loss 270 GHz couplers in a wire medium, not the conversion physics.

**Open uncertainties.**
1. No cryogenic 270 GHz loss measurements were found for WR-3, junctions or slotted lines; all are scaled from room temperature.
2. The port model is a 2D effective medium with line-sink couplers. η_ext depends on patch size: 0.51 / 0.63 / 0.84 for 0.5 / 1 / 2 mm at s = 10 mm. Grid convergence is +2–7%, and in small boxes higher modes add about 8%.
3. Series-fed couplers spanning the full wire length (up to 2 m, ~10³ slots each) are idealized, and their fabrication is not costed.
4. Stack tube cells are idealized: grazing reflection off resonant walls, reflector blockage, and in-field focal horns.
5. The background model assumes background ∝ N_rx. A single multimode detector (background independent of mode count) would favor the large-face stack configurations. Thermal emission and detectors operating in the field are not included.
6. The MADMAX usable length (1.3 vs 1–2 m) and the BREAD R/L choice are assumptions.
7. Fixed optics/phasing cannot track the daily halo-wind modulation; this is assumed compensated for both designs.
