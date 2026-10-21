# Wire-grid radiative plasma haloscope vs SQWARE: what differs once the readout is the same

Both devices use the same principle: an ENZ slab backed by a mirror radiates the axion-induced
field through its face, and external optics focus the emission onto a photon detector. The
readout principle is therefore the same. The difference is in what each converter offers the detector.

## 1. Figure of merit with identical detectors
Background-dominated counting with N_det detectors, each with dark rate D and efficiency eta:

    SNR^2  ∝  t * (eta * P_tot)^2 / (N_det * D),    P_tot ∝ eta_coll * B0^2 * A * beta^2 * (P_dish per unit B0^2 A)

Over a scan, beta^4/Q_r → the scan response I. Merit per unit time, with identical (eta, D):

    M = (eta_coll * B0^2 A)^2 * I / N_det,       g_reach ∝ M^(-1/4)

| comparison (identical detectors)        | m [meV] | g_SQWARE / g_ours (coll 1.0) | (coll 0.5) |
|-----------------------------------------|---------|------------------------------|------------|
| Config 2 (36 T, 9 cm^2) vs ours 1 m^2   | 1.53    | 2.17                         | 3.07       |
| Config 2 (36 T, 9 cm^2) vs ours 10 m^2  | 1.53    | 3.86                         | 5.45       |
| Config 3 (50 T, 25 cm^2) vs ours 1 m^2  | 1.07    | 0.54                         | 0.76       |
| Config 3 (50 T, 25 cm^2) vs ours 10 m^2 | 1.07    | 0.96                         | 1.35       |

(>1 means the wire-grid array reaches smaller couplings.) Source: sqware_headtohead.csv.
Ours: beta^2 = 1594 (1.53 meV) / 1702 (1.07 meV), Q_r = 219 / 309, B0^2 A = 100 / 1000 T^2 m^2, eta_ap = 0.65 and chi included.
SQWARE: beta^2 = 1.46e4 / 9.9e4, Q = 1.16e4 / 4.5e4 (dc tau assumed to set THz damping), B_par^2 A = 1.06 / 4.69 T^2 m^2.

Physical reason: SQWARE gains ~10x in beta^2 (and more in Q) from high-mobility electrons; the grids gain
100-1000x in B0^2 A from area. Merit goes as (B0^2 A)^2 * beta^4/Q_r, so area wins unless the
detector count is penalized — which is why we keep the total dark rate per tile comparable to theirs.

## 2. What an experimentalist gains with grids
1. Magnet: 10 T over ~1-10 m^2 (MADMAX-class dipole, BREAD-class solenoid) exists; 36 T over the sample exists only at
   hybrid-magnet facilities with limited time; 50 T steady state does not exist.
2. Converter: copper surface resistance at 4 K is measured (with an anomalous-skin-effect extrapolation);
   free-standing tungsten/copper grids with 5 μm wires are an established fabrication technique. SQWARE assumes that the THz damping equals the dc
   transport time (untested in this regime) and needs 2-20 mm of MBE growth (~3 months to 2 years at 1 μm/h).
3. Tuning: grid spacing (mechanical, as in polarization modulators, 4 K) at fixed full field along the wires,
   versus rotating the sample in the magnet (field component and chiral projection change).
4. Collection: per-tile reciprocity calculation including finite aperture and halo (eta_ap = 0.65, chi);
   SQWARE quotes infinite-slab emission focused without loss.
5. Polarization: linear, wires along B0 → standard polarizing optics; magnetoplasma emission is elliptical/circular.
6. Detectors outside the field; first science with existing heterodyne receivers (beats CAST to 0.96 meV with 1 m^2)
   and KID arrays (to 1.46 meV with NEP 3e-19), before single-photon counters exist.
7. Scales by tiling; each tile is a tested unit.

## 3. What SQWARE does better
Higher Q and beta^2 (if the dc lifetime holds at THz); a compact sample with a single detector; tuning without
moving the converter's internal geometry; magnetic-field tuning across a wide range.

## 4. What both still need
Single-photon counters at 0.2-2 THz with ~50% efficiency and ~1e-3 s^-1 dark rate per detector.
Quantum-dot counters reach the dark rate, but with an efficiency of 0.1-1%.
