# Numerical methods

The haloscope is described by four calculations of increasing fidelity, plus one independent cross-check, which together
validate the effective-medium picture. This file gives the equations each one integrates, its approximations, and what it was
checked against.
Units are SI. Fields are normalized to the field E0 that the axion induces in vacuum (E0 = −g_aγγ a B0 in natural units), with E0 = 1.
With this normalization |β|² is the power boost relative to an ideal mirror ("dish") of the same area.

## 1. Impedance-sheet model (`code/reproduce_final.py`)

**Geometry.** There is a perfect mirror at x = 0, and N grids at x_j = j d (j = 1..N), forming a slab of thickness L = N d. The axion drive
is confined to 0 < x < L. Each grid is a sheet of complex impedance per square

    Z_g = R_g − i ω L_s ,  L_s = (μ0 a / 2π) ln(a / 2π r) ,  R_g = R_s a / (2π r)

with pitch a, wire radius r, and copper surface resistance R_s = 83.6 mΩ (ν/250 GHz)^(2/3) (anomalous skin effect, RRR-limited).

**Field equations.** In each gap the total field is E = E_p + A_j e^{ikx} + B_j e^{−ikx}, where the particular solution E_p = −E0
is present only inside the driven region. At the mirror E(0) = 0. At each grid E is continuous, and the derivative jumps because of the
sheet current K = E/Z_g:

    ∂E/∂x |_{x_j^+} − ∂E/∂x |_{x_j^−} = −(i ω μ0 / Z_g) E(x_j)

Above the stack only the outgoing wave A_N e^{ikx} survives, and β ≡ A_N e^{ikL}. This gives a (2N+1)×(2N+1) linear system per
frequency.

**Derived quantities.**
- Peak boost |β_pk|² and response quality factor Q_r = ν_pk / Δν_FWHM, from the half-maximum points of |β(ν)|².
- Scan integral 𝓘 = ∫|β|⁴ dν/ν, which sets the scan rate at fixed SNR. The reach scales as 𝓘^{−1/4}.
- Material quality factor Q_m = ω_p/Γ, with ω_p² = 2πc²/(a² ln(a/2πr)) and Γ = R_s/(μ0 r ln(a/2πr)).
- Radiative and ohmic rates, γ_rad and γ_ohm, from Q_r together with the lossless Q_rad (R_s → 10⁻⁴ R_s). Their ratio is x = γ_rad/γ_ohm (Fig. 2a).

**Designs.** For each octave, N is chosen to maximize 𝓘, and a (with d = a, r = 0.03 a) is solved so that the peak sits at the
geometric band centre. Tuning is done by changing d at fixed a, r, N, along a trajectory that is interpolated in log m.

**Reach.** Results are given at SNR = 5 (local), with the dwell time from a log-uniform scan over each octave and one live year per octave.
- Heterodyne: noise 10 hν, bandwidth 10⁻⁶ ν.
- Counters and KIDs: dark rate D per tile plus the thermal rate of the converter at T. A KID enters as a counter with
  D = NEP²/2(hν)².
- Single-tile signal power: P1 = 7.33×10⁻²⁶ W × ½ η_ap ξ (A_tile/m²)(meV/m_a)², with η_ap = 0.65 and ξ the halo overlap of
  S2 (Gaussian mode, w = W/2).

`reach_from()` implements this. Tiles add incoherently, so g ∝ (A a_t T)^{−1/4} at a fixed background per detector.

## 2. Circular-wire full-wave model (`code/fullwave.py`, `code/fw_reconcile*.py`)

Each grid is an infinite periodic row of circular wires. The field scattered by a wire is expanded in cylindrical harmonics with |n| ≤ 3.
The boundary condition is a Leontovich surface impedance, Z_s = R_s(1 − i√3) (anomalous skin effect).
- **Within a grid:** coupling through lattice sums, using an Ewald-type split with asymptotic tails and Twersky's S0.
- **Between grids:** coupling through propagating and evanescent Floquet orders, including mirror images.
- **Drive:** the confined axion drive enters through the same particular solution as in §1, and the passive reflection subtracts the
  external-source term (`fw_conf`).

**Checks.**
- Harmonic cutoff |n| ≤ 3 → 4 gives negligible changes.
- An independent FDFD calculation (§5) agrees to < 1 %.
- At the sampled settings, the resonance shifts by −0.4 to −4.7 % and is retuned. After retuning, √(𝓘_wire/𝓘_sheet) = 0.92–1.05.
  This correction is interpolated in log m and folded into `data/reach_0p25_8_fw*.csv`.

## 3. Finite tile in the plane transverse to the wires (`code/finite_tile/`)

The tile is 2D in (x, z), with N_w wires per grid. Each wire is an outgoing cylindrical monopole that satisfies the impedance condition. The
finite mirror is represented by piecewise-constant surface currents (boundary-element method). Along the wires the treatment is paraxial.
The outputs are the forward power relative to the periodic model (≈ 0.96 for a 10 cm tile), the projection onto a Gaussian
receiver mode (0.74–0.78 of the periodic power), and the effects of tilts and the halo. In the periodic limit, the monopole-only
model is 0.3 % off in frequency relative to §2.

## 4. 3D FDTD: finite wire length and frames (`code/validation/fdtd3d/`)

This calculation tests the one approximation that §3 cannot: wires of finite length, terminated on frames.

**Grid.** Yee leapfrog, cubic cells of size D = a/6. The domain is periodic in x (one wire per grid per period). In y, along the wires,
the wires span W = 10 cm and are followed by free space (60 cells) and an absorbing layer. In z there is a PEC mirror at k = 0, grids at
k = 6(l+1), and free space with an absorbing layer above.

**Wires.** Each wire is a line of E_y nodes with conductivity σ = 1/(R′D²), where R′ = R_s/(2π r_eff) is the resistance per
length. The intrinsic FDTD wire radius r_eff = 0.18518 D (r_eff/a = 0.031) sets the inductance. It was calibrated against the sheet
model using the infinite-wire run.

**Drive.** A uniform current J_y(t), a Gaussian-modulated sinusoid with 4 % bandwidth, is applied for 0 < z < L, with half weight at the top node,
and over the wire span for finite wires. The free-space field of this source is the discrete particular solution
E_p(ω) = Δt J(ω) / (2iε0 sin(ωΔt/2)), and E0 = −E_p normalizes β.

**Output.** On-the-fly DFTs of E_y and H_x are taken 30 cells above the stack and averaged over x. From these:
- β²_tot = Poynting flux / W
- β²_G = |∫E u dy|² / (∫u² dy · W), with u = exp(−y²/(W/2)²). The ideal value for a uniform aperture is 0.89.
- β²_mid = |E(y = 0)|²

**Validation.**
- With a bare mirror, the normalization is reproduced to 1e-5.
- Infinite-wire FDTD vs the sheet model at the same numerical k (lossy, 1–2 meV): f = 1.0073 ν_c, β² 1707 vs 1693, Q_r 252 vs 246.
- Finite wires are always compared with the same-mesh infinite-wire run, so the discretization error cancels.
- Doubling the free-space margin between the frames and the absorber (60 → 120 cells, 0.25–0.5 meV, bars) leaves the peak,
  Q_r and the total-power scan integral unchanged to 1e-4. The Gaussian-collected value changes from 0.803 to 0.798.

**Result** (`data/fdtd3d/fdtd3d_finite_wire_summary.csv`)

| cartridge | W/λ | ends | Δf/ν_c | Q_r/Q_∞ | β²_tot/∞ | √(𝓘_G/𝓘_∞) |
|---|---|---|---|---|---|---|
| 1–2 meV | 114 | bars | <1e-4 | 0.98 | 0.97 | 0.90 |
| 1–2 meV | 114 | walls | <1e-4 | 1.00 | 1.00 | 0.92 |
| 1–2 meV | 114 | open | <1e-4 | 0.98 | 0.97 | 0.90 |
| 0.25–0.5 meV | 28.5 | bars | 2e-4 | 0.92 | 0.86 | 0.80 |
| 0.25–0.5 meV | 28.5 | walls | <1e-4 | 1.00 | 1.00 | 0.90 |
| 0.25–0.5 meV | 28.5 | open | 2e-4 | 0.92 | 0.84 | 0.78 |

Δf/ν_c is quoted on the 1e-4 frequency grid, so "<1e-4" means no shift was resolved. The ideal uniform-aperture value of the last column is 0.89 in one dimension. Conducting side walls make the
wire-direction response indistinguishable from infinite wires. With thin bars or open ends, the lowest cartridge loses a further
≈ 10 % in collected amplitude, which is within the η_ap = 0.65 margin.

## 5. FDFD cross-check (`code/validation/fdfd_check.py`, `run_fdfd_crosscheck.py`)

This is a 2D frequency-domain finite-difference solution with explicitly meshed perfectly conducting wires, a lossless three-grid stack,
r/a = 0.1, and Δx = a/200 and a/400. It solves the homogeneous Helmholtz equation for E_s = E + E0, with E_s = 1 on the conductors and a PML above.
Linear extrapolation in Δx agrees with the multipole solver to ~1e-6 in frequency and 0.4 % in peak β².

## 6. Tolerances (`code/flatness_mc.py`, `code/validation/gap_jitter.py`)

- **Independent grid displacements:** Monte Carlo over random grid positions (400 realizations), with the scan integral recomputed by
  `beta_pos()`. The rms error that raises g by 5 % is 6.0, 3.8, 2.4, 1.4 and 0.7 µm for the five cartridges.
- **Common-gap change:** a static change only shifts the resonance. The tabulated values (0.42–0.13 µm) are defined by a 10 % peak-power loss at
  fixed frequency. Gaussian jitter of rms σ_d = (0.5, 1, 2, 3, 5)× this value raises g by (0.3, 1.2, 4.0, 7.2, 13.4) %.
