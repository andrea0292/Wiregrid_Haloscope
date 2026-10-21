# SQWARE (Mehrani et al., arXiv:2509.14320v2 / PRL 136, 241001): number check

Main text and SM read in full. Recomputed with `sqware_check.py` from the paper's model: Drude CRI permittivity (m* = 0.067 m_e, ε_bg = 11.56, ε_B = 10.69, n_e = 3×10¹¹ cm⁻², d_QW = 30 nm), Maxwell-Garnett EMT, ω_c = eB cosθ/m*. Our own 1D boundary-value boost solution matches Eq. (3)/(S23) and (S24) to machine precision. ρ = 0.3 GeV/cm³ (as in SM I A).

## Parameters vs recomputation

| | Config 1 (HM) | Config 2 (LM) | Config 3 (LM-S) |
|---|---|---|---|
| D / W / d_B / B / τ | 2 mm / 3 cm / 90 nm / 36 T / 1.7 ns | 10 mm / 3 cm / 150 nm / 36 T / 1.7 ns | 20 mm / 5 cm / 90 nm / 50 T / 4 ns |
| η / Γ_dark / t_obs (t_tot = 300 d, T = 0.3 K) | 7% / 1 mHz / 84 min | 20% / 1 mHz / 34 min | 35% / 0.1 mHz / 7 min |
| Claimed range | 2–5 meV | 0.6–2 meV | 0.6–2 meV |
| ENZ limits, θ=90° / θ=0° (ħΩ = ħ[n_e e²/ε₀m*(ε_B d_B+ε_bg d_QW)]^½) | 6.87 / 0.75 meV | 5.63 / 0.50 meV | 6.87 / 0.54 meV |
| Tilt for range (ours) | 70–86° | 33–77° | 25–76° |
| β_max quoted (m_a, θ) | 70 (5.0 meV, 86°) | 115 (1.6, 72°) | 390 (1.1, 60°) |
| β at quoted θ (ours) | 74 (5.03 meV) | 121 (1.53) | 315 (1.07) |
| Global β_max, open / mirror (ours) | 75 / 106 | 128 / 175 | 316 / 416 |
| β² ("~10⁴ power") | 5.6×10³ | 1.6×10⁴ | 1.0×10⁵ |
| Q_max quoted / ours at quoted θ / ours global max | 6×10³ / 6.4×10³ / 8.2×10³ | 1×10⁴ / 1.16×10⁴ / 1.18×10⁴ | 5×10⁴ / 4.5×10⁴ / 4.9×10⁴ |
| ωτ over range | 5.2×10³–1.3×10⁴ | 1.5–5.1×10³ | 3.6×10³–1.2×10⁴ |
| CRI loss ε_d/ε″ ≈ (ω+ω_c)τ | 2.5–6.0×10⁴ | 4.1×10⁴–1.4×10⁵ | 1.4–4.8×10⁵ |
| Wells / MBE time at 1 µm/h | 1.7×10⁴ / 83 d | 5.6×10⁴ / 417 d | 1.7×10⁵ / 833 d |
| B²A (full B) | 1.17 T²m² | 1.17 T²m² | 6.25 T²m² |
| B∥²A = (B sinθ)²A over range | 1.03–1.16 | 0.35–1.11 | 1.16–5.86 |
| β²B²A over range | 0.9–6.6×10³ | 1.4×10³–1.9×10⁴ | 0.6–6.2×10⁵ T²m² |
| g limit (ours, open-slab β) vs Fig. 4 readout [GeV⁻¹] | 6.0/5.6/8.5×10⁻¹² at 2/3.5/5 meV vs 6.1/5.7/8.8 | 5.7/3.8/7.3×10⁻¹³ at 0.6/1/1.6 meV vs 5.8/3.7/6.4 | 1.4/1.0/5.2×10⁻¹³ (3-photon) at 0.6/1/2 meV vs ~0.8/0.7/5.6 |
| N steps × t_obs, step = m/Q_max | 210 d | 207 d | 309 d (claimed 300 d) |

Their worked example (θ = 72.437° → 1.555 meV, SM I F) is reproduced exactly; Configs 1–2 to ≤10%; Config 3 β_max is 20% below theirs.

**Signal (d)**, Config 2, 1.6 meV, KSVZ (6.3×10⁻¹³ GeV⁻¹), β = 115: P = β²g²B²ρA/m² = 1.3×10⁻²⁴ W; 1.0×10⁻³ detected photons/s, i.e. 2.1 signal vs 2.0 dark counts per 34 min. 2σ on KSVZ needs 64 min (36 min at ρ = 0.4); limit 7.4×10⁻¹³, so Config 2 touches KSVZ only near 1 meV. Config 3 (1.1 meV): DFSZ rate 0.024 s⁻¹, 3 photons in ~2 min < 7 min, so the DFSZ claim follows from its inputs.

## Issues (page/section)

1. **Error, Eq. (6), p.6:** its own ingredients (Γ = η|E|²W²/2ω, a₀ = √(2ρ)/m, 2√(Γ_d/t)) give prefactor 7.0×10⁻¹⁴, not 5×10⁻¹⁴. Fig. 4 matches 7.0×10⁻¹⁴.
2. **Error, SM I A:** "a₀ ≈ 1.52 eV (m_a/meV)⁻¹" is √ρ/m; √(2ρ)/m = 2.15 eV. Fig. 4 unaffected.
3. **Typo, p.6:** 30-d lines "at 1 and 3.5 meV for Configs 1 and 2" are swapped relative to the figure.
4. **Minor, Fig. S1:** quoted ENZ points 4.23/43.9 meV; we get 3.73/47.7 meV (43.9 < ħω_c = 44.0 meV is impossible for CRA).
5. **Inconsistency, SM III A:** "weeks to months" growth at 1 µm/h; actually 83–833 d. Stacks are 200–2000× thicker than the 0.01 mm state of the art (they propose wafer bonding as an alternative).
6. **Magnet, p.5/SM III B:** "1 ppm over a 3-cm sample region", but the 36 T SCH is specified at 1 ppm over a 10 mm DSV in a 40 mm bore (Bird et al., IEEE TAS 28, 8171763). A tilted 3×3×1 cm slab needs 3.2–4.3 cm of transverse clearance before the 0.3 K cryostat, lens and detector. The 50 T / 5 cm / 1 ppm magnet is a stated assumption. Scan steps are 0.002° (SM I F).
7. **Key assumption, p.4–5:** the THz damping time is taken to be the dc transport τ = μm*/e. Measured low-T cyclotron-resonance decay in ultrahigh-mobility GaAs 2DEGs is superradiance-dominated and uncorrelated with mobility (Zhang et al., PRL 113, 047601). The operating point sits at ν ≈ 0.3–5 (quantum Hall regime). No ns THz damping is demonstrated, and every β and Q depends on it.
8. **Assumptions, p.5/SM III C–D:** focusing is 100%. The projections use the infinite-area β, although Fig. S6 retains only ≥50% at W/D = 3 (up to 2× in g for Configs 2–3). The CRI emission is circular, so a linear bowtie collects ½ (our inference). η = 20–35% has not been demonstrated together with 0.1–1 mHz (1 mHz was shown at η ≈ 1%). The B-tuned QD detector would share the cavity's field.
9. **Minor, SM I F:** steps of one β-FWHM leave 25% of peak power at the step edges. The 2√B criterion with 2–5 dark counts is not a rigorous 95% CL.

## Verdicts

- **S1: correct with edits.** "SQWARE proposes an ENZ resonance of the cyclotron-resonance-inactive magnetoplasma mode of a mm–cm GaAs/AlGaAs multiple-quantum-well stack, read out by focusing its surface emission onto one THz quantum-dot photodetector; a backing mirror is optional and not used in the projections."
- **S2: correct only under their assumption.** "Taking the dc transport time (1.7–4 ns) as the THz damping time, SQWARE has ωτ ≈ 2×10³–10⁴ and boost-resonance Q ≈ 6×10³–5×10⁴, 10–100× our copper-grid Q of 170–900."
- **S3: incorrect, because tilting is a moving part.** "It tunes by rotating the stack in a fixed field (the existing 36 T NHMFL hybrid; a hypothetical 50 T magnet for Config 3), without changing the cavity geometry."
- **S4: correct for full B.** B²A = 1.2 (Configs 1–2) and 6.3 T²m² (Config 3). With B∥ the values are 0.35–1.2 and 1.2–5.9 T²m², i.e. 160–2900× below 1000 T²m². For fairness, note their β² ≈ 10³–10⁵.
- **S5: correct as a union.** "0.6–5 meV (0.15–1.2 THz) over three configurations (2–5 and 0.6–2 meV); only the long-term Config 3 reaches the QCD band."
- **S6: correct, as their stated assumption.** Add that they cite ≥90% achievable and ~70% demonstrated, and also use the infinite-area boost.
