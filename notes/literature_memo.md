# Literature memo: radiative wire-grid ENZ haloscope (0.25–4 meV)

Pages/equations refer to arXiv PDFs; keys to `refs_final.bib`.

## (a) What prior work already established

1. **Tunable wire-plasma conversion is Lawson.** Lawson:2019brd introduced resonant conversion at the wire-medium ω_p (Eq. 7, p3), noting that the enhancement is boundary-independent, "the same would occur for … one surrounded by vacuum" (p3). It proposes tuning through "a series of planes of wires, with spacing adjusted in a single direction", together with swappable inserts (p4). Readout: internal antennas/pickup loops (p4). **Plane-spacing tuning is therefore not new.**
2. **ALPHA identified loss-limited transport and the open-structure option.** Travelling waves decay over ~√(ω/Γ_p) λ_c, and single-antenna readout requires a size below ~√Q λ_c (ALPHA:2022rxj Sec. II.C, p6, Fig. 2). Very large systems "can require multiple antennas" (Eq. 17, p5). The paper also states that a sufficiently large plasma needs no cavity (Sec. VI, p18) and reports free-space wire-array measurements (Sec. IV.A). Bulk formulas: Eqs. (53)–(61), pp9–10. ALPHA Phase I (ALPHA:2026vwc) covers 10–20 GHz with antenna and JPA readout (Table I, p9). Future phases reach 50 GHz with ~100 µm wires at ~1 mm spacing (Sec. IV, pp9–10). No optical extraction is discussed.
3. **The radiating ENZ plasma with a backing mirror is SQWARE.** Mehrani:2025ryc proposes a 2DEG multiple-quantum-well ENZ magnetoplasma at 0.6–5 meV. Its emission leaves the surface and is focused onto one photodetector. An optional mirror redirects both faces (Fig. 1, p2; p3), and the emitted power is enhanced by ~10⁴ (p2).
4. **The Q^{4/3} optimum comes from the polariton papers.** In Schütte-Engel:2021bqm the loss-dominated resonant boost is Eq. (4.50) (p62) and the optimal thickness is d_opt ∝ Γ^{-1/3}, Eq. (4.51) (p64). Marsh:2022fmo Eq. (34) (p8) gives the same result. Together these imply β_max ∝ Γ^{-2/3}, i.e. power ∝ Q^{4/3}. This is the same escape–dissipation balance as the draft's Eq. (slab). It is implicit there, so write "follows from".
5. **The open mirror-backed multilayer read by external optics is MADMAX and DALI.** See Caldwell:2016dcw (Fig. 1, p1), Millar:2016cjp (Area Law, p6 and App. A) and MADMAX:2019pub (β² ≳ 10⁴, B²A = 100 T² m², pp10–12). Knirck (MADMAX:2021lxf, pp4–7) already separates total from Gaussian-beam power (w₀ ≈ ø/3), with up to 30% 3D loss. DALI (DeMiguel:2026ncy) is a mirror-backed magnetized layer stack coupled to a semi-open Fabry–Pérot resonator with a phased-array receiver.
6. **A large face with a single detector is established.** BREAD focuses a 10 m² barrel onto one mm-scale sensor, which must reach NEP < 10⁻²⁰ W/√Hz or DCR < 10⁻³ Hz (BREAD:2021tpx abstract; TeX l.325). Marsh's étendue bound A ≲ 10⁶ A′ (Eq. 33, p8) limits the area one detector can serve.
7. **Free-space wire-array slabs and mirror–grid tuning exist as hardware.** Wooten:2022vpj fits finite wire slabs in free space as Drude dielectric slabs (Eq. 3, p3). The CMB VPMs translate a mirror behind a free-standing grid: PIPER (Chuss:2014piper: 40 µm wire, 117 µm pitch, 39 cm, 1.5 K) and CLASS (Harrington:2018qpe: 50 µm, 160 µm, 60 cm).

## (b) What is plausibly new (hedged)

- **A free-standing metallic wire-grid ENZ converter at 60–967 GHz.** The axion literature I checked stops at ≲100 GHz for wire media: Lawson's soft limit is 100 GHz (p3) and ALPHA's plans reach 50 GHz. I did not survey the general THz-metamaterial literature, so "first wire-medium ENZ at sub-THz" should not be claimed. Claim only *first axion converter of this type in this band*.
- **Shallow mirror-backed operation in which the grid count N is set by the escape–absorption balance**, together with octave cartridges tuned by a common ~3× spacing change. The tuning principle is Lawson's. The calculated octave trajectory and its 0.2–0.4 µm tolerance are the new engineering content. Demonstrated tuning: 16–28% (Kowitt, Lindahl, Balafendiev).
- **Finite-aperture Gaussian collection from a wire-plasma face**, including the halo-coherence factor. The method parallels Knirck and SQWARE (SM Sec. II C); its application to a conductor plasma is new.
- **Trading ALPHA's internal-pickup transport problem at meV for few-receiver optical readout.** (Fixed-envelope caveat of the earlier check applies.)

## (c) Constraints for the reach plot

- **CAST 2024:** g < 5.8×10⁻¹¹ GeV⁻¹, 95% CL, for m_a ≲ 0.02 eV (CAST:2024eil abstract), so it applies over the whole band. (AxionLimits uses 5.67×10⁻¹¹.)
- **Globular clusters (R₂):** g < 0.47×10⁻¹⁰ GeV⁻¹ for semiconvective mixing (Dolan:2022kul abstract). This is not a formal CL: the bound is set where all 20 simulations fall below the observed 95% interval (p8).
- **NuSTAR solar (Ruz:2024gkl):** 7.3×10⁻¹² GeV⁻¹ at 95% CL for m_a ≲ 4×10⁻⁷ eV. It beats ground-based limits up to 3.4×10⁻⁴ eV (p2). The AxionLimits curve gives 4.7–5.1×10⁻¹¹ at 0.25–0.30 meV; this value is provisional. It should be plotted, or at least mentioned, at the bottom of the band.
- **Haloscopes: no exclusion exists in 0.25–4 meV.** The highest-mass haloscope result is ORGAN 1b, g ≥ 4×10⁻¹² excluded at 95% over 107.42–111.93 µeV (Quiskamp:2023ehr). Open-geometry results lie far below the band: MADMAX at ~77–80 µeV (~2×10⁻¹¹, 95%), DALI at 28.54 µeV (9.7×10⁻¹², 90%), GigaBREAD at 44–52 µeV. CADEx (0.33–0.46 meV) has not published a limit. (INSPIRE + AxionLimits search; not exhaustive.)

## (d) Draft claims that are wrong or unverifiable

1. **Shaikhaidarov:2016, "single-photon sensitivity near 0.3 and 0.8 THz": unverified.** The title concerns detecting coherent THz from an HTS Josephson junction; no abstract or text was accessible. Reword using Ikushima:2010, Komiyama:2000, Astafiev:2002.
2. **Ikushima:2010, "low dark rates in this detector family": misleading.** The demonstrated dark rate is 10⁻³ s⁻¹ at 70 mK (0.1 s⁻¹ at 0.3 K), but at 1.5–3 THz, with η = 0.1–1% (pp. 445–450). The draft assumes η_d = 0.5 at 0.8–4 meV. A double-QD device near 0.5 THz exists (p3), but its dark rate is not given.
3. **VerduAndres:2024 is unidentified.** It is replaced by Chambers:1952. The free-electron relation ρℓ = ħk_F/(ne²) gives 6.58×10⁻¹⁶ Ω m² (computed here); Chambers' abstract confirms that the measured σ/ℓ of Cu agrees with theory. Chambers' printed number was not read. Cross-check: 2.78 mΩ at 3.9 GHz vs Posen's 2.9 mΩ at 6 T (p8).
4. **R_∞ prefactor:** the draft labels [√3 μ₀²ω²ρℓ/16π]^{1/3} as the *diffuse*-reflection limit. Reuter–Sondheimer give different specular and diffuse prefactors. Confirm which applies; I could not access RS or the Calatroni note to check.
5. **Calatroni:2020ase:** the report number is confirmed only through Posen's Ref. [39]; the title and content are unverified.
6. **Kerr:2013alma** now points to IEEE TTHz 4, 201 (2014), which gives SSB 60 K in Band 6. That is about 5hν/k at 250 GHz, so N_sys = 10hν is conservative by roughly ×2 (convention caveats as in the draft).
7. **The SQWARE citation understates the overlap.** The draft says "radiating plasma converters have also been proposed". The Letter should state that SQWARE already combines ENZ, a backing mirror, face emission and a single photodetector in the same mass decade.
8. **Lawson's plane-spacing tuning and ALPHA's "no cavity needed" remark (p18) should be credited explicitly** in the main text, not only in the appendix.
9. **Ramanathan:2024hsf** proposes qubit-based phonon detectors with meV thresholds. It is not a THz photon-counting demonstration. Relevant 2024–2026 work: Rydberg-vapour THz single-photon detection, NEP 9.5×10⁻¹⁹ W/√Hz (Li:2024jbz); a ~100 GHz transmon (Anferov:2024ljv); proposals Devlin:2026czx and Banerjee:2025jss. I found no demonstrated low-dark-rate counter at 0.2–1 THz.
