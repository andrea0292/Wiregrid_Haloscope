# MADGRID: layered wire-grid haloscope for millielectronvolt axions and dark photons

Code, data and paper source for the Letter *"Layered Wire-Grid Haloscope for Millielectronvolt Axions and Dark Photons"*
(A. Caputo). MADGRID (**M**illielectronvolt **A**xions and **D**ark photons with layered **GRID**s) is built from
mirror-backed stacks of free-standing copper wire grids, with wires parallel to B0, tuned by their spacing. Each stack acts
as an epsilon-near-zero slab that radiates through its face into external optics. Five octave cartridges cover 0.25–8 meV.
The reference array is 1 m² at 10 T, with a 0.1 K converter and one live year per octave. Detector benchmarks in the
Letter: heterodyne N_sys = 10ħω (up to 4 meV), KIDs with NEP = 6×10⁻¹⁹ W/√Hz per channel (the white amplitude-readout
noise of Baselmans et al. 2017), and photon counters with
η_d = 0.5 and D_dark = 1 mHz per channel.

Everything here runs on a laptop. The heaviest calculations are the finite-tile solver (~45 min, 4 GB) and the 3D
FDTD check (~35 min per cartridge on 10 cores).

## Layout

```
paper/        main.tex (current Letter + End Matter + SM; BibTeX file bib.bib) and the figures it includes:
                fig1.pdf, fig2.pdf, fig_axion_1m2_context.png (Fig. 3), fig_dp_1m2_context_proj.png (Fig. 4),
                fig_overlap_modes.png (Fig. S1); other figure files there are variants not included in the Letter
code/         solvers, figure scripts and checks (run scripts from code/ unless stated otherwise)
  reproduce_final.py      impedance-sheet solver, tuning trajectories, halo factor chi_their(), reach formula
  make_reach_tables.py    END-TO-END driver: five cartridges + wire correction + halo + detectors -> reach tables of Fig. 3
  nsweep_fig2_fig4.py     optimal number of grids (Fig. 2a data: data/nsweep.csv)
  figures_tuning.py       Fig. 2 (paper/fig2.pdf) from data/nsweep.csv and data/fig2b_tuning_paths.csv
  fig_madgrid_context.py  Fig. 3: axion reach over existing bounds, 0.01–20 meV (paper/fig_axion_1m2_context*.png|pdf)
  fig_madgrid_dp_context.py  Fig. 4: dark-photon reach over existing bounds (paper/fig_dp_1m2_context_proj.* is the Letter's)
  fig_madgrid.py          the same curves on a 0.25–8 meV frame (paper/fig_axion_1m2.*, fig_dp_1m2.*; not in the Letter)
  figures_reach_1m2.py    superseded 1 m² figures with the earlier counter benchmark (4 mHz per 10 cm tile)
  figures_reach.py        superseded figures with 10 m² as the reference (shared plotting helpers still imported)
  figures_slab.py         slab boost vs thickness (paper/fig5.*; not included in the current Letter)
  fullwave.py             circular-wire multipole + lattice-sum solver (S1)
  fw_reconcile.py         wire vs sheet comparison at the sampled settings, 0.25–4 meV (S1)
  fw_reconcile_4_8.py     same for the 4–8 meV cartridge
  flatness_mc.py          grid-flatness Monte Carlo (tolerance table)
  finite_tile/            wire-resolved finite tile: radiation, Gaussian collection, halo, tilts (S2)
  comparator/             fixed-envelope port-coupled comparator model (supporting study; not quoted in the current Letter)
  finite_collection/      independent flanged-enclosure mode-matching solver (S2 collection 0.72-0.78) and distributed-port
                          comparison incl. volcheck2.py (equal-volume comparison); see its README.txt / MODEL.txt
  dark_photon/            dark-photon mapping, polarization Monte Carlo, prototype, BREAD normalization
  sqware_check.py         independent recomputation of the SQWARE benchmarks
  sqware_poisson_check.py SQWARE vs wire grid at identical detectors
  validation/
    run_fdfd_crosscheck.py + fdfd_check.py   independent finite-difference check of the multipole solver
    gap_jitter.py                             coupling penalty from random common-gap jitter
    overlap_receiving_mode.py, overlap_modes_figure.py   boost as overlap with the receiving mode (SM S3)
    overlap_modes_figure_im.py                 Fig. S1 of the Letter (real and imaginary parts of the receiving modes)
    endmatter_checks.py                        Q_u and l_amp (Introduction), per-element boost and dish comparison (End Matter B)
    fdtd3d/                                   3D FDTD check of finite wire length and frame termination
data/         every number quoted in the paper, as CSV; data/limits/ holds external limits and projections
data/limits/axionlimits/  AxionLimits files used for the context figures (commit 7d375f4), with proj/, dp/ and dp/proj/
data/fdtd3d/  reduced 3D FDTD results (spectra, aperture profiles, summary table)
results/      figures from the validation runs that are not in the Letter
notes/        literature tables, SQWARE and dark-photon cross-checks, reference check, mechanism drawing
docs/         numerical_methods.md: the equations each solver integrates and how they are validated
```

## Quick start

```bash
conda env create -f environment.yml && conda activate wiregrid
cd code
python reproduce_final.py          # Table I designs + sheet-model reach              (~6 s)
python make_reach_tables.py --compare   # full reach tables of Fig. 3, compared with data/ (~40 s)
python fig_madgrid_context.py      # Fig. 3 -> ../paper/fig_axion_1m2_context*.png|pdf  (~2 s)
python fig_madgrid_dp_context.py   # Fig. 4 -> ../paper/fig_dp_1m2_context*.png|pdf     (~2 s)
python figures_tuning.py           # Fig. 2 -> ../paper/fig2.*                         (~1 s)
python figures_slab.py             # slab-boost figure -> ../paper/fig5.* (not in the Letter)  (~1 s)
cd dark_photon && python dp_tools.py && cd ..   # dark-photon table -> ../../data/dp_reach_check.csv
cd validation && python overlap_modes_figure_im.py && python endmatter_checks.py   # Fig. S1, End Matter numbers
```
or simply `bash run_all.sh` (fast items only; add `--slow` for the full-wave, FDFD, finite-tile and FDTD runs).

To compile the paper: `cd paper && pdflatex main && bibtex main && pdflatex main && pdflatex main`.

## Where each result in the paper comes from

| Paper element | Script | Data |
|---|---|---|
| Fig. 3 (axion reach, 1 m² reference; 10 m² faint) | `code/make_reach_tables.py` → `code/fig_madgrid_context.py` | `data/reach_0p25_8_fw.csv`, `data/reach_0p25_8_fw_Rs2.csv`, `data/limits/axionlimits/` |
| Fig. 4 (dark-photon reach) | `code/fig_madgrid_dp_context.py` (`_proj` version) | `data/dp_reach.csv` (columns from `dark_photon/dp_tools.py`), `data/limits/dp_projections_digitized.csv`, `data/limits/axionlimits/dp/` |
| Abstract / Sensitivity: crossings (heterodyne 1.0 meV; KIDs 1.2 and 1.4 meV; counters DFSZ-I 0.8–1.2, KSVZ 0.8–3.9 meV), Eq. (5) 1.30×10⁻¹³ GeV⁻¹, χ down to 3.4×10⁻¹⁶; End Matter B: KIDs about four times weaker than BREAD near 1 meV | `code/quoted_numbers.py` (the Fig. 3/4 curves; see "Benchmark scaling" below) | as above |
| Fig. 2(a) optimal number of grids | `code/nsweep_fig2_fig4.py`, `code/figures_tuning.py` | `data/nsweep.csv` |
| Fig. 2(b) tuning paths | `code/figures_tuning.py` | `data/fig2b_tuning_paths.csv` |
| Table I (cartridge designs) | `code/reproduce_final.py` (writes `design_check.csv`) | |
| Fig. 1 concept (`paper/fig1.pdf`) | drawn interactively; no generating script | |
| Introduction: Q_u ≈ 510 (m_a/meV)^(−2/3), l_amp = 6.3 / 0.4 mm | `code/validation/endmatter_checks.py` | printed |
| S1: sheet vs wire, 0.4–4.7 % shift, sqrt(I_wire/I_sheet) = 0.92–1.05 | `code/fw_reconcile.py`, `code/fw_reconcile_4_8.py` | `data/fullwave_reconciliation*.csv` |
| S1: independent FDFD agreement | `code/validation/run_fdfd_crosscheck.py` (~15 min) | `data/fdfd_crosscheck.csv` |
| S1: dark-photon source outside the stack | `code/dark_photon/dp_tools.py` | `data/dp_reach.csv` |
| S2: finite 10 cm tile, 96 % forward power, 74–78 % Gaussian collection | `code/finite_tile/run_finite_tile.py` | `data/finite_tile_results.csv` |
| S2: finite wire length and frame termination (3D) | `code/validation/fdtd3d/` | `data/fdtd3d/` |
| S2: flanged-enclosure collection 0.72–0.78 | `code/finite_collection/` | `code/finite_collection/results/` |
| S2: halo factor ξ | `reproduce_final.chi_their` (see docs §8) | |
| S3 / Fig. S1: overlap with the receiving mode | `code/validation/overlap_receiving_mode.py`, `overlap_modes_figure_im.py` | `paper/fig_overlap_modes.png` |
| Table II (gap stability, rms grid error) | `code/flatness_mc.py` | `data/flatness_tolerance.csv` |
| Gap-jitter penalty (supporting; not in the current Letter) | `code/validation/gap_jitter.py` | printed |
| End Matter B: SQWARE comparison | `code/sqware_check.py`, `code/sqware_poisson_check.py` | `data/sqware_results.csv`, `data/sqware_headtohead.csv` |
| End Matter B: per-element boost (2.5 per sapphire disk, 2.6 / 1.3 per grid), BREAD comparison (2.8–6) | `code/validation/endmatter_checks.py`, `code/dark_photon/dp_tools.py` | printed |
| Port-coupled comparator (supporting; not in the current Letter) | `code/comparator/` | `data/comparator_results.csv`, `data/final_comparison_equal_volume.csv` |
| Detector scenarios | `code/reproduce_final.py` | `data/detector_scenarios.csv` |

### Benchmark scaling of the stored reach tables
`data/reach_0p25_8_fw*.csv` and `data/dp_reach.csv` are computed for a 10 m² array (columns `lin10`, `kid10`, `target`;
`het`, `kid`, `cnt`), with counters at 4 mHz per 10 cm tile below 4 meV and 1 mHz per 5 cm tile above. The Letter's
1 m² curves follow from the exact scalings g, χ ∝ (N_tile T)^(−1/4) and ∝ (D/η_d²)^(1/4) (thermal term negligible):
all curves × 10^(1/4), and counters below 4 meV additionally × (1/4)^(1/4) for D = 1 mHz per channel. `fig_madgrid*.py`
apply these factors (`F` and `cfac()` in `fig_madgrid.py`); the faint 10 m² curves are the same tables with the counter
factor only. The KID columns of the stored tables assume NEP = 3×10⁻¹⁹ W/√Hz per channel; since KIDs enter as a dark
rate D = NEP²/2(hν)² and the thermal term is negligible, g, χ ∝ NEP^(1/2), and the figures scale every KID curve
(including the R_s × 2 bands, the faint 10 m² lines and the one-tile curve) by `KIDF` = (KID_NEP/3×10⁻¹⁹)^(1/2) = √2 for
the Letter's KID_NEP = 6×10⁻¹⁹ W/√Hz (set in `fig_madgrid.py`). `code/quoted_numbers.py` prints the reach numbers quoted
in the Letter from the same scaled curves.

### Reproduction status
- **Figures:** Figs. 3, 4 and S1 regenerate pixel-identically to the files included in the Letter; Fig. 2 regenerates
  pixel-identically to `paper/fig2.*` from `data/nsweep.csv` and `data/fig2b_tuning_paths.csv`.
- **Reach tables (Fig. 3):** `make_reach_tables.py` rebuilds them from scratch (designs, tuning, sheet response, circular-wire
  correction, halo factor, detectors), and the stored `data/reach_0p25_8_fw*.csv` are now its exact output. At the 4 meV
  junction one cartridge contributes on each side, and the grid point m = 4.000 meV belongs to the 4–8 meV cartridge in every
  column (earlier tables mixed the two cartridges at that single point). The heterodyne curves are plotted only below 4 meV.
  `--overlap-at-4` would let both cartridges contribute near 4 meV, as at the other junctions; this is not used because the
  tile size and counter benchmark change there.
- **Dark photon (Fig. 4):** `dark_photon/dp_tools.py` maps the axion tables to χ and writes `data/dp_reach_check.csv`; the
  columns of `data/dp_reach.csv` used in the figures (het, kid, cnt, their R_s × 2 versions, proto_kid, dish_*) are copied from it.
- **Halo factor:** all forecasts use `chi_their()`: Gaussian mode truncated at the tile, mean wind in the tile plane
  (docs §8), as described in SM S2.
- **Not distributed:** the primary PDFs used by `dark_photon/build_dp_tables.py` for digitization (outputs are in `data/`).

## 3D FDTD check of finite wire length (`code/validation/fdtd3d/`)

This check is a numba Yee solver. It is periodic across the wires and finite (10 cm) along them. The model has a PEC mirror,
resistive wire lines with the paper's copper R_s, and graded absorbing layers. The axion drive is a uniform J_y confined to the
stack. Three wire terminations are compared with infinite wires (y-periodic): thin conducting bars, conducting side walls of
stack height, and open ends.

```bash
cd code/validation/fdtd3d
python run_cases.py A          # 1–2 meV cartridge (N = 11), ~35 min
python run_cases.py B          # 0.25–0.5 meV cartridge (N = 15), ~30 min
python reduce.py               # -> data/fdtd3d/*.npz and fdtd3d_finite_wire_summary.csv
python make_figure.py          # -> results/fig_fdtd3d_finite_wires.(pdf|png)
```
Raw field records (~200 MB) go to `raw/` and are not version-controlled.
The infinite-wire FDTD agrees with the sheet model to 1–3 % in β² and Q_r at the same mesh (Nx = 6, r_eff/a = 0.031).
The geometry is periodic across the wires (it is not a complete square tile), the wires are resistive mesh lines with an effective
radius, the drive is uniform (no halo velocity), and only two cartridge centres at a single mesh resolution were run.
The quoted √(𝓘_G/𝓘_∞) is an effective power-collection efficiency along the wires (ideal uniform aperture: 0.89).

## External data
`data/limits/` contains CAST, globular-cluster, NuSTAR and dark-photon bounds from the AxionLimits compilation
(C. O'Hare, https://github.com/cajohare/AxionLimits) and projections digitized from the cited papers
(`data/dp_projections_digitized.csv`). Please cite the original experiments when reusing them.
