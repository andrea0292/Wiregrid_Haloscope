#!/usr/bin/env bash
# Regenerate the Letter's figures and tables. Fast items by default (~1 min); --slow adds the full-wave,
# FDFD, finite-tile and 3D FDTD checks (several hours on a laptop).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/code"
python reproduce_final.py            # design_check.csv, reach_check.csv (sheet model)
python make_reach_tables.py --compare   # data/reach_0p25_8_fw*_repro.csv vs stored tables
python fig_madgrid_context.py         # Fig. 3: paper/fig_axion_1m2_context*
python fig_madgrid_dp_context.py      # Fig. 4: paper/fig_dp_1m2_context*
python fig_madgrid.py                 # 0.25-8 meV variants: paper/fig_axion_1m2.*, paper/fig_dp_1m2.*
python figures_tuning.py             # paper/fig2.*
python figures_slab.py               # paper/fig5.*
(cd dark_photon && python dp_tools.py)   # data/dp_reach_check.csv
(cd validation && python gap_jitter.py && python overlap_receiving_mode.py && python overlap_modes_figure.py \
   && python overlap_modes_figure_im.py && python endmatter_checks.py)   # Fig. S1, End Matter numbers
if [[ "${1:-}" == "--slow" ]]; then
  python nsweep_fig2_fig4.py                                     # nsweep.csv
  python flatness_mc.py                                          # flatness tolerances
  (cd validation && python run_fdfd_crosscheck.py)               # ~15 min
  (cd finite_tile && python run_finite_tile.py)                  # ~45 min, 4 GB
  (cd finite_collection && python reproduce_all.py)              # cached results in finite_collection/results
  (cd validation/fdtd3d && python run_cases.py A && python run_cases.py B && python reduce.py && python make_figure.py)
fi
echo "done"
