#!/usr/bin/env bash
# Regenerate the Letter's figures and tables. Fast items by default (~1 min); --slow adds the full-wave,
# FDFD, finite-tile and 3D FDTD checks (several hours on a laptop).
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
ROOT="$(cd "$(dirname "$0")" && pwd)"
cd "$ROOT/code"
python reproduce_final.py            # design_check.csv, reach_check.csv (sheet model)
python figures_reach_1m2.py          # paper/fig_axion_1m2.*, paper/fig_dp_1m2.*
python figures_tuning.py             # paper/fig2.*
python figures_slab.py               # paper/fig5.*
(cd dark_photon && python dp_tools.py)   # data/dp_reach_check.csv
(cd validation && python gap_jitter.py)
if [[ "${1:-}" == "--slow" ]]; then
  python nsweep_fig2_fig4.py                                     # nsweep.csv
  python flatness_mc.py                                          # flatness tolerances
  (cd validation && python run_fdfd_crosscheck.py)               # ~15 min
  (cd finite_tile && python run_finite_tile.py)                  # ~45 min, 4 GB
  (cd validation/fdtd3d && python run_cases.py A && python run_cases.py B && python reduce.py && python make_figure.py)
fi
echo "done"
