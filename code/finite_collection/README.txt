INDEPENDENT FINITE-COLLECTION AND PORT COMPARISON

Read finite_collection_comparison.pdf first. MODEL.txt specifies the port
operator, power normalization, coherent receiver, halo and loss assumptions.
The optional comparison_appendix_insert.tex is a compact insertion for the
existing paper. The existing paper and absolute coupling curves are unchanged.

Dependencies: Python 3, NumPy >=2, SciPy, Matplotlib; LaTeX for the note.
Run from this directory, preferably with OPENBLAS_NUM_THREADS=1 and
OMP_NUM_THREADS=1:

  python reproduce_all.py
  pdflatex -interaction=nonstopmode -halt-on-error finite_collection_comparison.tex
  pdflatex -interaction=nonstopmode -halt-on-error finite_collection_comparison.tex

Radiative and port result files are cached. To force complete recalculation,
move the corresponding results/radiative_*.json and results/ports_*.json
files aside before running. The delivered package contains the final cached
results, raw radiative spectra and verification output. Intermediate diagnostic
versions, temporary files and Library transfer helpers are not included.

The new port implementation is model/distributed_port.py. Its reference is
the axion-driven Drude plasma, not the earlier fixed-envelope comparator (code/comparator/). It uses
one passive coherent output even with many pickup elements. Source power,
ohmic loss and extracted power balance explicitly. validate_comparison.py
contains the independent finite-difference check and convergence studies.

The finite-aperture solver and impedance-sheet solver are copied from the
previously delivered, validated manuscript package. The added work computes
new finite collection benchmarks, the distributed port family, consistent
halo overlaps, available-depth comparisons and a copper-guide loss threshold.

The comparison is a reduced-model calculation, not a claim that these
pickup networks, optics or detector targets are built.
