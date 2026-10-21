from pathlib import Path
import runpy
from run_radiative import run as radiative
from run_ports import run as ports
P=Path(__file__).resolve().parent
for mass in [.5,1.,2.,4.]:
    r=radiative(mass)
    ports(mass,r['best']['L'],'same_active_volume')
    ports(mass,.03,'30mm_available_depth')
for script in ['validate_comparison.py','highmass_check.py','assess_resources.py','make_figures.py']:
    runpy.run_path(str(P/script),run_name='__main__')
