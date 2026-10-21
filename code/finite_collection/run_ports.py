from pathlib import Path
import sys,json,time
import numpy as np
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from distributed_port import PlasmaPort,grid_shapes,single_mode_validation

def run(m,L,label):
    name=P/'results'/f'ports_{m:g}_{label}.json'
    if name.exists():
        cached=json.loads(name.read_text())
        if cached.get('model_version')==5:return cached
    obj=PlasmaPort(m,L);rows=[]
    shapes={K:grid_shapes(K) for K in [1,4,16,64]}
    for K,ss in shapes.items():
        candidates=[]
        for nx,ny in ss:
            out,b,u=obj.optimize(nx,ny)
            candidates.append(out)
        best={**max(candidates,key=lambda x:x['I']),'candidate_geometries':candidates}
        best['radiative_I']=json.loads((P/'results'/f'radiative_{m:g}.json').read_text())['finite']['I_gaussian']
        best['R_equal_downstream_noise']=best['radiative_I']/best['I']
        rows.append(best)
        print('PORT',m,label,K,best['nx'],best['ny'],'I',best['I'],'R',best['R_equal_downstream_noise'],flush=True)
    full,b,u=obj.optimize(uniform=True);full['R_equal_downstream_noise']=rows[0]['radiative_I']/full['I']
    out=dict(model_version=5,mass_meV=m,geometry=label,L_mm=L*1e3,rows=rows,continuous_profile=full)
    name.write_text(json.dumps(out,indent=2));return out

if __name__=='__main__':
    print('single_mode',single_mode_validation(),flush=True)
    for m in [.5,1.,2.]:
        r=json.loads((P/'results'/f'radiative_{m:g}.json').read_text())
        run(m,r['best']['L'],'same_active_volume')
        run(m,.03,'30mm_available_depth')
