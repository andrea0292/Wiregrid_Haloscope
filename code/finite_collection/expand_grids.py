"""Complete all regular-grid aspect ratios without repeating cached solutions."""
from pathlib import Path
import json,sys
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from distributed_port import PlasmaPort,grid_shapes
for path in sorted((P/'results').glob('ports_*.json')):
    data=json.loads(path.read_text())
    if data.get('model_version')==5:continue
    obj=PlasmaPort(data['mass_meV'],data['L_mm']*1e-3);new=[]
    for row in data['rows']:
        candidates=row['candidate_geometries'];known={(x['nx'],x['ny']) for x in candidates}
        for nx,ny in grid_shapes(row['pickups']):
            if (nx,ny) not in known:
                rr,b,u=obj.optimize(nx,ny);candidates.append(rr)
        best={**max(candidates,key=lambda x:x['I']),'candidate_geometries':candidates,
              'radiative_I':row['radiative_I']}
        best['R_equal_downstream_noise']=best['radiative_I']/best['I'];new.append(best)
        print(path.name,row['pickups'],'old',row['nx'],row['ny'],'new',best['nx'],best['ny'],'I ratio',best['I']/row['I'],flush=True)
    data['rows']=new;data['model_version']=5;path.write_text(json.dumps(data,indent=2))
