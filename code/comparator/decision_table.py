import pandas as pd, json, numpy as np
df=pd.read_csv('comparator_results.csv')
def Rmax(d,scn,es,fR=1.0): 
    x=d[(d.port_scenario==scn)&(d.eta_s==es)&(d.fR_surface_resistance==fR)].R
    return float(x.max()) if len(x) else np.nan
rows=[]
for (env,b),d in df.groupby(['envelope','budget'],sort=False):
    r=dict(envelope=env,budget=b,
      R_nom=Rmax(d,'port-nominal',0.68), R_opt_eta0p6=Rmax(d,'port-optimistic',0.6), R_pess_eta0p75=Rmax(d,'port-pessimistic',0.75),
      R_opt=Rmax(d,'port-optimistic',0.68), R_pess=Rmax(d,'port-pessimistic',0.68),
      R_nom_eta0p3=Rmax(d,'port-nominal',0.3), R_nom_eta0p9=Rmax(d,'port-nominal',0.9),
      R_nom_Rs0p5=Rmax(d,'port-nominal',0.68,0.5), R_nom_Rs2=Rmax(d,'port-nominal',0.68,2.0))
    r['gate']='PASS (robust)' if r['R_opt_eta0p6']>=2 else ('FAIL (robust)' if r['R_pess_eta0p75']<2 else 'UNDECIDED')
    r['winner_nominal']='stack' if r['R_nom']>1 else 'port'
    rows.append(r)
t=pd.DataFrame(rows); t.to_csv('decision_table.csv',index=False)
pd.set_option('display.width',260); print(t.round(3).to_string())
