import sys, os, json, numpy as np, time
sys.path.append(os.getcwd()); import wm_extraction as wm
lam_mm=1.095876665964563; Qu=473
def eta_cell(s, a_p, h, ndelta=121):
    n=int(round(s/h)); hh=s/n
    mask=wm.port_mask(n,hh,[(s/2,s/2)],a_p)
    f=lambda gm: wm.solve_F(s,hh,lam_mm,Qu,mask,gm,bc='N',ndelta=ndelta)[0]
    Fb,gb=wm.optimise_gamma(f,lo=-3,hi=2.5,n=9)
    return float(np.sqrt(Fb/wm.F_ref(s*s,Qu))), float(gb)
out=[]
for a_p in [1.0,0.5,2.0]:
    for s in [2,3,4,5,7.5,10,15,20,30,40]:
        if s<=a_p*1.5: continue
        h=0.25 if s<=20 else 0.5
        if a_p==0.5: h=min(h,0.125) if s<=10 else (0.25 if s<=20 else None)
        if h is None: continue
        t=time.time(); e,gb=eta_cell(s,a_p,h)
        out.append(dict(a_p=a_p,s=s,h=h,eta_ext=e,gamma=gb,t=time.time()-t)); print(out[-1],flush=True)
        json.dump(out,open("cells.json","w"),indent=1)
