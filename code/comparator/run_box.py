import sys, os, json, numpy as np, time
sys.path.append(os.getcwd()); import wm_extraction as wm
lam_mm=1.095876665964563; Qu=473; h=0.5; a_p=1.0
out=[]
def run(W,nside):
    f=lambda gm: wm.solve_F_box_sym(W,h,lam_mm,Qu,nside,a_p,gm,ndelta=101)[0]
    Fb,gb=wm.optimise_gamma(f,lo=-3,hi=2.5,n=8)
    return float(np.sqrt(Fb/wm.F_ref(W*W,Qu))), float(gb)
# uniform loading reference per W (zeta=2)
for W in [5,10,20,30,50,75,100]:
    Fu=wm.solve_F_box_sym(W,h,lam_mm,Qu,1,a_p,2/Qu,ndelta=101,uniform=True)[0]
    out.append(dict(kind="uniform",W=W,nside=0,eta_ext=float(np.sqrt(Fu/wm.F_ref(W*W,Qu)))));print(out[-1],flush=True)
cases=[]
for W in [5,10,15,20,30,50,75,100]:
    for ns in [1,2,3,4,6,8,10,13,16,20,25,33]:
        s=W/ns
        if s<3 or (W>=75 and ns<3 and ns!=1): continue
        cases.append((W,ns))
for W,ns in cases:
    t=time.time(); e,gb=run(W,ns)
    out.append(dict(kind="lattice",W=W,nside=ns,N=ns*ns,s=W/ns,eta_ext=e,gamma=gb,t=time.time()-t)); print(out[-1],flush=True)
    json.dump(out,open("box.json","w"),indent=1)
