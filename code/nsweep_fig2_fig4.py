"""N-sweep (Fig. 2a), slab boost (Fig. 4) and Fig. 2 plotting for Letter v11.
Requires reproduce_final.py in the same directory (impedance-sheet solver).
Outputs nsweep.csv, fig2.pdf, fig4.pdf. Fig. 2(b) is drawn from data/fig2b_tuning_paths.csv."""
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from scipy.optimize import brentq
import reproduce_final as rf
def design_band_r(mc_meV,N,rfix=None):
    nuc=rf.nu_of_m(mc_meV*1e-3); rr=(lambda a: 0.03*a) if rfix is None else (lambda a: rfix)
    g=lambda a: rf.peak_info(a,rr(a),a,N,nuc,0.3,601)[0]/nuc-1
    l=np.log(1/(2*np.pi*0.03)); a0=0.873*rf.c0*np.sqrt(2*np.pi/l)/(2*np.pi*nuc)
    a=brentq(g,0.7*a0,1.3*a0,xtol=1e-4*a0); f0,pk,Q,I=rf.peak_info(a,rr(a),a,N,nuc); return a,rr(a),f0,pk,Q,I
def Qwidth(a,r,d,N,f0,scale):
    old=rf.RS_SCALE; rf.RS_SCALE=scale
    try:
        f=np.linspace(0.97,1.03,1201)*f0; b=rf.spec_conf(a,r,d,N,f); i=b.argmax()
        res=rf.minimize_scalar(lambda x:-abs(rf.beta_conf(x,a,r,d,N,rf.Rs_p))**2,bracket=(f[i-1],f[i],f[i+1]),tol=1e-10); fp=res.x; P=-res.fun
        gg=lambda x: abs(rf.beta_conf(x,a,r,d,N,rf.Rs_p))**2-P/2
        return fp/(brentq(gg,fp,fp*1.03,xtol=1e-11*fp)-brentq(gg,fp*0.97,fp,xtol=1e-11*fp))
    finally: rf.RS_SCALE=old
if __name__=="__main__":
    bands=[((0.25,0.5),None),((0.5,1),None),((1,2),None),((2,4),None),((4,8),2.5e-6)]
    rows=[]
    for (lo,hi),rfix in bands:
        for N in [3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,23,26,30]:
            a,r_,f0,pk,Q,I=design_band_r(np.sqrt(lo*hi),N,rfix); Qrad=Qwidth(a,r_,a,N,f0,1e-4)
            rows.append(dict(band=f"{lo}-{hi}",N=N,a_um=a*1e6,L_mm=N*a*1e3,beta2=pk,Q_r=Q,I=I,Q_rad=Qrad))
    NS=pd.DataFrame(rows); NS['Q_abs']=1/(1/NS.Q_r-1/NS.Q_rad); NS['x']=NS.Q_abs/NS.Q_rad; NS.to_csv('nsweep.csv',index=False)
    # slab boost, Fig. 4
    def beta_slab(Lol,eps):
        n=np.sqrt(eps+0j); phi=2*np.pi*n*Lol; return (1-np.cos(phi))/(eps*np.cos(phi)-1j*n*np.sin(phi))
    Q=300.; Lg=np.linspace(0.01,8,9000); uopt=(np.pi/(4*Q))**(2/3)
    fig,ax=plt.subplots(figsize=(5.6,4.3)); ax.plot(Lg,abs(beta_slab(Lg,1.0))**2,label='eps=1')
    for u in [0.05,uopt,0.006]: ax.plot(Lg,abs(beta_slab(Lg,u+1j/Q))**2,label=f'u={u:.3f}')
    ax.set_yscale('log'); ax.set_ylim(0.5,4000); ax.set_xlabel('L/lambda'); ax.set_ylabel('beta^2'); ax.legend(); fig.savefig('fig4.pdf')
