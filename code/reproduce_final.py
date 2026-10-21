"""Reproduces the cartridge designs (Table I), tuning curves (Fig. 2, sheet model) and reach curves (Fig. 3)
of "A Radiative Wire-Grid Plasma Haloscope for Millielectronvolt Axion Dark Matter".
Impedance-sheet model with confined axion drive; Rs = RS_SCALE * 83.6 mOhm (nu/250 GHz)^(2/3); r/a = 0.03.
Usage: python reproduce_final.py [--rs2]   -> writes design_check.csv, reach_check.csv
Runtime: a few minutes on a laptop."""
import sys, numpy as np, pandas as pd
from scipy.optimize import brentq, minimize_scalar
c0=299792458.0; mu0=1.2566370614359173e-06; hP=6.62607015e-34; h=hP; kB=1.380649e-23; eV=1.602176634e-19; yr=3.156e7
RS_SCALE=2.0 if '--rs2' in sys.argv else 1.0
def nu_of_m(m_eV):
    return m_eV * eV / h

def grid_impedance(nu, a, r, Rs):
    w = 2 * np.pi * nu
    Ls = mu0 * a / (2 * np.pi) * np.log(a / (2 * np.pi * r))
    Rg = Rs * a / (2 * np.pi * r)
    return -1j * w * Ls + Rg, Ls, Rg
def beta_conf(nu,a,r,d,N,Rs_fun=lambda nu: 83.6e-3*(nu/250e9)**(2/3)):
    """confined drive: B only for 0<x<L=N d; sheets at x_j=j d, j=1..N; mirror at 0. E0=1."""
    k=2*np.pi*nu/c0; Zg=grid_impedance(nu,a,r,Rs_fun(nu))[0]; w=2*np.pi*nu
    xs=d*np.arange(1,N+1)
    # unknowns: A_j,B_j for regions j=0..N-1 (inside), A_N outgoing  -> 2N+1
    M=np.zeros((2*N+1,2*N+1),complex); rhs=np.zeros(2*N+1,complex)
    def col(j): return 2*j
    # mirror: -1 + A0 + B0 = 0
    M[0,0]=1; M[0,1]=1; rhs[0]=1
    row=1
    for j in range(N):
        x=xs[j]; e=np.exp(1j*k*x); ei=1/e
        pin=-1.0; pout=-1.0 if j<N-1 else 0.0
        # left region j: E_L = pin + A_j e + B_j ei ; right region j+1
        # continuity: E_L - E_R = 0
        M[row,2*j]=e; M[row,2*j+1]=ei
        if j<N-1: M[row,2*(j+1)]=-e; M[row,2*(j+1)+1]=-ei
        else: M[row,2*N]=-e
        rhs[row]=-(pin-pout); row+=1
        # jump: E'_R - E'_L = -(i w mu0/Zg) E(x) ; E(x)=E_L
        g=-1j*w*mu0/Zg
        if j<N-1: M[row,2*(j+1)]=1j*k*e; M[row,2*(j+1)+1]=-1j*k*ei
        else: M[row,2*N]=1j*k*e
        M[row,2*j]+=-1j*k*e; M[row,2*j+1]+=1j*k*ei
        M[row,2*j]-=g*e; M[row,2*j+1]-=g*ei; rhs[row]=g*pin; row+=1
    sol=np.linalg.solve(M,rhs); return sol[2*N]

def spec_conf(a,r,d,N,nus): return np.array([abs(beta_conf(f,a,r,d,N,Rs_p))**2 for f in nus])

def peak_fast(a,r,d,N,nu_g,span=0.12,n=241):
    f=np.linspace(1-span,1+span,n)*nu_g; b=spec_conf(a,r,d,N,f); i=b.argmax()
    res=minimize_scalar(lambda x:-abs(beta_conf(x,a,r,d,N,Rs_p))**2,bracket=(f[i-1],f[i],f[i+1]),tol=1e-9); return res.x,-res.fun

def peak_info(a,r,d,N,nu_g,span=0.25,n=801):
    f=np.linspace(1-span,1+span,n)*nu_g; b=spec_conf(a,r,d,N,f); i=b.argmax()
    # refine
    f2=np.linspace(f[max(i-2,0)],f[min(i+2,n-1)],401); b2=spec_conf(a,r,d,N,f2); j=b2.argmax(); f0=f2[j]; pk=b2[j]
    # FWHM and scan integral on a fine grid around peak
    ff=np.linspace(f0*(1-0.04),f0*(1+0.04),4001); bb=spec_conf(a,r,d,N,ff)
    above=ff[bb>=pk/2]; fw=above.max()-above.min()
    I=np.trapezoid(bb**2,ff)/f0
    return f0,pk,f0/fw,I

def design_band(mc_meV,N):
    nuc=nu_of_m(mc_meV*1e-3)
    def g(a): return peak_info(a,0.03*a,a,N,nuc*1.0,0.3,601)[0]/nuc-1
    # initial guess: omega_c ~0.87 omega_p, omega_p = c sqrt(2pi/l)/a
    l=np.log(1/(2*np.pi*0.03)); a0=0.873*c0*np.sqrt(2*np.pi/l)/(2*np.pi*nuc)
    a=brentq(g,0.85*a0,1.15*a0,xtol=1e-4*a0)
    f0,pk,Q,I=peak_info(a,0.03*a,a,N,nuc)
    return a,f0,pk,Q,I

def Q_exact(D,dv,f0,pk):
    g=lambda f: abs(beta_conf(f,D['a'],D['r'],dv,D['N'],Rs_p))**2-pk/2
    fl=brentq(g,f0*(1-0.02),f0,xtol=1e-9*f0); fr=brentq(g,f0,f0*(1+0.02),xtol=1e-9*f0); return f0/(fr-fl),fl,fr

def build_traj(DES, edges):
    tr_={}
    for (lo,hi),D in DES.items():
        lo_g=lo/1.015 if lo==edges[0] else lo; hi_g=hi*1.015 if hi==edges[-1] else hi
        ms=np.geomspace(lo_g,hi_g,15); ds=[]; mc=np.sqrt(lo*hi)
        for m in ms:
            dg=D['a']*(mc/m)**(1/0.64) if not ds else ds[-1]*(ms[len(ds)-1]/m)**(1/0.64)
            tgt=nu_of_m(m*1e-3); g=lambda dv: peak_fast(D['a'],D['r'],dv,D['N'],tgt*(dg/dv)**0.64)[0]/tgt-1
            ds.append(brentq(g,dg*0.85,dg*1.15,xtol=1e-6*dg))
        tr_[(lo,hi)]=dict(m=ms,d=np.array(ds),lo=lo_g,hi=hi_g)
    return tr_

def J_both(m_a, DES, TRJ):
    nua=nu_of_m(m_a*1e-3); out=[]
    for key,tr in TRJ.items():
        if not (tr['lo']/1.05 <= m_a <= tr['hi']*1.05): continue
        D=DES[key]; ms=np.clip(m_a,tr['lo'],tr['hi']); dv=np.exp(np.interp(np.log(ms),np.log(tr['m']),np.log(tr['d'])))
        fpk,pk=peak_fast(D['a'],D['r'],dv,D['N'],nu_of_m(ms*1e-3),0.03,121)
        w=fpk/D['Q']; y=np.linspace(-8*w,8*w,1601); b=spec_conf(D['a'],D['r'],dv,D['N'],fpk+y)
        half=y[b>=pk/2]; Qr=fpk/(half.max()-half.min())
        nlo,nhi=nu_of_m(tr['lo']*1e-3),nu_of_m(tr['hi']*1e-3); mask=((nua-y)>=nlo)&((nua-y)<=nhi)
        Jl=np.trapezoid(np.where(mask,b**2,0),y)/nua; mc=mask&(np.abs(y)<fpk/(2*Qr)); Jc=np.trapezoid(np.where(mc,b**2,0),y)/nua
        out.append(dict(Jl=Jl,Jc=Jc,Qr=Qr,fpk=fpk,b=pk))
    return out

def chi_their(m_meV,w=0.05,L=0.10,wind=232e3,sigv=220e3/np.sqrt(2)):
    k=2*np.pi*nu_of_m(m_meV*1e-3)/c0
    xx=np.linspace(-L/2,L/2,2001); u=np.exp(-xx**2/w**2); n0=abs(np.trapezoid(u,xx))**2
    tg_,wg_=np.polynomial.hermite.hermgauss(40)
    def c1(vm):
        v=vm+np.sqrt(2)*sigv*tg_; return np.sum(wg_*np.array([abs(np.trapezoid(u*np.exp(1j*k*vv/c0*xx),xx))**2 for vv in v]))/np.sqrt(np.pi)/n0
    return c1(wind)*c1(0.0)

def reach_from(m_a,terms,mode,A_tot,n_t,T_live,Dd=None,T=None,eta_ap=0.65,chi=None):
    nua=nu_of_m(m_a*1e-3); rho=T_live*yr/(4*np.log(2)); chi=chi_their(m_a) if chi is None else chi
    P1=7.33e-26*0.5*eta_ap*chi*(A_tot/n_t)/m_a**2
    if mode=='lin':
        s=n_t*rho*sum(t['Jl'] for t in terms)*P1**2/((10*hP*nua)**2*1e-6*nua)
    else:
        s=n_t*rho*sum(t['Jc']/(Dd+0.5*(t['fpk']/t['Qr'])/np.expm1(hP*t['fpk']/(kB*T))) for t in terms)*(P1/(hP*nua))**2
    return 1e-12*(25/s)**0.25

Rs_p=lambda nu: RS_SCALE*83.6e-3*(nu/250e9)**(2/3)
if __name__=="__main__":
    bands=[(0.25,0.5),(0.5,1),(1,2),(2,4)]; Ns=[15,13,11,9]; des={}; rows=[]
    for (lo,hi),N in zip(bands,Ns):
        a,f0,pk,Q,I=design_band(np.sqrt(lo*hi),N)
        des[(lo,hi)]=dict(N=N,a=a,r=0.03*a,b=pk,Q=Q,I=I)
        Qr=Q_exact(des[(lo,hi)],a,f0,pk)[0]  # exact half-maximum width (Table I)
        rows.append(dict(band=f"{lo}-{hi}",N=N,a_um=a*1e6,two_r_um=0.06*a*1e6,L_mm=N*a*1e3,beta2=pk,Q_r=Qr,Q_resp_grid=Q))
    pd.DataFrame(rows).to_csv('design_check.csv',index=False); print(pd.DataFrame(rows))
    traj=build_traj(des,[0.25,0.5,1,2,4])
    ms=[0.25,0.5,0.707107,0.802923,1.0,1.414214,2.0,2.828427,4.0]; out=[]
    for m in ms:
        t=J_both(m,des,traj)
        out.append(dict(m_meV=m,lin_ref=reach_from(m,t,'lin',1.0,100,1.0),
            cnt_ref=reach_from(m,t,'cnt',1.0,100,1.0,Dd=4.0,T=0.5),
            kid_ref=reach_from(m,t,'cnt',1.0,100,1.0,Dd=(3e-19)**2/(2*(hP*nu_of_m(m*1e-3))**2),T=0.5),  # KID array NEP 3e-19 W/rtHz as counter with D=NEP^2/2(h nu)^2
            target=reach_from(m,t,'cnt',10.0,1000,4.0,Dd=0.004,T=0.1)))
    df=pd.DataFrame(out); df.to_csv('reach_check.csv',index=False); print(df)
