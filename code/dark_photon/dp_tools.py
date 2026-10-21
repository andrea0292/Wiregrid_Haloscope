"""Dark-photon reach of the wire-grid stack (Letter v13, Fig. 4 and Appendix D).
Usage: python dp_tools.py   (reads ../../data/reach_0p25_8_fw.csv, writes ../../data/dp_reach_check.csv)
"""
import numpy as np, pandas as pd
rng=np.random.default_rng(1)
h=6.62607015e-34; eV=1.602176634e-19
B0ev=1953.5                       # 10 T in eV^2 (natural units)
rhoc=0.4e9*eV*1e6*2.99792458e8    # rho_X c in W/m^2 for rho=0.4 GeV/cm^3 (P/A = kappa chi^2 rho c)
Om=np.pi*(7.338e-4)**2            # solid angle of the halo-velocity spot (sigma_v-based), ~1.69e-6 sr
def hnu(m): return m*1e-3*eV

def chi_from_g(g,m_meV,kappa=1/3):
    """axion reach g [GeV^-1] at B0=10 T -> dark-photon reach chi; chi cos(theta)=g B0/m."""
    return g*1e-9*B0ev/(m_meV*1e-3)/np.sqrt(kappa)

def pol_factors(Tw_days, n=200000, split=True):
    # Earth axis = z. Polarization n uniform on sphere. Wire A: along z (parallel to axis): c1=nz^2 constant.
    # Wire B: perpendicular to axis, rotating at sidereal rate: cos^2 = (nx cos(phi)+ny sin(phi))^2 = s2 cos^2(phi-phi_n)
    nz=rng.uniform(-1,1,n); s2=1-nz**2; ph0=rng.uniform(0,2*np.pi,n)
    w=2*np.pi*Tw_days/0.99727  # sidereal
    # average of cos^2(ph0+t) over window length w
    avg=0.5+ (np.sin(2*(ph0+w))-np.sin(2*ph0))/(4*w)
    c2=s2*avg; c1=nz**2
    keff_split=np.sqrt(0.5*(c1**2+c2**2)); keff_perp=c2
    return np.percentile(keff_split,5), np.percentile(keff_perp,5), keff_split.mean(), np.percentile(c1,5)


def dish_chi(m,A=10,T_yr=5,eta=0.5,kap=2/3,d=4e-3,NEP=None,Nd=None):
    """ideal broadband dish, SNR=5; Nd detectors (default: number of halo-spot modes, >=1)."""
    lam=1.23984e-3/m; Nd=np.maximum(1,A*Om/lam**2) if Nd is None else Nd
    dd=d if NEP is None else NEP**2/(2*hnu(m)**2)
    T=T_yr*3.156e7
    return (25*Nd*dd*hnu(m)**2/(T*eta**2*kap**2*A**2*rhoc**2))**0.25


def dish_het(m,A=10,T_yr=5,eta=0.5,kap=2/3,BIF=8e9):
    lam=1.23984e-3/m; Nd=np.maximum(1,A*Om/lam**2); nu=m*241.799e9
    tm=T_yr*3.156e7*BIF/(1.93e12-0.0605e12)
    Nsys=10*6.62607015e-34*nu; dnu=1e-6*nu
    return (25*Nd*Nsys**2*dnu/(tm*eta**2*kap**2*A**2*rhoc**2))**0.25



def proto_chi(NEP=3e-19,t=86400.,A=0.01,eta=0.313,beta2=1648.,kappa=1/3,snr=5):
    Pmin=snr*NEP/np.sqrt(2*t)
    return np.sqrt(Pmin/(eta*A*kappa*rhoc*beta2))

if __name__=="__main__":
    A=pd.read_csv("../../data/reach_0p25_8_fw.csv")
    B=pd.read_csv("../../data/reach_0p25_8_fw_Rs2.csv")
    ntile=np.where(A.m<4.0,1000,4000)   # 10 cm tiles below 4 meV, 5 cm above
    out=pd.DataFrame({"m_meV":A.m,"het":np.where(A.m<=4.0001,chi_from_g(A.lin10,A.m),np.nan),"kid":chi_from_g(A.kid10,A.m),"cnt":chi_from_g(A.target,A.m),
                      "het_Rs2":np.where(A.m<=4.0001,chi_from_g(B.lin10,A.m),np.nan),"kid_Rs2":chi_from_g(B.kid10,A.m),"cnt_Rs2":chi_from_g(B.target,A.m),
                      "cnt_fixed":chi_from_g(A.target,A.m,kappa=0.299),
                      # one tile with its own KID, one month per octave: g ~ (A a_t T)^(-1/4)
                      "proto_kid":chi_from_g(A.kid10,A.m)*(ntile*12)**0.25,
                      "dish_cnt":dish_chi(A.m),"dish_kid":dish_chi(A.m,NEP=3e-19)})
    out.to_csv("../../data/dp_reach_check.csv",index=False)
    for m in [0.25,0.5,1,2,4,8]:
        r=out.iloc[(out.m_meV-m).abs().argmin()]; print(m, "%.2e %.2e %.2e"%(r.het,r.kid,r.cnt))
    print("constant", chi_from_g(1e-12,1.0))
    print("prototype", proto_chi())
    for Tw in [0.9,1.4,2.2,3.6,4.3]: print("pol",Tw, np.round(pol_factors(Tw),3))
    print("BREAD KID 1000 d check (rescaled to BREAD rho=0.45)", np.sqrt(0.4/0.45)*dish_chi(np.array([1.0]),A=10,T_yr=1000/365.25,eta=0.5,kap=2/3,NEP=2e-19,Nd=1)[0], '(BREAD: 3.1e-14)')
