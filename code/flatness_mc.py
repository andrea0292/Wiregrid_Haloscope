"""Grid-flatness Monte Carlo: average the complex impedance-sheet response over random
independent out-of-plane displacements of each grid (patch-to-patch), then compare the
scan response I = int beta^4 dnu/nu with the ideal stack. Uses functions from reproduce_final.py."""
import numpy as np
exec(open('reproduce_final.py').read().split('if __name__')[0])
def beta_pos(nu,a,r,xs,Rs_fun=Rs_p):
    """confined drive, sheets at arbitrary positions xs (ascending), mirror at 0."""
    N=len(xs); k=2*np.pi*nu/c0; Zg=grid_impedance(nu,a,r,Rs_fun(nu))[0]; w=2*np.pi*nu; g=-1j*w*mu0/Zg
    M=np.zeros((2*N+1,2*N+1),complex); rhs=np.zeros(2*N+1,complex)
    M[0,0]=1; M[0,1]=1; rhs[0]=1; row=1
    for j in range(N):
        x=xs[j]; e=np.exp(1j*k*x); ei=1/e; pout=-1.0 if j<N-1 else 0.0
        M[row,2*j]=e; M[row,2*j+1]=ei
        if j<N-1: M[row,2*(j+1)]=-e; M[row,2*(j+1)+1]=-ei
        else: M[row,2*N]=-e
        rhs[row]=-(-1.0-pout); row+=1
        if j<N-1: M[row,2*(j+1)]=1j*k*e; M[row,2*(j+1)+1]=-1j*k*ei
        else: M[row,2*N]=1j*k*e
        M[row,2*j]+=-1j*k*e-g*e; M[row,2*j+1]+=1j*k*ei-g*ei; rhs[row]=g*(-1.0); row+=1
    return np.linalg.solve(M,rhs)[2*N]

def flat_I(D,sig,f0,Qr,nreal=300,seed=1):
    rng=np.random.default_rng(seed); a,r,N=D['a'],D['r'],D['N']; xs0=a*np.arange(1,N+1)
    fs=f0*(1+np.linspace(-6,6,161)/Qr); eps=rng.normal(0,sig,(nreal,N))
    B=np.array([[beta_pos(f,a,r,xs0+e) for f in fs] for e in eps])  # nreal x nf
    bmean=B.mean(0); b0=np.array([beta_pos(f,a,r,xs0) for f in fs])
    I=lambda b: np.trapezoid(np.abs(b)**4,fs)/f0
    return I(bmean)/I(b0), np.max(np.abs(bmean)**2)/np.max(np.abs(b0)**2)

if __name__=="__main__":
    # example: 4-8 meV cartridge (N=8, a=63.75 um, 2r=5 um), centre resonance computed internally
    a,r,N=63.74776988601645e-6,2.5e-6,8
    D=dict(N=N,a=a,r=r)
    f0,pk,Q,I=peak_info(a,r,a,N,nu_of_m(5.657e-3))
    D.update(b=pk,Q=Q,I=I); Qr=Q_exact(D,a,f0,pk)[0]
    for s in [0.25e-6,0.5e-6,1e-6,2e-6]: print(s*1e6, flat_I(D,s,f0,Qr,nreal=400))
