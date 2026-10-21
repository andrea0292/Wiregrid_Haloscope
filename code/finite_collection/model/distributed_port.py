"""Finite Drude plasma with a passive, spatially distributed coherent port.

PEC rectangle in x,y, wires and E along z, uniform in z. The port profile is
a sum of Gaussian pickup columns coupled to ONE ideal coherent receiver.
This is a homogenized network benchmark, not a wire-resolved antenna design.
All modal source overlaps, losses and port backaction are retained.
"""
import numpy as np
from scipy.constants import c,h,e,mu_0
from scipy.special import erf
from scipy.optimize import differential_evolution, minimize, minimize_scalar,brentq
from alpha_comparison import F_SQUARE, resistance

def gaussian_sine(m,L,centers,w):
    q=np.pi*m/L; result=np.zeros(len(m))
    for xc in centers:
        v=np.sqrt(np.pi)*w/2*np.exp(1j*q*xc-(q*w/2)**2)
        v*=erf((L-xc)/w-1j*q*w/2)-erf(-xc/w-1j*q*w/2)
        result+=np.sqrt(2/L)*v.imag
    return result

def frequency_grid(Q,span=.15,step=.05):
    # Dimensionless detuning x=2Q(nu/nu0-1); all poles have intrinsic width ~2.
    extent=2*Q*span
    xx=np.unique(np.r_[np.linspace(-extent,-20,max(2,int(np.ceil((extent-20)/(.2*step/.05)))+1)),np.arange(-20,20+step/2,step),np.linspace(20,extent,max(2,int(np.ceil((extent-20)/(.2*step/.05)))+1))])
    xx=xx[xx>-1.8*Q]
    return 1+xx/(2*Q)

class PlasmaPort:
    def __init__(self,mass,L,W=.1,xi=.03,pickup_width_pitches=2.,basis_scale=1.,span=.15,step=.05):
        self.mass,self.L,self.W=mass,L,W
        self.nu=mass*1e-3*e/h; self.k=2*np.pi*self.nu/c
        self.pp2=1-((np.pi/L)**2+(np.pi/W)**2)/self.k**2
        self.base_pp2=self.pp2
        if self.pp2<=0: raise ValueError('Vacuum cutoff exceeds target frequency')
        ell=np.log(1/(2*np.pi*xi))+F_SQUARE
        self.pitch=2*np.pi/(self.k*np.sqrt(self.pp2)*np.sqrt(2*np.pi*ell))
        self.radius=xi*self.pitch
        self.Gamma=resistance(self.nu)/(mu_0*self.radius*ell)
        self.Q=2*np.pi*self.nu/(self.Gamma*self.pp2)
        self.w=pickup_width_pitches*self.pitch
        mx=max(31,int(np.ceil(6*L/(np.pi*self.w))))
        my=max(31,int(np.ceil(6*W/(np.pi*self.w))))
        self.mx=np.arange(1,int(mx*basis_scale)+1,2)
        self.my=np.arange(1,int(my*basis_scale)+1,2)
        mm,nn=np.meshgrid(self.mx,self.my,indexing='ij')
        self.q2=((mm*np.pi/L)**2+(nn*np.pi/W)**2).ravel()/self.k**2
        self.f=(8/(np.pi**2*mm*nn)).ravel()
        self.z=frequency_grid(self.Q,span,step)
        self.span,self.step=span,step

    def diag(self,z):
        z=np.asarray(z)
        gamma=self.Gamma/(2*np.pi*self.nu)*z**(2/3)
        return self.Q*(self.q2[None,:]-z[:,None]**2+self.pp2*z[:,None]/(z[:,None]+1j*gamma[:,None]))

    def profile(self,nx,ny,uniform=False):
        if uniform: return self.f/np.linalg.norm(self.f)
        ux=gaussian_sine(self.mx,self.L,(np.arange(nx)+.5)*self.L/nx,self.w)
        uy=gaussian_sine(self.my,self.W,(np.arange(ny)+.5)*self.W/ny,self.w)
        u=np.outer(ux,uy).ravel()
        return u/np.linalg.norm(u)

    def greens(self,u,z=None):
        if z is None:z=self.z
        G=np.zeros(len(z),complex);H=G.copy()
        for start in range(0,len(z),128):
            sl=slice(start,start+128);zz=z[sl]
            inv=1/self.diag(zz)
            G[sl]=inv@abs(u)**2
            H[sl]=zz**2*(inv@(u.conj()*self.f))
        return G,H

    def spectrum(self,G,H,alpha,zeta,z=None):
        if z is None:z=self.z
        t=H/(1+z*(alpha-1j*zeta)*G)
        return self.k*self.L*self.Q*zeta*abs(t)**2

    def optimize(self,nx=1,ny=1,uniform=False,reactive=True):
        self.pp2=self.base_pp2
        u=self.profile(nx,ny,uniform);G,H=self.greens(u)
        G0=self.greens(u,np.array([1.]))[0][0]
        inv=1/G0; am=-inv.real; zm=-inv.imag
        def unpack(v): return (am+v[1]*zm if reactive else 0.,zm*np.exp(v[0]))
        def integral(v):
            alpha,zeta=unpack(v)
            b=self.spectrum(G,H,alpha,zeta)
            return float(np.trapezoid(b*b,np.log(self.z)))
        bounds=[(-5.,5.),(-12.,12.)] if reactive else [(-5.,5.)]
        opt=differential_evolution(lambda v:-integral(v),bounds,popsize=10,maxiter=120,tol=2e-7,seed=71,polish=True)
        alpha,zeta=unpack(opt.x)
        # Retune the material plasma frequency, keeping wire damping and the
        # chosen passive load fixed, so the loaded peak is at the target mass.
        def at_target(pp):
            self.pp2=pp;gg,hh=self.greens(u,np.array([1.]))
            return self.spectrum(gg,hh,alpha,zeta,np.array([1.]))[0]
        extent=.04
        while True:
            ps=np.linspace(self.base_pp2-extent,self.base_pp2+extent,int(4000*extent)+1)
            powers=np.array([at_target(p) for p in ps]);jj=int(np.argmax(powers))
            if jj not in (0,len(ps)-1):break
            extent*=2
            if extent>.16:raise RuntimeError('Retuning bracket insufficient')
        pp=minimize_scalar(lambda p:-at_target(p),bounds=(ps[jj-1],ps[jj+1]),method='bounded',options={'xatol':1e-13}).x
        def derivative(p):
            self.pp2=p;zz=np.array([1-1e-6,1+1e-6]);gg,hh=self.greens(u,zz)
            bb=self.spectrum(gg,hh,alpha,zeta,zz)
            return (bb[1]-bb[0])/(2e-6)
        self.pp2=brentq(derivative,pp-.001,pp+.001,xtol=1e-13)
        G,H=self.greens(u);b=self.spectrum(G,H,alpha,zeta)
        j=int(np.argmax(b));zpk=float(self.z[j])
        # Exact source-work identity at the numerical peak.
        D=self.diag(np.array([zpk]))[0];beta=zpk*(alpha-1j*zeta)
        t=H[j]/(1+beta*G[j]);field=(zpk*zpk*self.f-beta*u*t)/D
        loss=float(np.sum(-D.imag*abs(field)**2));extracted=float(zpk*zeta*abs(t)**2)
        work=float(np.imag(np.vdot(zpk*zpk*self.f,field)))
        # Receiver contribution to field equation has nonnegative damping.
        source_norm=float(np.sum(self.f**2))
        meta=dict(mass_meV=self.mass,L_mm=self.L*1e3,W_mm=self.W*1e3,
            pitch_um=self.pitch*1e6,radius_um=self.radius*1e6,Q_unloaded=2*np.pi*self.nu/(self.Gamma*self.pp2),Q_reference=self.Q,
            plasma_ratio_squared=self.pp2,base_plasma_ratio_squared=self.base_pp2,
            plasma_frequency_retune_fraction=float(np.sqrt(self.pp2/self.base_pp2)-1),
            pickup_width_um=self.w*1e6,pickups=nx*ny,nx=nx,ny=ny,
            uniform_profile=uniform,reactive_tuning=reactive,alpha=alpha,zeta=zeta,
            I=integral(opt.x),peak_boost=float(b[j]),peak_over_nominal=zpk,
            source_norm=source_norm,basis_modes=len(self.f),basis_shape=[len(self.mx),len(self.my)],
            energy_balance_relative=(work-loss-extracted)/work,
            port_fraction_of_dissipation=extracted/(loss+extracted),
            optimizer_success=bool(opt.success),parameter_boundary=bool(any(abs(v-a)<.03 or abs(v-bb)<.03 for v,(a,bb) in zip(opt.x,bounds))),**spectrum_summary(self.z,b))
        return meta,b,u

def grid_shapes(K):
    return [(n,K//n) for n in range(1,K+1) if K%n==0]

def spectrum_summary(z,b):
    j=int(np.argmax(b));half=b[j]/2
    lo=j
    while lo>0 and b[lo]>half:lo-=1
    hi=j
    while hi<len(z)-1 and b[hi]>half:hi+=1
    left=float(np.interp(half,[b[lo],b[lo+1]],[z[lo],z[lo+1]]))
    right=float(np.interp(half,[b[hi],b[hi-1]],[z[hi],z[hi-1]]))
    xx=np.r_[left,z[lo+1:hi],right];yy=np.r_[half,b[lo+1:hi],half]
    return dict(I_FWHM=float(np.trapezoid(yy**2,np.log(xx))),Q_response=float(z[j]/(right-left)),FWHM_bounds=[left,right])

def single_mode_validation():
    x=np.linspace(-200,200,200001)
    from scipy.optimize import minimize_scalar
    def integral(zeta):return np.trapezoid((zeta/(x*x+(1+zeta)**2))**2,x)
    z=minimize_scalar(lambda z:-integral(z),bounds=(.1,6),method='bounded').x
    exact=np.pi*2**2/(2*3**3)
    return dict(optimal_zeta=z,expected_zeta=2.,integral_at_zeta2=integral(2),exact=exact,relative_error=integral(2)/exact-1)
