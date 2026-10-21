"""Finite flanged rectangular impedance-grid enclosure, exp(-i omega t).
x normal, y across wires, z along wires. PEC back/side/end walls. N anisotropic
sheets Jz=Ez/Zg at x=j*d; finite W by H aperture in an infinite ideal PEC flange.
This vector mode-matching model resolves macroscopic fields, not circular wires.
Exterior admittance uses Galerkin scalar-Green integrals, including evanescence.
Default parity: l odd, n even. Both tangential polarizations retained for n>0.
"""
import numpy as np
from scipy.constants import c,mu_0,epsilon_0
from scipy.fft import fft2,ifft2,next_fast_len
from scipy.special import roots_legendre,roots_hermite
Z0=np.sqrt(mu_0/epsilon_0)

def basis_fourier(m,q,L,kind):
 q=np.asarray(q)
 if m==0:return np.sqrt(L)*np.sinc(q*L/(2*np.pi)) if kind=='cos' else np.zeros(q.shape,complex)
 minus=np.exp(.5j*m*np.pi)*np.sinc(q*L/(2*np.pi)-m/2)
 plus=np.exp(-.5j*m*np.pi)*np.sinc(q*L/(2*np.pi)+m/2)
 return np.sqrt(L/2)*(minus+plus if kind=='cos' else (minus-plus)/1j)

def scalar_kernel(k,ny,nz,hy,hz,order=4):
 """Double-cell-averaged Green, divided by one cell area; Duffy self integral."""
 iy=np.arange(-(ny-1),ny);iz=np.arange(-(nz-1),nz)
 y=iy[:,None]*hy;z=iz[None,:]*hz;g=np.zeros((len(iy),len(iz)),complex)
 u,w=roots_legendre(order);u=(u+1)/2;w=w/2
 nodes=np.r_[-u[::-1],u];weights=np.r_[w[::-1],w]*(1-abs(nodes))
 for a,wa in zip(nodes,weights):
  for b,wb in zip(nodes,weights):
   R=np.sqrt((y+a*hy)**2+(z+b*hz)**2);R=np.where(R==0,1.,R)
   g+=wa*wb*np.exp(1j*k*R)/(4*np.pi*R)*hy*hz
 u,w=roots_legendre(40);u=(u+1)/2;w=w/2
 t=u[:,None];s=u[None,:];ww=w[:,None]*w[None,:];val=0j
 for aa,bb in [(hy,hz),(hz,hy)]:
  rr=np.sqrt(aa*aa+(bb*s)**2)
  val+=np.sum(ww*(1-t)*(1-t*s)*np.exp(1j*k*t*rr)/rr)*hy*hz/np.pi
 g[ny-1,nz-1]=val
 u,w=roots_legendre(28);u=(u+1)/2;w=w/2
 nodes=np.r_[-u[::-1],u];weights=np.r_[w[::-1],w]*(1-abs(nodes))
 for dy in [-1,0,1]:
  for dz in [-1,0,1]:
   if dy==dz==0:continue
   R=np.sqrt(((dy+nodes[:,None])*hy)**2+((dz+nodes[None,:])*hz)**2)
   val=np.sum(weights[:,None]*weights[None,:]*np.exp(1j*k*R)/(4*np.pi*R))*hy*hz
   g[dy+ny-1,dz+nz-1]=val
 shape=(next_fast_len(2*ny-1),next_fast_len(2*nz-1));emb=np.zeros(shape,complex)
 emb[np.ix_(iy%shape[0],iz%shape[1])]=g
 return fft2(emb,workers=4)

class Aperture:
 def __init__(self,W,H,lmax=15,nmax=4,ny=160,nz=None):
  self.W,self.H=W,H;self.ny=ny;self.nz=nz if nz is not None else round(ny*H/W)
  self.hy=W/self.ny;self.hz=H/self.nz;self.groups=[];self.modes=[]
  for l in range(1,lmax+1,2):
   for n in range(0,nmax+1,2):
    idx=[]
    if n>0:idx.append(len(self.modes));self.modes.append((l,n,0))
    idx.append(len(self.modes));self.modes.append((l,n,1));self.groups.append((l,n,idx))
  self.M=len(self.modes);self.y=(np.arange(self.ny)+.5)*self.hy;self.z=(np.arange(self.nz)+.5)*self.hz
  self.F=[];self.C=[];self.drive=np.zeros(self.M)
  for i,(l,n,pol) in enumerate(self.modes):
   sy=np.sqrt(2/W)*np.sin(l*np.pi*self.y/W);cy=np.sqrt(2/W)*np.cos(l*np.pi*self.y/W)
   cn=np.sqrt((1 if n==0 else 2)/H)*np.cos(n*np.pi*self.z/H)
   sn=np.sqrt((1 if n==0 else 2)/H)*np.sin(n*np.pi*self.z/H)
   f=cy[:,None]*sn[None,:] if pol==0 else sy[:,None]*cn[None,:]
   curl=-(n*np.pi/H)*cy[:,None]*cn[None,:] if pol==0 else (l*np.pi/W)*cy[:,None]*cn[None,:]
   self.F.append(f.ravel());self.C.append(curl.ravel())
   if pol==1 and n==0:self.drive[i]=2*np.sqrt(2*W*H)/(l*np.pi)
  self.F=np.array(self.F);self.C=np.array(self.C);self.pol=np.array([x[2] for x in self.modes]);self.area=W*H

 def radiation_admittance(self,nu,kernel_order=4):
  if self.ny <= max(m[0] for m in self.modes) or self.nz <= max(m[1] for m in self.modes):
   raise ValueError('Aperture quadrature must resolve the retained basis; interior-only objects may use a small grid.')
  k=2*np.pi*nu/c;K=scalar_kernel(k,self.ny,self.nz,self.hy,self.hz,kernel_order)
  nfft=K.shape;dA=self.hy*self.hz;GF=np.empty((self.M,self.M),complex);GC=np.empty_like(GF)
  for j in range(self.M):
   f=self.F[j].reshape(self.ny,self.nz);cur=self.C[j].reshape(self.ny,self.nz)
   conv=ifft2(fft2(f,s=nfft,workers=4)*K,workers=4)[:self.ny,:self.nz].ravel()
   convc=ifft2(fft2(cur,s=nfft,workers=4)*K,workers=4)[:self.ny,:self.nz].ravel()
   GF[:,j]=(self.F@conv)*dA*(self.pol==self.pol[j]);GC[:,j]=(self.C@convc)*dA
  self.Yrad=(-2j/k)*(k*k*GF-GC)
  self.reciprocity_residual=float(np.linalg.norm(self.Yrad-self.Yrad.T)/np.linalg.norm(self.Yrad))
  self.Yrad=(self.Yrad+self.Yrad.T)/2;self.nu_rad=nu
  return self.Yrad

 def interior(self,nu,a,r,N,d=None,rs_factor=1.8,contact_ohm=0):
  if d is None:d=a
  k=2*np.pi*nu/c;omega=2*np.pi*nu;Rs=rs_factor*2.9e-3*(nu/3.9e9)**(2/3)
  Ls=mu_0*a/(2*np.pi)*np.log(a/(2*np.pi*r));Rg=Rs*a/(2*np.pi*r)
  Rg+=2*contact_ohm*a/self.H # uniform equivalent contact-loss proxy, not resolved contacts
  Zg=Rg-1j*omega*Ls;ygrid=Z0/Zg
  Yin=np.zeros((self.M,self.M),complex);src=np.zeros(self.M,complex);records=[]
  for l,n,idx in self.groups:
   q=np.array([l*np.pi/self.W,n*np.pi/self.H]);v=np.array([-q[1],q[0]])
   A=np.eye(2)-np.outer(q,q)/k**2;B=np.eye(2)-np.outer(v,v)/k**2;kx=np.sqrt(k*k-q@q+0j)
   ct=np.cos(kx*d);ss=np.sin(kx*d)/kx
   T=np.block([[ct*np.eye(2),1j*k*ss*A],[1j*k*ss*B,ct*np.eye(2)]])
   S=np.eye(4,dtype=complex);S[3,1]=-ygrid;sv=np.zeros(2,complex)
   for ii in idx:sv[self.modes[ii][2]]=self.drive[ii]
   Ep=-np.linalg.solve(B,sv);pp=np.r_[Ep,0j,0j];f=(np.eye(4)-T)@pp
   TT=np.eye(4,dtype=complex);ff=np.zeros(4,complex)
   for j in range(N):TT=S@T@TT;ff=S@(T@ff+f)
   U=TT[:2,2:];V=TT[2:,2:];yy=V@np.linalg.inv(U);bb=ff[2:]-yy@ff[:2]
   pol=[self.modes[ii][2] for ii in idx]
   Yin[np.ix_(idx,idx)]=yy[np.ix_(pol,pol)];src[idx]=bb[pol]
   records.append((idx,pol,T,S,f,TT,ff,Ep,sv,A,B,ygrid))
  return Yin,src,records,Zg

 def solve(self,nu,a,r,N,d=None,rs_factor=1.8,contact_ohm=0,receiver_w=None,frozen_radiation=False):
  if not hasattr(self,'Yrad') or (not frozen_radiation and abs(nu/self.nu_rad-1)>1e-12):self.radiation_admittance(nu)
  Yin,src,rec,Zg=self.interior(nu,a,r,N,d,rs_factor,contact_ohm)
  D=self.Yrad-Yin;E=np.linalg.solve(D,src);Pr=float(np.real(np.vdot(E,self.Yrad@E)))
  absorption,source_work,wall_info=self.budget(nu,E,rec,N,a if d is None else d)
  result=dict(radiated_boost=Pr/self.area,absorbed_boost=absorption/self.area,source_work_boost=source_work/self.area,
   balance_relative=(source_work-Pr-absorption)/max(abs(source_work),1e-100),
   linear_residual=float(np.linalg.norm(D@E-src)/np.linalg.norm(src)),E=E,Yin=Yin,src=src,records=rec,wall_info=wall_info)
  if receiver_w:
   v,norm=self.gaussian_port(nu,receiver_w);p=abs(np.vdot(v,E))**2/norm
   Epass=np.linalg.solve(D,2*v);apass=self.passive_absorption(nu,Epass,rec,N,a if d is None else d)
   result.update(gaussian_boost=float(p/self.area),gaussian_fraction=float(p/Pr),gaussian_absorptivity=float(apass/norm),receiver_w=receiver_w)
  self.last=result;return result

 def _fields(self,E,rec,N):
  fields=[]
  for idx,pol,T,S,f,TT,ff,Ep,sv,A,B,ygrid in rec:
   ev=np.zeros(2,complex);ev[pol]=E[idx];h0=np.linalg.solve(TT[:2,2:],ev-ff[:2]);state=np.r_[0j,0j,h0];sheets=[];left=[]
   for j in range(N):left.append(state.copy());state=T@state+f;sheets.append(state[:2].copy());state=S@state
   fields.append((np.array(left),np.array(sheets),Ep,sv,A,B,ygrid,h0))
  return fields

 def budget(self,nu,E,rec,N,d):
  k=2*np.pi*nu/c;absorb=0.;work=0.;wall_back=0.;x,w=roots_legendre(18);x=(x+1)*d/2;w=w*d/2
  for group,field in zip(self.groups,self._fields(E,rec,N)):
   l,n,idx=group;left,es,Ep,sv,A,B,yg,h0=field;absorb+=yg.real*np.sum(abs(es)**2);wall_back+=np.sum(abs(h0)**2)
   q2=(l*np.pi/self.W)**2+(n*np.pi/self.H)**2;kx=np.sqrt(k*k-q2+0j)
   for state in left:
    e0,h0=state[:2],state[2:]
    val=Ep[None,:]+np.cos(kx*x)[:,None]*(e0-Ep)[None,:]+1j*k*(np.sin(kx*x)/kx)[:,None]*(A@h0)[None,:]
    work+=k*np.imag(np.sum(w*(val@sv.conj())))
  return float(absorb),float(work),dict(back_mirror_h_norm=float(wall_back))

 def passive_absorption(self,nu,E,rec,N,d):
  ab=0.
  for idx,pol,T,S,f,TT,ff,Ep,sv,A,B,yg in rec:
   ev=np.zeros(2,complex);ev[pol]=E[idx];h0=np.linalg.solve(TT[:2,2:],ev);state=np.r_[0j,0j,h0]
   for j in range(N):state=T@state;ab+=yg.real*np.sum(abs(state[:2])**2);state=S@state
  return float(ab)

 def fourier_fields(self,qy,qz,E):
  qy,qz=np.broadcast_arrays(qy,qz);out=np.zeros((2,)+qy.shape,complex);cachey={};cachez={}
  for coef,(l,n,pol) in zip(E,self.modes):
   ky='cos' if pol==0 else 'sin';kz='sin' if pol==0 else 'cos'
   if (l,ky) not in cachey:cachey[l,ky]=basis_fourier(l,qy,self.W,ky)
   if (n,kz) not in cachez:cachez[n,kz]=basis_fourier(n,qz,self.H,kz)
   out[pol]+=coef*cachey[l,ky]*cachez[n,kz]
  return out

 def gaussian_port(self,nu,w,nquad=64):
  k=2*np.pi*nu/c;u,wt=roots_hermite(nquad);qy=2*u[:,None]/w;qz=2*u[None,:]/w;weight=wt[:,None]*wt[None,:]
  qy,qz=np.broadcast_arrays(qy,qz);kx=np.sqrt(k*k-qy*qy-qz*qz+0j);prop=qy*qy+qz*qz<k*k
  yz=np.where(prop,qy*qz/(k*kx),0.);zz=np.where(prop,(k*k-qy*qy)/(k*kx),0.);v=np.zeros(self.M,complex)
  for j,(l,n,pol) in enumerate(self.modes):
   ff=basis_fourier(l,qy,self.W,'cos' if pol==0 else 'sin')*basis_fourier(n,qz,self.H,'sin' if pol==0 else 'cos')
   v[j]=np.sum(weight*ff.conj()*(yz if pol==0 else zz))/np.pi
  qy=np.sqrt(2)*u[:,None]/w;qz=np.sqrt(2)*u[None,:]/w;kx=np.sqrt(k*k-qy*qy-qz*qz+0j)
  yy=np.where(qy*qy+qz*qz<k*k,(k*k-qy*qy)/(k*kx),0.);norm=float(np.real(w*w/2*np.sum(weight*yy)))
  return v,norm

 def farfield_power(self,nu,E,ntheta=240,nphi=512):
  k=2*np.pi*nu/c;u,w=roots_legendre(ntheta);th=(u+1)*np.pi/4;wt=w*np.pi/4;ph=2*np.pi*np.arange(nphi)/nphi
  st=np.sin(th)[:,None];ct=np.cos(th)[:,None];qy=k*st*np.cos(ph)[None,:];qz=k*st*np.sin(ph)[None,:]
  F=self.fourier_fields(qy,qz,E);para=np.cos(ph)[None,:]*F[0]+np.sin(ph)[None,:]*F[1]
  integrand=k*k*st*(ct*ct*np.sum(abs(F)**2,axis=0)+st*st*abs(para)**2)/(2*np.pi)**2
  return float(np.sum(wt[:,None]*integrand)*2*np.pi/nphi)

 def fast_response(self,nu,a,r,N,d=None,rs_factor=1.8,Yrad=None,port=None):
  """Spectrum using caller-supplied radiative admittance and receiving mode.
  No implicit reuse at a different frequency; caller controls interpolation.
  """
  if Yrad is None:
   if not hasattr(self,'Yrad') or abs(nu/self.nu_rad-1)>1e-12:self.radiation_admittance(nu)
   Yrad=self.Yrad
  Yin,src,rec,Zg=self.interior(nu,a,r,N,d,rs_factor)
  E=np.linalg.solve(Yrad-Yin,src)
  rad=float(np.real(np.vdot(E,Yrad@E)))/self.area
  collected=None if port is None else float(abs(np.vdot(port[0],E))**2/(port[1]*self.area))
  return rad,collected,E,Yin,src

 def localized_port(self,nu,waist,nquad=120):
  """L2-normalized Gaussian aperture-E constraint, PEC on orthogonal modes.
  Ideal modal transformer benchmark, not a resolved physical antenna.
  """
  t,w=roots_legendre(nquad);y=(t+1)*self.W/2;z=(t+1)*self.H/2
  gy=np.exp(-((y-self.W/2)/waist)**2);gz=np.exp(-((z-self.H/2)/waist)**2)
  p=np.zeros(self.M)
  for j,(l,n,pol) in enumerate(self.modes):
   if pol==0:continue
   iy=np.sum(w*gy*np.sqrt(2/self.W)*np.sin(l*np.pi*y/self.W))*self.W/2
   iz=np.sum(w*gz*np.sqrt((1 if n==0 else 2)/self.H)*np.cos(n*np.pi*z/self.H))*self.H/2
   p[j]=iy*iz
  return p/np.linalg.norm(p)

 def receiving_profile(self,nu,port,rec,N,d,Yin):
  """Depth-integrated reciprocal Ez coefficients; passive illumination.
  x phase neglected here, error scale (L/ell_a)^2 for halo applications.
  """
  ev=np.linalg.solve(self.Yrad-Yin,2*port[0]);coef=np.zeros(self.M,complex)
  k=2*np.pi*nu/c
  for l,n,idx in self.groups:
   record=next(rr for rr in rec if rr[0]==idx)
   _,pol,T,S,f,TT,ff,Ep,sv,A,B,yg=record
   es=np.zeros(2,complex);es[pol]=ev[idx];h0=np.linalg.solve(TT[:2,2:],es);state=np.r_[0j,0j,h0]
   kx=np.sqrt(k*k-(l*np.pi/self.W)**2-(n*np.pi/self.H)**2+0j)
   inte=np.zeros(2,complex)
   for j in range(N):
    inte+=np.sin(kx*d)/kx*state[:2]+1j*k*(1-np.cos(kx*d))/kx**2*(A@state[2:])
    state=S@T@state
   for ii in idx:
    if self.modes[ii][2]==1:coef[ii]=inte[1]
  return coef

 def perturbative_wall_loss(self,nu,E,rec,N,d,wall_Rs):
  """Ohmic power evaluated on PEC fields: not self-consistent lossy-wall solve."""
  k=2*np.pi*nu/c;xx,ww=roots_legendre(18);xx=(xx+1)*d/2;ww=ww*d/2
  fields=self._fields(E,rec,N);back=sum(np.sum(abs(f[-1])**2) for f in fields)
  sides=0.;ends=0.
  for j in range(N):
   hn={};hl={};xn={};xl={}
   for (l,n,idx),field in zip(self.groups,fields):
    left,es,Ep,sv,A,B,yg,h0=field;e0=left[j,:2];h0=left[j,2:];kx=np.sqrt(k*k-(l*np.pi/self.W)**2-(n*np.pi/self.H)**2+0j)
    ee=Ep[None,:]+np.cos(kx*xx)[:,None]*(e0-Ep)[None,:]+1j*k*(np.sin(kx*xx)/kx)[:,None]*(A@h0)[None,:]
    hh=np.cos(kx*xx)[:,None]*h0[None,:]+1j*k*(np.sin(kx*xx)/kx)[:,None]*(B@(e0-Ep))[None,:]
    hx=((l*np.pi/self.W)*ee[:,1]-(n*np.pi/self.H)*ee[:,0])/(1j*k)
    sy=np.sqrt(2/self.W);sz=np.sqrt((1 if n==0 else 2)/self.H)
    hn[n]=hn.get(n,0)+sy*hh[:,0];xn[n]=xn.get(n,0)+sy*hx
    hl[l]=hl.get(l,0)-sz*hh[:,1];xl[l]=xl.get(l,0)+sz*hx
   sides+=2*sum(np.sum(ww*abs(x)**2) for x in [*hn.values(),*xn.values()])
   ends+=2*sum(np.sum(ww*abs(x)**2) for x in [*hl.values(),*xl.values()])
  return {'back':float(wall_Rs/Z0*back/self.area),'sides':float(wall_Rs/Z0*sides/self.area),'ends':float(wall_Rs/Z0*ends/self.area)}
