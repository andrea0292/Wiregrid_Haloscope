from pathlib import Path
import sys,json,hashlib
import numpy as np
from scipy.constants import c,h,e,mu_0,epsilon_0
from numpy.polynomial.hermite import hermgauss
P=Path(__file__).resolve().parent;sys.path.insert(0,str(P/'model'))
from distributed_port import PlasmaPort,spectrum_summary
from aperture_model import Aperture,basis_fourier
from alpha_comparison import resistance

def sine_source(m,L,q):
    arg=q[:,None]*L/(2*np.pi);mm=m[None,:]
    return np.sqrt(2)/(2j)*(np.exp(1j*mm*np.pi/2)*np.sinc(arg+mm/2)-np.exp(-1j*mm*np.pi/2)*np.sinc(arg-mm/2))

def port_halo(obj,row,order=16):
    obj.pp2=row['plasma_ratio_squared']
    t,w=hermgauss(order);w=w/np.sqrt(np.pi);sigma=220e3/np.sqrt(2)
    vx=np.sqrt(2)*sigma*t;vy=232e3+np.sqrt(2)*sigma*t;vz=vx
    fx=sine_source(obj.mx,obj.L,obj.k*vx/c)
    fy=sine_source(obj.my,obj.W,obj.k*vy/c)
    chiz=float(np.sum(w*np.sinc(obj.k*vz/c*.1/(2*np.pi))**2))
    u=obj.profile(row.get('nx',1),row.get('ny',1),row.get('uniform_profile',False))
    vals=[]
    for z in [row['FWHM_bounds'][0],row['peak_over_nominal'],row['FWHM_bounds'][1]]:
        D=obj.diag(np.array([z]))[0]
        coeff=(u/D).reshape(len(obj.mx),len(obj.my))
        amp=fx@coeff@fy.T;norm=abs(np.sum(u*obj.f/D))**2
        vals.append(float(np.sum(w[:,None]*w[None,:]*abs(amp)**2)/norm*chiz))
    return dict(chi_peak=vals[1],chi_FWHM_edges=[vals[0],vals[2]])

def rad_halo(m,r,order=16):
    cachepath=P/'results/radiative_halo_cache.json'
    cache=json.loads(cachepath.read_text()) if cachepath.exists() else {}
    key=hashlib.sha256(json.dumps(r,sort_keys=True).encode()).hexdigest()+f':{order}'
    if key in cache:return cache[key]
    nu=m*1e-3*e/h;g=r['best'];v=r['finite']
    o=Aperture(.1,.1,lmax=23,nmax=4,ny=v['ny']);Y=o.radiation_admittance(nu)
    f=nu*v['peak_over_nominal'];Yin,src,rec,Zg=o.interior(f,g['a'],g['r'],g['N'],d=g['d'])
    port=o.gaussian_port(nu,v['waist_m'])
    co=o.receiving_profile(f,port,rec,g['N'],g['d'],Yin)
    t,w=hermgauss(order);w=w/np.sqrt(np.pi);sigma=220e3/np.sqrt(2)
    k=2*np.pi*nu/c;qy=k*(232e3+np.sqrt(2)*sigma*t)/c;qz=k*np.sqrt(2)*sigma*t/c
    amps=np.zeros((order,order),complex)
    for coeff,(l,n,pol) in zip(co,o.modes):
        if pol==1:amps+=coeff*np.outer(basis_fourier(l,qy,.1,'sin'),basis_fourier(n,qz,.1,'cos'))
    norm=abs(np.sum(co*o.drive))**2
    chi=float(np.sum(w[:,None]*w[None,:]*abs(amps)**2)/norm)
    cache[key]=chi;cachepath.write_text(json.dumps(cache,indent=2))
    return chi

def guide_loss(m,ratio=1.7,rs_scale=1.):
    nu=m*1e-3*e/h;lam=c/nu;aa=ratio*lam/2;bb=aa/2
    rs=rs_scale*resistance(nu);Z0=np.sqrt(mu_0/epsilon_0)
    alpha=rs/(Z0*bb)*(1+(2*bb/aa)/ratio**2)/np.sqrt(1-1/ratio**2)
    return dict(a_mm=aa*1e3,b_mm=bb*1e3,surface_resistance_ohm=rs,
        attenuation_amplitude_Np_m=alpha,power_loss_dB_m=float(20/np.log(10)*alpha))

out={'halo_corrected':[],'loss_thresholds':[],
    'halo_assumptions':dict(sigma_km_s=220/np.sqrt(2),wind_y_km_s=232,wire_length_m=.1,
    method='Reciprocal profile at the peak; multiplication of scan integrals by chi_peak squared. Port variation across FWHM is recorded. Radiative source phase across the thin depth is neglected.')}
for m in [.5,1.,2.,4.]:
    r=json.loads((P/'results'/f'radiative_{m:g}.json').read_text())
    chi=rad_halo(m,r);a=np.load(P/'results'/f'radiative_spectrum_{m:g}.npz')
    rad_count=spectrum_summary(a['z'],a['gaussian'])['I_FWHM']
    for label in ['same_active_volume','30mm_available_depth']:
        path=P/'results'/f'ports_{m:g}_{label}.json'
        if not path.exists():continue
        d=json.loads(path.read_text());obj=PlasmaPort(m,d['L_mm']*1e-3)
        for family in d['rows']+[d['continuous_profile']]:
            alternatives=[]
            for row in family.get('candidate_geometries',[family]):
                hp=port_halo(obj,row)
                R=r['finite']['I_gaussian']/row['I'];Rc=rad_count/row['I_FWHM']
                rr=dict(mass_meV=m,geometry=label,pickups=row['pickups'],nx=row['nx'],ny=row['ny'],uniform_profile=row['uniform_profile'],
                    R=R,R_count=Rc,chi_radiative=chi,**hp,
                    R_halo=R*(chi/hp['chi_peak'])**2,R_count_halo=Rc*(chi/hp['chi_peak'])**2)
                alternatives.append(rr)
            out['halo_corrected'].append(min(alternatives,key=lambda q:q['R_halo']))
        # Strongest tested one-output port, including the continuous profile.
        rows=[x for x in out['halo_corrected'] if x['mass_meV']==m and x['geometry']==label]
        strong=min(rows,key=lambda x:x['R_halo']);R=strong['R_halo'];Rc=min(x['R_count_halo'] for x in rows)
        gg=guide_loss(m);eta_opt=.8
        threshold=-10*np.log10(eta_opt*np.sqrt(R))
        row=dict(mass_meV=m,geometry=label,R_best_port=R,R_count_best_port=Rc,
            optics_power_efficiency=eta_opt,port_loss_for_parity_dB=float(threshold),
            guide_length_for_parity_cm=float(100*threshold/gg['power_loss_dB_m']),guide=gg,
            conditional_R_20cm=float(R*(eta_opt/10**(-gg['power_loss_dB_m']*.2/10))**2),
            port_receiver_count=1)
        out['loss_thresholds'].append(row)
    print('MASS',m,'halo rad',chi,'threshold',out['loss_thresholds'][-1],flush=True)
(P/'results/resource_assessment.json').write_text(json.dumps(out,indent=2))
print(json.dumps(out['loss_thresholds'],indent=2),flush=True)
