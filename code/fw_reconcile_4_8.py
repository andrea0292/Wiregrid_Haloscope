
import numpy as np, json, sys
from scipy.optimize import brentq, minimize_scalar
exec(open('fullwave.py').read())
exec(open('reproduce_final.py').read().split('if __name__')[0])
def fw_conf(st,f,phase):
    Zs=Rs_p(f)*(1-1j*phase); k=2*np.pi*f/c0; L=st.z[-1]
    be=st.solve(f,Zs,'axion'); r0=st.solve(f,Zs,'passive')
    return be-(1-r0*np.exp(2j*k*L))*np.exp(-1j*k*L)/2
def peak(fun,fg,Qg,n=25):
    fs=np.linspace(fg*(1-4/Qg),fg*(1+4/Qg),n); vs=[abs(fun(x))**2 for x in fs]; i=int(np.argmax(vs))
    if i in (0,n-1): return peak(fun,fs[i],Qg/2,n)
    res=minimize_scalar(lambda f:-abs(fun(f))**2,bracket=(fs[i-1],fs[i],fs[i+1]),tol=1e-9); return res.x,-res.fun
def spec_I(fun,fg,Qg):
    f0,pk=peak(fun,fg,Qg); h=lambda f: abs(fun(f))**2-pk/2; w=f0/Qg
    fl=brentq(h,f0-3*w,f0,xtol=1e-9*f0); fr=brentq(h,f0,f0+3*w,xtol=1e-9*f0); Q=f0/(fr-fl)
    y=np.r_[np.linspace(-12,-2,40,endpoint=False),np.linspace(-2,2,161,endpoint=False),np.linspace(2,12,41)]*(f0/Q)
    b2=np.array([abs(fun(f0+x))**2 for x in y]); I=np.trapezoid(b2**2,f0+y)/f0
    return dict(f0=f0,pk=pk,Q=Q,I=I)
def run(args):
    band,setting,m,N,a,da,phase,r=args; d=da*a; ft=nu_of_m(m*1e-3)
    s=spec_I(lambda f: beta_conf(f,a,r,d,N,Rs_p),ft,300)
    st=Stack(a,r,d*np.arange(1,N+1),M=3); w=spec_I(lambda f: fw_conf(st,f,phase),s['f0']*0.98,s['Q'])
    # retune spacing so the wire resonance sits at the sheet resonance
    dd=d*(w['f0']/s['f0'])**(1/0.64)
    for it in range(6):
        st2=Stack(a,r,dd*np.arange(1,N+1),M=3); f0,_=peak(lambda f: fw_conf(st2,f,phase),s['f0'],s['Q'])
        if abs(f0/s['f0']-1)<2e-5: break
        dd=dd*(f0/s['f0'])**(1/0.64)
    wr=spec_I(lambda f: fw_conf(st2,f,phase),s['f0'],s['Q'])
    return dict(band=band,setting=setting,m=m,N=N,a=a,d_over_a=da,phase=phase,sheet=s,fw_same=w,fw_retuned=wr,d_retuned_over_d=dd/d,
        sqrtI_same=float(np.sqrt(w['I']/s['I'])),sqrtI_retuned=float(np.sqrt(wr['I']/s['I'])),b_same=float(w['pk']/s['pk']),b_retuned=float(wr['pk']/s['pk']),
        fshift_same=float(w['f0']/s['f0']-1),fres_retuned=float(wr['f0']/s['f0']-1))
if __name__=='__main__':
    import os
    from multiprocessing import Pool
    a5=float(sys.argv[1]); r5=float(sys.argv[2]); N5=int(sys.argv[3])
    jobs=[]
    for setting,m,da in [('lower',4.0,float(sys.argv[4])),('centre',5.656854,1.0),('upper',8.0,float(sys.argv[5]))]:
        for ph in [1.0,np.sqrt(3)]: jobs.append(('4-8',setting,m,N5,a5,da,float(ph),r5))
    with Pool(6) as p: out=p.map(run,jobs)
    json.dump(out,open('fw_reconcile_4_8.json','w'),indent=1,default=float); print('done',len(out))
