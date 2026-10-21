"""Regenerates comparator_R_vs_envelope.png and R_vs_Nel.json (requires house-plot-style apply_house_style in namespace)."""
import json, sys, os, numpy as np, matplotlib.pyplot as plt
sys.path.append(os.path.dirname(os.path.abspath(__file__))); import comparator_model as cm
try: apply_house_style()
except NameError: pass
Nels=[10,30,100,300,1e3,3e3,1e4,3e4,1e5,1e6,1e7]; R0=cm.equal_volume(2066,313)["R0"]; curves={}
for env in cm.envelopes():
    for scn,sc in cm.SCEN.items():
        cfgs=cm.stack_configs(env,sc,1.0); Rl=[]
        for Ne in Nels:
            pb=cm.port_best(env,sc,1.0,100,Ne); best=None
            for cfg in cfgs:
                sb=cm.stack_best(cfg,0.68,100,Ne)
                if sb is None: continue
                R=R0*(cfg["A"]*cfg["d"])**2*sb[0]/pb["fig"]; best=R if best is None or R>best else best
            Rl.append(best)
        curves[(env["name"],scn)]=Rl
json.dump({f"{k[0]}|{k[1]}":v for k,v in curves.items()}|{"Nel":Nels},open("R_vs_Nel.json","w"),indent=1)
box=json.load(open("box.json"))
fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13.5,5.6)); cols=plt.rcParams['axes.prop_cycle'].by_key()['color']
for j,N in enumerate([1,4,16,64,256]):
    pts=sorted([(r["W"],r["eta_ext"]) for r in box if r["kind"]=="lattice" and r["N"]==N])
    W_,e_=zip(*pts); ax1.plot(W_,e_,'o-',color=cols[j],label=f"$N_{{\\rm port}}={N}$",ms=5)
ax1.plot([r["W"] for r in box if r["kind"]=="uniform"],[r["eta_ext"] for r in box if r["kind"]=="uniform"],'k--',lw=1.8,label="uniform loading")
ax1.axhline(1,color='0.5',lw=1,ls=':'); ax1.axvline(cm.ell_em0*2e3,color='0.6',lw=1,ls='-.'); ax1.text(11.0,1.9,r"$2\ell_{\rm em}$",fontsize=14,color='0.4')
ax1.set_xscale('log'); ax1.set_yscale('log'); ax1.set_ylim(2e-3,3); ax1.set_xlim(4,120)
ax1.set_xlabel(r"transverse size $W$ [mm]"); ax1.set_ylabel(r"extraction efficiency $\eta_{\rm ext}$")
ax1.legend(fontsize=12.5,loc='lower left'); ax1.set_title("(a) port-coupled wire medium, 273.6 GHz",fontsize=15)
envn=[e["name"] for e in cm.envelopes()]; lab=dict(zip(envn,["10 cm cube (solenoid)","10 cm cube (dipole)","MADMAX-like dipole","BREAD-like solenoid"]))
for j,en in enumerate(envn):
    n,o,p=[np.array(curves[(en,s)],float) for s in ("port-nominal","port-optimistic","port-pessimistic")]
    ax2.fill_between(Nels,np.maximum(o,1e-4),p,color=cols[j],alpha=0.18,lw=0); ax2.plot(Nels,np.maximum(n,1e-4),'-',color=cols[j],label=lab[en])
ax2.axhline(2,color='k',lw=1.6,ls='--'); ax2.text(1.3e5,2.6,"gate $R=2$",fontsize=14); ax2.axhline(1,color='0.5',lw=1,ls=':')
ax2.set_xscale('log'); ax2.set_yscale('log'); ax2.set_ylim(1e-3,1e5); ax2.set_xlim(10,1e7)
ax2.set_xlabel(r"max. coupling elements $N_{\rm el}$ (both designs)"); ax2.set_ylabel(r"scan-rate ratio $R=\dot\nu_s/\dot\nu_p$")
ax2.set_title(r"(b) fixed envelope, $\eta_s=0.68$, $N_{\rm rx}\leq100$",fontsize=15); ax2.legend(fontsize=12.5,loc='upper right')
fig.tight_layout(); fig.savefig("comparator_R_vs_envelope.png",dpi=150)
