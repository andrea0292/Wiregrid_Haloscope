"""Figure scripts for Figs. 2 and 3 (house style). Requires arrays produced by the reach/design scripts."""
def make_fig2_5(path):
    fig,axs=plt.subplots(1,2,figsize=(10.5,4.2))
    allc=list(curves.items())+[((4,8),curves5)]
    for (key,v),cc in zip(allc,cols5):
        axs[0].plot(v['m'],v['d_over_a'],color=cc); axs[1].plot(v['m'],v['b'],color=cc)
        sub=Rp[Rp.band==f'{key[0]}-{key[1]}'] if key!=(4,8) else Rp5
        axs[0].plot(sub.m_fw,sub.d_over_a,'o',mfc='w',mec=cc,mew=1.8,ms=7)
        axs[1].plot(sub.m_fw,sub.b_fw_same,'o',mfc='w',mec=cc,mew=1.8,ms=7)
        Nk=des[key]['N'] if key in des else N5
        axs[0].text(key[0]*1.14,1.55,rf'${Nk}$',color=cc,fontsize=15)
    axs[0].text(0.27,1.40,r'$N=$',fontsize=15)
    for ax in axs:
        ax.set_xscale('log'); ax.set_xlim(0.24,8.4); ax.set_xlabel(r'$m_a$ [meV]')
        ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8'])
    axs[0].set_ylabel(r'$d/a$'); axs[1].set_ylabel(r'$\beta^2$'); axs[1].set_yscale('log'); axs[1].set_ylim(380,8000)
    axs[1].plot([],[],'-',color='k',label=r'impedance sheets'); axs[1].plot([],[],'o',mfc='w',mec='k',mew=1.8,ms=7,label=r'circular wires')
    axs[1].legend(fontsize=14,loc='upper right')
    for a_,l_ in zip(axs,['(a)','(b)']): a_.text(0.04,0.06,l_,transform=a_.transAxes,fontsize=18)
    fig.tight_layout(); fig.savefig(path+'.pdf',bbox_inches='tight'); fig.savefig(path+'.png',dpi=150,bbox_inches='tight'); return fig

def make_fig3_8(path,o_,o2_):
    fig,ax=plt.subplots(figsize=(5.6,4.7))
    mq=np.geomspace(0.25,8,200)
    ax.fill_between(mq,2.04e-13*0.75*mq,2.04e-13*1.92*mq,color='#E8C547',alpha=0.55,lw=0,zorder=1)
    ax.text(0.29,4.6e-14,r'QCD axion',rotation=14,fontsize=15,color='#8a6d00')
    ax.fill_between([0.25,8],[5.8e-11]*2,[1e-8]*2,color='0.78',lw=0,zorder=2); ax.text(0.27,2e-10,r'CAST',fontsize=15,color='0.25')
    ax.axhline(4.7e-11,color='0.35',ls='--',lw=1.6,zorder=3); ax.text(2.0,3.0e-11,r'globular clusters',fontsize=13,color='0.3')
    m_=o_.m.values
    def band(col,color,msk,ls='-',alpha=0.22):
        ax.fill_between(m_[msk],o_[col].values[msk],o2_[col].values[msk],color=color,alpha=alpha,lw=0,zorder=4)
        ax.plot(m_[msk],o_[col].values[msk],color=color,lw=2.6,ls=ls,zorder=5)
    band('lin_ref',C[1],m_<=4.0)
    band('cnt_ref',C[3],m_>=0.8)
    band('target',C[2],m_>=0.8); band('target',C[2],m_<=0.8,ls=(0,(4,2.5)),alpha=0.12)
    ax.text(0.30,2.0e-12,r'$1\,{\rm m}^2$, linear',color=C[1],fontsize=14)
    ax.text(0.85,1.0e-11,r'$1\,{\rm m}^2$, counting',color=C[3],fontsize=14)
    ax.text(1.4,3.0e-14,r'$10\,{\rm m}^2$, counting',color=C[2],fontsize=14)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.25,8); ax.set_ylim(5e-15,1e-9)
    ax.set_xlabel(r'$m_a$ [meV]'); ax.set_ylabel(r'$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]')
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8'])
    sx=ax.secondary_xaxis('top',functions=(lambda m: m*241.8,lambda f: f/241.8)); sx.set_xlabel(r'$\nu$ [THz]',fontsize=17)
    sx.set_xticks([100,200,500,1000,1900]); sx.set_xticklabels(['0.1','0.2','0.5','1','1.9']); sx.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout(); fig.savefig(path+'.pdf',bbox_inches='tight'); fig.savefig(path+'.png',dpi=150,bbox_inches='tight'); return fig


# v6: Fig. 3 with the KID-array reference (kid_ref column) replacing the 1 m^2 counting reference
def make_fig3_9(path,o_,o2_):
    fig,ax=plt.subplots(figsize=(5.6,4.7))
    mq=np.geomspace(0.25,8,200)
    ax.fill_between(mq,2.04e-13*0.75*mq,2.04e-13*1.92*mq,color='#E8C547',alpha=0.55,lw=0,zorder=1)
    ax.text(0.29,4.6e-14,r'QCD axion',rotation=14,fontsize=15,color='#8a6d00')
    ax.fill_between([0.25,8],[5.8e-11]*2,[1e-8]*2,color='0.78',lw=0,zorder=2); ax.text(0.27,2e-10,r'CAST',fontsize=15,color='0.25')
    ax.axhline(4.7e-11,color='0.35',ls='--',lw=1.6,zorder=3); ax.text(2.0,3.0e-11,r'globular clusters',fontsize=13,color='0.3')
    m_=o_.m.values
    def band(col,color,msk,ls='-',alpha=0.22):
        ax.fill_between(m_[msk],o_[col].values[msk],o2_[col].values[msk],color=color,alpha=alpha,lw=0,zorder=4)
        ax.plot(m_[msk],o_[col].values[msk],color=color,lw=2.6,ls=ls,zorder=5)
    band('lin_ref',C[1],m_<=4.0)
    band('kid_ref',C[3],m_>=0.25)
    band('target',C[2],m_>=0.8); band('target',C[2],m_<=0.8,ls=(0,(4,2.5)),alpha=0.12)
    ax.text(0.30,2.0e-12,r'$1\,{\rm m}^2$, linear',color=C[1],fontsize=14)
    ax.text(1.05,1.3e-11,r'$1\,{\rm m}^2$, KID array',color=C[3],fontsize=14)
    ax.text(0.9,3.0e-14,r'$10\,{\rm m}^2$, counters (required)',color=C[2],fontsize=14)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.25,8); ax.set_ylim(5e-15,1e-9)
    ax.set_xlabel(r'$m_a$ [meV]'); ax.set_ylabel(r'$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]')
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8'])
    sx=ax.secondary_xaxis('top',functions=(lambda m: m*241.8,lambda f: f/241.8)); sx.set_xlabel(r'$\nu$ [THz]',fontsize=17)
    sx.set_xticks([100,200,500,1000,1900]); sx.set_xticklabels(['0.1','0.2','0.5','1','1.9']); sx.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout(); fig.savefig(path+'.pdf',bbox_inches='tight'); fig.savefig(path+'.png',dpi=150,bbox_inches='tight'); return fig


# v7: Fig. 3 adds NuSTAR solar bound (AxionLimits limit_data/AxionPhoton/NuSTAR_Sun.txt as nm [meV], ng)
def make_fig3_9(path,o_,o2_):
    fig,ax=plt.subplots(figsize=(5.6,4.7))
    mq=np.geomspace(0.25,8,200)
    ax.fill_between(mq,2.04e-13*0.75*mq,2.04e-13*1.92*mq,color='#E8C547',alpha=0.55,lw=0,zorder=1)
    ax.text(0.29,4.6e-14,r'QCD axion',rotation=14,fontsize=15,color='#8a6d00')
    ax.fill_between([0.25,8],[5.8e-11]*2,[1e-8]*2,color='0.78',lw=0,zorder=2); ax.text(0.27,2e-10,r'CAST',fontsize=15,color='0.25')
    ax.axhline(4.7e-11,color='0.35',ls='--',lw=1.6,zorder=3); ax.text(2.0,3.0e-11,r'globular clusters',fontsize=13,color='0.3')
    ax.plot(nm[sn],ng[sn],color='0.2',ls=':',lw=2.0,zorder=3); ax.text(0.56,1.05e-10,r'NuSTAR',fontsize=11,color='0.2')
    m_=o_.m.values
    def band(col,color,msk,ls='-',alpha=0.22):
        ax.fill_between(m_[msk],o_[col].values[msk],o2_[col].values[msk],color=color,alpha=alpha,lw=0,zorder=4)
        ax.plot(m_[msk],o_[col].values[msk],color=color,lw=2.6,ls=ls,zorder=5)
    band('lin_ref',C[1],m_<=4.0)
    band('kid_ref',C[3],m_>=0.25)
    band('target',C[2],m_>=0.8); band('target',C[2],m_<=0.8,ls=(0,(4,2.5)),alpha=0.12)
    ax.text(0.30,2.0e-12,r'$1\,{\rm m}^2$, linear',color=C[1],fontsize=14)
    ax.text(1.05,1.3e-11,r'$1\,{\rm m}^2$, KID array',color=C[3],fontsize=14)
    ax.text(0.9,3.0e-14,r'$10\,{\rm m}^2$, counters (required)',color=C[2],fontsize=14)
    ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlim(0.25,8); ax.set_ylim(5e-15,1e-9)
    ax.set_xlabel(r'$m_a$ [meV]'); ax.set_ylabel(r'$|g_{a\gamma\gamma}|$ [GeV$^{-1}$]')
    ax.xaxis.set_minor_formatter(NullFormatter()); ax.set_xticks([0.25,0.5,1,2,4,8]); ax.set_xticklabels(['0.25','0.5','1','2','4','8'])
    sx=ax.secondary_xaxis('top',functions=(lambda m: m*241.8,lambda f: f/241.8)); sx.set_xlabel(r'$\nu$ [THz]',fontsize=17)
    sx.set_xticks([100,200,500,1000,1900]); sx.set_xticklabels(['0.1','0.2','0.5','1','1.9']); sx.xaxis.set_minor_formatter(NullFormatter())
    fig.tight_layout(); fig.savefig(path+'.pdf',bbox_inches='tight'); fig.savefig(path+'.png',dpi=150,bbox_inches='tight'); return fig
