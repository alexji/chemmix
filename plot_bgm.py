import numpy as np
import pylab as plt

from optparse import OptionParser
import backgroundmodel as bgm
import util

def plot_one_bgm(axarr,envname,**kwargs):
    tarr = bgm.gettarr(envname); tmax=np.max(tarr)
    Mhalo,zvir,vturb,lturb,nSN,trecovery,Dt,uSN,E51 = bgm.envparams(envname)
    th,rho,Vmix,Vmixdot = util.load_Vmix(envname,get_thrho=True,getVmixdot=True)
    ((ax0,ax1),(ax2,ax3)) = axarr
    ax0.plot(tarr,rho(tarr),**kwargs)
    ax0.plot([100,100],[10**-27,10**-20],'c:')
    ax0.plot([100+trecovery,100+trecovery],[10**-27,10**-20],'c:')
    ax0.set_ylabel(r'$\rho$')
    #ax0.set_xscale('log'); 
    ax0.set_yscale('log')
    ax0.set_xlabel('t')
    ax0.set_xlim((.1,np.max(tarr)))
    ax0.set_ylim((10**-27,10**-20))

    ax1.plot(tarr,Vmix(tarr),**kwargs)
    ax1.set_ylabel(r'$V_{\rm mix}$')
    ax1.set_xlabel('t')
    ax1.set_ylim((0,ax1.get_ylim()[-1]))

    ax3.plot(tarr,Vmixdot(tarr),**kwargs)
    ax3.set_ylabel(r'$\dot{V}_{\rm mix}$')
    ax3.set_xlabel('t')

    mmg = util.load_mmixgrid(envname)
    mmg = np.log10(mmg) #plot log
    if 'minihalo' in envname:
        levels = np.arange(1,5.5,.2)
    if 'atomiccoolinghalo' in envname:
        levels = np.arange(1,7.5,.2)
    CS = ax2.contour(mmg,levels,extent=(0,tmax,0,tmax),aspect=1)
    ax2.clabel(CS,levels,inline=1,fontsize=6,fmt="%1.1f")
    ax2.set_ylabel('t'); ax2.set_xlabel(r'$\tau$')
    #ax2.text(.5*tmax,.3*tmax,envname,fontsize=10)
    ax2.text(.75*tmax,.15*tmax,r'$M_{\rm mix}(t,\tau)$',fontsize=10)

if __name__=="__main__":
    #envlist = ['minihalo_lowE','atomiccoolinghalo','minihalo','atomiccoolinghalo_lowDt','minihalo_lowDt']
    envlist = ['minihalo_z17','minihalo_z17_lowE','minihalo_z20','minihalo_z20_lowE']
    for envname in envlist:
        fig,axarr = plt.subplots(2,2)
        fig.subplots_adjust(hspace=.25,wspace=.3)
        plot_one_bgm(axarr,envname)
        plt.savefig('PLOTS/bgm_'+envname+'.png',bbox_inches='tight')
