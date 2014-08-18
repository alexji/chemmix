import numpy as np
import pylab as plt

from optparse import OptionParser
import backgroundmodel as bgm
import util

def plot_one_sfr(axarr,envname,sfrname,**kwargs):
    tarr = bgm.gettarr(envname)
    uII,uIII,muII,muIII,wII,wIII = util.load_sfr(envname,sfrname)
    ((ax0,ax1),(ax2,ax3)) = axarr
    ax0.plot(tarr,uII(tarr),**kwargs)
    ax0.set_title('uII')
    ax0.set_ylim((0,ax0.get_ylim()[-1]))
    ax1.plot(tarr,uIII(tarr),**kwargs)
    ax1.set_ylim((0,ax1.get_ylim()[-1]))
    ax1.set_title('uIII')
    ax2.plot(tarr,muII(tarr),**kwargs)
    ax2.set_title('muII')
    ax3.plot(tarr,muIII(tarr),**kwargs)
    ax3.set_title('muIII')
    #ax0.set_yscale('log')
    #ax1.set_yscale('log'); ax1.set_ylim((10**-4,10**1.5))
    #ax2.set_yscale('log')
    #ax3.set_yscale('log')

if __name__=="__main__":
    #parser = OptionParser()
    #options,args = parser.parse_args()
    #envname = args[0]
    #sfrlist = args[1:]

    #envname='atomiccoolinghalo'
    #sfrlist = ['burstachsimTS500','burstachsloTS500']
    #labellist = ['ACH burst 500','ACH burstlo 500']
        #sfrlist = ['fixTS500','10xTS500']

    #envname='atomiccoolinghalo_lowDt'
    #sfrlist = ['fixTS500','10xTS500']
    #labellist = ['ACH fix 500','ACH 10x 500']

    #fig,axarr = plt.subplots(2,2,sharex=True)
    #for sfrname in sfrlist:
    #    plot_one_sfr(axarr,envname,sfrname)
    #plt.legend(labellist,loc='best',fontsize=10)
    #plt.savefig("PLOTS/sfr_"+envname+".png",bbox_inches='tight')

    envlist = ['minihalo','minihalo_lowE']
    sfrlist = ['betalogitTS200']
    sfrlist = ['betalogitVmaxTS200']
#    sfrlist = ['mhsimTS150']
    labellist = ['MHbl200','MHlEbl200']#['MH','MH E49']
    labellist = ['MHVmbl200','MHlEVmbl200']#['MH','MH E49']
#        #sfrlist = ['fixTS150']
    fig,axarr = plt.subplots(2,2,sharex=True)
    for envname in envlist:
        for sfrname in sfrlist:
            plot_one_sfr(axarr,envname,sfrname)
    plt.legend(labellist,loc='best',fontsize=10)
    #plt.savefig("PLOTS/sfr_minihalo.png",bbox_inches='tight')
#    plt.savefig("PLOTS/sfr_minihalo_bl.png",bbox_inches='tight')
    plt.savefig("PLOTS/sfr_minihalo_Vm_bl.png",bbox_inches='tight')
    #plt.show()

    #envname='minihalo_lowDt'
    #sfrlist = ['fixTS150']
    #labellist = ['MH fix 150']

