import numpy as np
import pylab as plt
import h5py

from scipy import stats
import karlsson
from tophat import TopHat


if __name__=="__main__":
    RHOP2 = False

    Mhalo = 10**8; fb = 0.1551; Mmax = Mhalo*fb
    #filenames = ['atomiccoolinghalo','minihalo']
    #filenames = ['atomiccoolinghalo','minihalo4']
    #filenames = ['atomiccoolinghalo','k08']
    #filenames = ['atomiccoolinghalo','atomiccoolinghalo4']
    #filenames = ['atomiccoolinghalo','atomiccoolinghalolate']
    #filenames = ['atomiccoolinghalo','atomiccoolinghalolate_lowvturb']
    #filenames = ['atomiccoolinghalo','atomiccoolinghalo_lowmass']
    #filenames = ['atomiccoolinghalolate_lowvturb','atomiccoolinghalo_lowmass']
    filenames = ['k08','k084']
    lefttaillist = []; righttaillist = []
    q1list = []; medlist = []; q3list = []
    karr = np.arange(1,21)
    offset = [0,.1]
    colors = ['k','r']

    fig = plt.figure(); ax = fig.gca()
    for i,filename in enumerate(filenames):
        Mplot = np.load(karlsson.get_Mplot_filename(filename,rhop2=RHOP2))
        minarr = []; maxarr = []
        meanarr = []
        lefttailarr = []; righttailarr = []
        q1arr = []; medarr = []; q3arr = []
        for k in karr:
            fMk = np.load(karlsson.get_fMk_filename(filename,k,rhop2=RHOP2))
            distr = stats.rv_discrete(name=str(k),values=(Mplot,fMk))
            minarr.append(distr.ppf(0))
            lefttailarr.append(distr.ppf(.05))
            q1arr.append(distr.ppf(.25))
            meanarr.append(distr.mean())
            medarr.append(distr.median())
            q3arr.append(distr.ppf(.75))
            righttailarr.append(distr.ppf(.95))
            maxarr.append(distr.ppf(1))
        minarr = np.array(minarr)
        lefttailarr  = np.array(lefttailarr)
        q1arr        = np.array(q1arr)
        meanarr = np.array(meanarr)
        medarr       = np.array(medarr)
        q3arr        = np.array(q3arr)
        righttailarr = np.array(righttailarr)
        maxarr = np.array(maxarr)
        #ax.plot(karr+offset[i],minarr,colors[i]+':')
        ax.plot(karr+offset[i],lefttailarr,colors[i]+'--')
        ax.plot(karr+2*offset[i],meanarr,colors[i]+'s')#,label=filename+' mean')
        #ax.plot(karr+offset[i],q1arr,colors[i]+':')
        ax.errorbar(karr+offset[i],medarr,yerr=[medarr-q1arr,q3arr-medarr],fmt=colors[i]+'o-',lw=1.8,label=filename+' med')
        #ax.plot(karr+offset[i],q3arr,colors[i]+':')
        ax.plot(karr+offset[i],righttailarr,colors[i]+'--')
        #ax.plot(karr+offset[i],maxarr,colors[i]+':')
        lefttaillist.append(lefttailarr); righttaillist.append(righttailarr)
        q1list.append(q1arr); medlist.append(medarr); q3list.append(q3arr)
    ax.set_xlim((0,21)); ax.set_xlabel(r'$k$')
    ax.set_yscale('log'); ax.set_ylabel(r'M$_{\rm mix}$ (q1,med,q3) $M_\odot$')
    ax.set_ylim((10**4,10**8))
    plt.plot([1,20],[Mmax,Mmax],'c:',lw=2)
    plt.plot([1,20],[Mmax,Mmax],'c:',lw=2)
    plt.legend(loc='best')
    plt.savefig('PLOTS/fMk_boxplots_'+filenames[0]+'_'+filenames[1]+'.png',bbox_inches='tight')


    #f1 = h5py.File(karlsson.get_mmixgrid_filename('minihalo'),'r')
    #tarr1 = np.array(f1['tarr']); f1.close()
    #Mhalo,zvir,vt,lt,nsn,trec = karlsson.params_minihalo()
    #th1 = TopHat(Mhalo=Mhalo,zvir=zvir)
    #rho1 = th1.get_rho_of_t_fn()
    
    #f2 = h5py.File(karlsson.get_mmixgrid_filename('atomiccoolinghalo'),'r')
    #tarr2 = np.array(f2['tarr']); f2.close()
    #Mhalo,zvir,vt,lt,nsn,trec = karlsson.params_atomiccoolinghalo()
    #th2 = TopHat(Mhalo=Mhalo,zvir=zvir)
    #rho2 = th2.get_rho_of_t_fn()
    
    #fig2 = plt.figure()
    #ax2 = fig2.gca()
    
    #ax2.plot(tarr1,rho1(tarr1),'k',lw=2)
    #ax2.plot(tarr2,rho2(tarr2),'r',lw=1)
    #ax2.set_xscale('log')
    #ax2.set_yscale('log')
    #plt.savefig('PLOTS/tophat.png',bbox_inches='tight')

    plt.show()
