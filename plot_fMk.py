import numpy as np
import pylab as plt
import h5py

from scipy import stats
import karlsson
from tophat import TopHat


if __name__=="__main__":
    filenames = ['minihalo','atomiccoolinghalo']
    lefttaillist = []; righttaillist = []
    q1list = []; medlist = []; q3list = []
    karr = np.arange(1,21)
    offset = [0,.1]
    colors = ['k','r']

    fig = plt.figure(); ax = fig.gca()
    for i,filename in enumerate(filenames):
        Mplot = np.load(karlsson.get_Mplot_filename(filename))
        minarr = []; maxarr = []
        meanarr = []
        lefttailarr = []; righttailarr = []
        q1arr = []; medarr = []; q3arr = []
        for k in karr:
            fMk = np.load(karlsson.get_fMk_filename(filename,k))
            distr = stats.rv_discrete(name=str(k),values=(Mplot,fMk))
            minarr.append(distr.ppf(0))
            lefttailarr.append(distr.ppf(.05))
            q1arr.append(distr.ppf(.25))
            meanarr.append(distr.mean())
            medarr.append(distr.median())
            q3arr.append(distr.ppf(.75))
            righttailarr.append(distr.ppf(.95))
            maxarr.append(distr.ppf(1))
        #ax.plot(karr+offset[i],minarr,colors[i]+':')
        ax.plot(karr+offset[i],lefttailarr,colors[i]+'--')
        ax.plot(karr+2*offset[i],meanarr,colors[i]+'s',label=filename+' mean')
        ax.errorbar(karr+offset[i],medarr,yerr=[q1arr,q3arr],fmt=colors[i]+'o-',lw=1.8,label=filename+' med')
        ax.plot(karr+offset[i],righttailarr,colors[i]+'--')
        #ax.plot(karr+offset[i],maxarr,colors[i]+':')
        lefttaillist.append(lefttailarr); righttaillist.append(righttailarr)
        q1list.append(q1arr); medlist.append(medarr); q3list.append(q3arr)
    ax.set_xlim((0,21)); ax.set_xlabel(r'$k$')
    ax.set_yscale('log'); ax.set_ylabel(r'M$_{\rm mix}$ (q1,med,q3) $M_\odot$')
    #ax.set_ylim((10**-2.5,10**8.5))
    plt.legend(loc='best')



    f1 = h5py.File(karlsson.get_mmixgrid_filename('minihalo'),'r')
    tarr1 = np.array(f1['tarr']); f1.close()
    th1 = TopHat(Mhalo=10**6,nvir=0.1,fb=0.1551,mu=1.4)
    rho1 = th1.get_rho_of_t_fn()

    f2 = h5py.File(karlsson.get_mmixgrid_filename('atomiccoolinghalo'),'r')
    tarr2 = np.array(f2['tarr']); f2.close()
    th2 = TopHat(Mhalo=10**8,nvir=0.1,fb=0.1551,mu=1.4)
    rho2 = th2.get_rho_of_t_fn()

    fig2 = plt.figure()
    ax2 = fig2.gca()
    ax2.plot(tarr1,rho1(tarr1),'k',lw=2)
    ax2.plot(tarr2,rho2(tarr2),'r',lw=1)
    ax2.set_xscale('log')
    ax2.set_yscale('log')

    plt.show()
