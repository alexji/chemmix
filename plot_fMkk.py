import numpy as np
import pylab as plt
import karlsson
import util
from scipy import stats
from optparse import OptionParser

def plot_boxplot(karr,xdistr,ydistr,name,color,ax=None,dx=.1):
    if ax==None: ax = plt.gca()
    minarr = []; maxarr = []; meanarr = []
    leftarr = []; q1arr = []; medarr = []; q3arr = []; rightarr = []
    for k in karr:
        ik = k-1
        distr = stats.rv_discrete(name=name,values=(xdistr,ydistr[ik,:]))
        minarr.append(distr.ppf(0))
        leftarr.append(distr.ppf(.05))
        q1arr.append(distr.ppf(.25))
        meanarr.append(distr.mean())
        medarr.append(distr.median())
        q3arr.append(distr.ppf(.75))
        rightarr.append(distr.ppf(.95))
        maxarr.append(distr.ppf(1))
    minarr  = np.array(minarr)
    leftarr = np.array(leftarr)
    q1arr   = np.array(q1arr)
    meanarr = np.array(meanarr)
    medarr  = np.array(medarr)
    q3arr   = np.array(q3arr)
    rightarr= np.array(rightarr)
    maxarr  = np.array(maxarr)
    
    ax.plot(karr,leftarr,linestyle='--',color=color)
    ax.plot(karr+dx,meanarr,marker='s',linestyle='None',color=color)
    ax.errorbar(karr,medarr,yerr=[medarr-q1arr,q3arr-medarr],
                marker='o',ls='-',lw=1.8,label=name,color=color)
    ax.plot(karr,rightarr,linestyle='--',color=color)
    return ax

if __name__=="__main__":
    parser = OptionParser()
    options,args = parser.parse_args()

    envname,sfrname = args
    if 'atomiccoolinghalo' in envname: 
        Mhalo = 10**8
        fb = 0.1551; Mmax = Mhalo*fb
        xlim=((10**5,Mmax*10**0.1))
        ylim=(0,.07)
    if 'minihalo' in envname: 
        Mhalo = 10**6
        fb = 0.1551; Mmax = Mhalo*fb
        ylim=(0,.1)
        xlim=((10**4.0,10**5.5))
    
    k2max,k3max,cgrid = util.load_ckk(envname,sfrname)
    Mplot = np.load(util.fname_Mplotkk(envname,sfrname))
    fMkk2 = np.load(util.fname_fMkkII(envname,sfrname))
    fMkk3 = np.load(util.fname_fMkkIII(envname,sfrname))
    fig,axes = plt.subplots(k2max+1,k3max+1,sharex=True,sharey=True,figsize=(20,20))
    fig.subplots_adjust(hspace=0,wspace=0)
    for kII in range(k2max+1):
        for kIII in range(k3max+1):
            ax = axes[kII,kIII]
            if kII==0: ax.set_title(r'$k_{III}='+str(kIII)+r'$')
            if kIII==0: ax.set_ylabel(r'$k_{II}='+str(kII)+r'$')
            if kII+kIII==0: 
                continue
            ax.plot([Mmax,Mmax],[0,1.],'c:',lw=2)
            ax.plot(Mplot,fMkk2[kII,kIII,:])
            ax.plot(Mplot,fMkk3[kII,kIII,:])
            ax.set_xscale('log')
            ax.set_xlim(xlim)
            ax.set_ylim(ylim); ax.set_yticklabels([])
            if kII==k2max: ax.set_xlabel('Mmix')
            #print kII,kIII,np.nansum(fMkk2[kII,kIII,:]),np.nansum(fMkk3[kII,kIII,:])
    plt.savefig('PLOTS/fMkkgrid_'+envname+'_'+sfrname+'.png')

#    for kp in range(kpmax+1):
#        karr = np.arange(kp+1,kmax+1)
#        ax = plot_boxplot(karr,Mplot,fMkkp[:,kp,:],str(kp),'black')
#        ax.set_xlim((0,kmax+1)); ax.set_xlabel(r'$k$')
#        ax.set_yscale('log')
#        ax.set_ylim((10**4,10**8))
#        #ax.plot([1,kmax],[Mmax,Mmax],'c:',lw=2)
#        filename=util.default_filename(envname,sfrname)+'_kp'+str(kp)+'.png'
#        plt.savefig('PLOTS/'+filename)
#        plt.clf()
