import numpy as np
import pylab as plt
import karlsson
import util
from scipy import stats

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
    envname = 'atomiccoolinghalo'
    sfrname = 'fidTS400'
    Mhalo = 10**8; fb = 0.1551; Mmax = Mhalo*fb
    
    kmax,ckkp = util.load_ckkp(envname,sfrname)
    Mplot = np.load(karlsson.get_Mplotkkp_filename(envname,sfrname))
    fMkkp = np.load(karlsson.get_fMkkp_filename(envname,sfrname))
    plt.figure()
    for kp in range(kmax+1):
        karr = np.arange(kp+1,kmax+1)
        ax = plot_boxplot(karr,Mplot,fMkkp[:,kp,:],str(kp),'black')
        ax.set_xlim((0,kmax+1)); ax.set_xlabel(r'$k$')
        ax.set_yscale('log')
        ax.set_ylim((10**4,10**8))
        ax.plot([1,kmax],[Mmax,Mmax],'c:',lw=2)
        filename=util.default_filename(envname,sfrname)+'_kp'+str(kp)+'.png'
        plt.savefig('PLOTS/'+filename)
        plt.clf()
