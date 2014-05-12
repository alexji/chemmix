import numpy as np
import pylab as plt
from karlsson import *

if __name__=="__main__":
    tmin = 0.0; tmax = 1000.0; dt = 0.1
    Nstars = 10**6
    print "tmax,dt:",tmax,dt
    print "Nstars:",Nstars
    fileprefix="Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    kmax=10
    chemarr = np.load(fileprefix+'_chemgrid_N'+str(Nstars)+'.npy')
    chemarr = chemarr[:,:,0:kmax]
    mufn = get_mufn_K05A()
    psi_K05 = uSN_K05A

    bins = np.arange(-6.5,0.05,.05)
    plt.subplots_adjust(hspace=0.25)
    plt.figure(figsize=(8,10))
    plt.title("[Fe/H]")

    for this_k in range(kmax):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore') #ck gives a convergence warning which isn't important
            ck = calc_ck(this_k,wISM_K05,mufn,psi_K05)
        this_chemarr = np.zeros(chemarr.shape)
        for k in range(this_k):
            this_chemarr[:,:,k] = chemarr[:,:,k] * ck[k]
        this_chemarr = np.log10(np.sum(this_chemarr,2))
        for k in range(1,this_k+1):
            plt.subplot(5,2,k)
            plt.hist(this_chemarr[:,(k-1)],bins=bins,normed=True)
            plt.xticks([-6,-5,-4,-3,-2,-1,0])
            plt.yticks([0,1,2,3,4]); plt.ylim((0,4))
    plt.show()
