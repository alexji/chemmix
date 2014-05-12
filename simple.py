import numpy as np
import tophat
from karlsson08 import * #Vmix_karlsson08,get_sigE_fn,nk,mu_avgsn
from scipy.integrate import quad
from scipy.interpolate import interp1d
import time

def simple_uII(t):
    """ Using Greif et al 2008, assume 10 Pop III SN over 300 Myr in (5kpc)^3"""
    return 10./(300. * 5**3)

def run_karlsson08():
    TH = tophat.TopHat(Mhalo=10**8,nvir=0.1,verbose=False)
    tophat_rho_of_t = TH.get_rho_of_t_fn(Nres=200)
    tmin = 2*np.min(tophat_rho_of_t.x)

    import pylab as plt
    tarr = np.linspace(10,700,300)
    #plt.plot(tarr,tophat_rho_of_t(tarr))
    #plt.gca().set_yscale('log')
    #plt.xlabel(r'$t$'); plt.ylabel(r'$\rho(t)$ (g/cm$^3$)')
    #plt.show()

    tarr = np.linspace(tmin,10000,500)
    sigEfn = get_sigE_fn(tophat_rho_of_t,10**5)
    Vmix = lambda t: Vmix_karlsson08(t,sigEfn)
    start = time.time()
    muarr = [mu_avgsn(t,Vmix,simple_uII,tmin=tmin) for t in tarr]
    mufn = interp1d(tarr,muarr)
    print "Computed mufn,",(time.time()-start)
    def wISM(k,t):
        return wII(k,t,mufn)
    #plt.plot(tarr,muarr)
    #plt.xlabel(r'$t$'); plt.ylabel(r'$\mu(t)$')
    #plt.show()
    
    print "calculating fMdistr"
    k=10
    Mbins = [0.1,1,10,100,1000,10000,10**6,10**7,10**8,10**9]
    start = time.time()
    dVdtgrid,Mgrid,fMdistr = Mmix_distr(Mbins,k,Vmix,wISM,simple_uII,tophat_rho_of_t,Nres_Mgrid=100,tmin=tmin)
    print "took:",(time.time()-start)

if __name__ == "__main__":
    rho_of_t = lambda t: 2.05 * 10**-25 #g cm^-3
    uconst = lambda t: 0.25 #kpc^-3 Myr^-1
    Vmix = Vmix_karlsson05
    tarr = np.linspace(0,10000,500)
    muarr = [mu_avgsn(t,Vmix,uconst) for t in tarr]
    mufn = interp1d(tarr,muarr)
    def wISM(k,t):
        return wII(k,t,mufn)
    
    k=10
    Mbins = np.linspace(0,6,301)*(10**6)
    #dVdtgrid,Mgrid,fMdistr = Mmix_distr(Mbins,k,Vmix,wISM,uconst,rho_of_t,Nres_Mgrid=100)
    fMk1 = Mmix_distr_constrho(wISM,k=11,Mbins=Mbins)
    fMk10 = Mmix_distr_constrho(wISM,k=10,Mbins=Mbins)
    fMk100 = Mmix_distr_constrho(wISM,k=100,Mbins=Mbins)

    import pylab as plt
    plt.plot(Mbins[:-1],fMk1,label='k=1')
    plt.plot(Mbins[:-1],fMk10,label='k=10')
    plt.plot(Mbins[:-1],fMk100,label='k=100')
    plt.legend(loc='best')
    plt.show()
    #plt.plot(Mbins[1:],fMdistr)
    #plt.show()

    ## TODO: figure out how to calculate the mixing mass distribution
    ## Then need a function to draw from that distribution
    ## Don't worry about the density profile right now

    #K = 10 #max number of SN enriching a star
    #star_arrays = []
    #bins = np.logspace(-6,1,71)
    #for k in range(1,K+1):
    #    star_lists.append(conditional_k_star_list(k,bins,Nstars=100))
