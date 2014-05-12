import numpy as np
from math import factorial
from scipy.integrate import quad,trapz
from scipy.interpolate import interp1d

def calc_sigE(t,rho_of_t,Mdil):
    """ Calculate initial dilution radius for Vmix_karlsson08 (kpc) """
    return ((Mdil*3.*1.989*10**33)/(rho_of_t(t)*4*np.pi))**(1./3) / (3.086 * 10**21)

def get_sigE_fn(rho_of_t,Mdil=10**5.):
    """ Calculate initial dilution radius for Vmix_karlsson08 (kpc) """
    return lambda t: calc_sigE(t,rho_of_t,Mdil)

def Vmix_karlsson05(t,sigma_mix=7*10**-4):
    """ Calculate Vmix from Karlsson 2005
        sigma_mix = 7 * 10**-4 kpc^2/Myr """
    return 4*np.pi/3 * (sigma_mix * t)**1.5

def Vmix_karlsson08(t,sigEfn,Dt=3.3*10**-4):
    """ Calculate Vmix from Karlsson et al 2008
        Default Dt = 3.3 * 10**-4 kpc^2/Myr (turbulent velocity ~ 10 km/s)
        Use get_sigE_fn to get sigE function """
    return 4*np.pi/3 * (6 * Dt * t + sigEfn(t)**2)**1.5

def wII(k,t,mu_of_t):
    """ Normalized pop II weighting function wII (Poisson) """
    mu = mu_of_t(t)
    return np.exp(-1.*mu) * mu**k / factorial(k)

def wIII(k,t,mu_of_t):
    """ UNNORMALIZED pop III weighting function wIII (modified from Poisson) """
    mu = mu_of_t(t)
    return np.exp(-1.*mu) * mu**k / math.factorial(k) * np.exp(-1.*(mu-k)**2)

def fIII(k,kp,uII,tmin=0,tHubble=10000):
    """ UNNORMALIZED star counts, used to weight conditional
        metallicity distributions (Karlsson et al 2008) """
    return quad(lambda t: wIII(kp,t) * wII(k-kp,t) * uII(t),tmin,tHubble)[0]

def nk(k,uII,tmin=0,tHubble=10000):
    """ UNNORMALIZED star counts for simple random SN,
        used to weight conditional metallicity distributions (Karlsson 05) """
    return quad(lambda t: wII(k,t) * uII(t), tmin, tHubble)[0]

def mu_avgsn(t,Vmix,uSN,tmin=0):
    """ mu(t) = Integrate[Vmix(t-tp) * u_SN(tp),{tp,tmin,t}] """
    if t==tmin:
        return 0
    return quad(lambda tp: Vmix(t-tp)*uSN(tp),tmin,t)[0]

#def Vmix_distr(t,tau,sigEfn):
#    if t < tau:
#        return 0
#    else:
#        return Vmix_karlsson08(tau,sigEfn)

def Mmix_distr(Mbins,k,Vmix,wISM,uSN,rho_of_t,tmin=0,tHubble=10**4.,Nres_dVdt=200,Nres_Mgrid=1000,returnMgrid=False):
    ## TODO
    
    #numerically calculate function for dVmix/dt; calculate Mmix
    tarr = np.linspace(tmin,tHubble,Nres_dVdt)
    Vgrid = [Vmix(t) for t in tarr]
    dVdtgrid = (np.array(Vgrid[2:Nres_dVdt])-np.array(Vgrid[0:Nres_dVdt-2]))/(tarr[2]-tarr[0])
    dVdtgrid = np.concatenate(([(Vgrid[1]-Vgrid[0])/(tarr[1]-tarr[0])],dVdtgrid,[(Vgrid[-1]-Vgrid[-2])/(tarr[-1]-tarr[-2])]))
    dVdt = interp1d(tarr,dVdtgrid)
    def Mmix(t,tau): #units of Msun
        return ((3.086 * 10**21)**3./(1.989*10**33)) * quad(lambda tp: dVdt(max(tp-t+tau,tmin))*rho_of_t(tp),max(t-tau,tmin),t)[0]
    #calculate Mmix on t-tau grid
    tarr = np.linspace(tmin,tHubble,Nres_Mgrid); dt = tarr[1]-tarr[0]
    tauarr = np.linspace(tmin,tHubble,Nres_Mgrid); dtau = tauarr[1]-tauarr[0]
    Mgrid = np.zeros((len(tarr),len(tauarr)))
    for i,t in enumerate(tarr):
        for j,tau in enumerate(tauarr):
            if t > tau:
                Mgrid[i,j] = Mmix(t,tau)
    #return distribution, using Mbins
    nout = len(Mbins)-1
    outarr = np.zeros(nout)
    for i in xrange(nout):
        Mlo = Mbins[i]; Mhi = Mbins[i+1]
        iix,iiy = np.where(np.logical_and(Mlo < Mgrid, Mgrid < Mhi))
        for it,itau in zip(iix,iiy):
            t = tarr[it]; tau = tauarr[itau]
            outarr[i] += dt*dtau*Vmix(tau)*wISM(k-1,t)*uSN(t-tau)*uSN(t)
    return dVdtgrid,Mgrid,outarr

def Mmix_distr_constrho(wISM,Mbins=np.linspace(0,6,601)*(10**6),uSN=0.25,k=10,tHubble=10**4.):
    def tauv(M,rho=2.05*10**-25,sigma=7*10**-4):
        return (3*M*2*10**33/(4*np.pi*rho*(3.086*10**21)**3))**1.5 /sigma
    fMk = np.zeros(len(Mbins)-1)
    for i,M in enumerate(Mbins[:-1]):
        tau = tauv(M)
        #fMk[i] = M*quad(lambda t: wISM(k-1,t)*uSN*uSN,tau,tHubble)[0]
        fMk[i] = M*trapz([wISM(k-1,t) for t in np.arange(tau,tHubble)]
    return fMk

def conditional_k_star_list(k,bins,Nstars = 10**6):
    star_array = []
    for n in xrange(Nstars):
        ## TODO generate k SN, and calculate A/H for this star (for all A)
        ## then add the A/H list to star_array
        pass
    return star_array
