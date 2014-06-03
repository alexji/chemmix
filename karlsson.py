import numpy as np
from math import factorial
from scipy.integrate import quad,trapz
from scipy.interpolate import interp1d
import numpy.random as random

from tophat import TopHat
import yields

import time
from multiprocessing import Pool
import functools
import warnings

def MbinsMplot(Mmin,Mmax,dM):
    Mbins = np.arange(Mmin,Mmax+dM,dM)
    Mplot = np.arange(Mmin,Mmax,dM)+dM/2.
    return Mbins,Mplot
def logMbinsMplot(logMmin,logMmax,logdM):
    Mbins = 10.**np.arange(logMmin,logMmax+logdM/2.,logdM)
    Mplot = 10.**(np.arange(logMmin,logMmax,logdM)+logdM/2.)
    return Mbins,Mplot
    
def params_minihalo():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**6; zvir = 25; nSN = 2; trecovery = 10
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
    #return 10.**6,25,2.,.01,2,10
def params_atomiccoolinghalo():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 10; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
    #return 10.**8,10,10.,.1,10,300

def get_mmixgrid_filename(filename,mmixgrid_foldername = "MMIXGRIDS"):
    return mmixgrid_foldername+'/'+filename+'_mmixgrid.hdf5'
def get_fMk_filename(filename,k,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+filename+'_fM'+str(k)+'.npy'
def get_Mbins_filename(filename,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+filename+'_Mbins.npy'
def get_Mplot_filename(filename,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+filename+'_Mplot.npy'
def get_chemgrid_filename(filename,chemgrid_foldername = "CHEMGRIDS",postfix=''):
    return chemgrid_foldername+'/'+filename+'_chemgrid'+postfix+'.hdf5'

def calc_sigE(t,rho_of_t,Mdil):
    """ Calculate initial dilution radius for Vmix_karlsson08 (kpc) """
    return ((Mdil*3.*1.989*10**33)/(rho_of_t(t)*4*np.pi))**(1./3) / (3.086 * 10**21)

def get_sigEfn_K08(rho_of_t,Mdil=10**5):
    """ Calculate initial dilution radius for Vmix_karlsson08 (kpc) """
    return lambda t: calc_sigE(t,rho_of_t,Mdil)

def Vmix_K08(t,sigEfn,Dt=3.3*10**-4):
    """ Calculate Vmix from Karlsson et al 2008
        Default Dt = 3.3 * 10**-4 kpc^2/Myr (turbulent velocity ~ 10 km/s)
        Use get_sigE_fn to get sigE function """
    return 4*np.pi/3 * (6 * Dt * t + sigEfn(t)**2)**1.5

def get_Vmixfn_K08(rho_of_t,Mdil=10**5,Dt=3.3*10**-4,const_sigE=True):
    """ Calculate Vmix from Karlsson et al 2008
        Default Dt = 3.3 * 10**-4 kpc^2/Myr (turbulent velocity ~ 10 km/s) """
    if const_sigE: #corresponds to mu=4/3, n=0.1/cc
        sigEfn = lambda t: 0.1936*(Mdil/10**5.)**(1./3)
    else:
        sigEfn = get_sigEfn_K08(rho_of_t,Mdil)
    return lambda t: Vmix_K08(t,sigEfn,Dt=Dt)

def uSN_G08_const(t):
    """ 10 SN in ~400 Myr over ~5 kpc^3 = .005 SN/kpc^3/Myr """
    return .005

def wISM_K05(k,mu):
    return np.exp(-1.*mu)*mu**k/factorial(k)

def Vmix_K05(t,sigma = 7*10**-4):
    """
    t in Myr, sigma in kpc^2/Myr
    returns volume in kpc^3
    """
    return (4*np.pi/3) * (sigma * t)**1.5

def Vmixdot_K05(t,sigma = 7*10**-4):
    """
    t in Myr, sigma in kpc^2/Myr
    returns volume derivative in kpc^3/Myr
    """
    return 2*np.pi * sigma * (sigma * t)**0.5

def rho_K05A(t):
    """ g/cm^3 (corresponds to n = 0.1/cm^3) """
    try:
        return 2.05 * 10**-25 + np.zeros(len(t))
    except:
        return 2.05 * 10**-25

def uSN_K05A(t):
    """ Number/(kpc^3 Myr) """
    return 0.25

def get_mufn_K05A(tmin=0,tmax=2000,dt=1.0):
    tarr = np.arange(tmin,tmax,dt)
    muarr = [quad(lambda tp: Vmix_K05(t-tp)*uSN_K05A(tp),tmin,t)[0] for t in tarr]
    return interp1d(tarr,muarr)

def get_mufn(Vmix,uSN,tmin=0,tmax=2000,dt=1.0):
    tarr = np.arange(tmin,tmax,dt)
    muarr = [quad(lambda tp: Vmix(t-tp)*uSN(tp),tmin,t)[0] for t in tarr]
    return interp1d(tarr,muarr)

def Mmix(t,tau,Vmixdot,rho):
    """ Mixing mass in Msun """
    if t <= tau:
        return 0
    else:
        return ((3.086 * 10**21)**3./(1.989*10**33)) * quad(lambda tp: Vmixdot(tp-t+tau)*rho(tp),t-tau,t)[0]

def calc_Mmix_grid(Vmixdot,rho,tmin=0,tmax=1000,dt=1.0,numprocs=1):
    """ Returns (tarr,tauarr,Mmixgrid) """
    if (tmax-tmin)/dt > 10**6:
        print "Warning: length of tarr is %i (Mmixgrid is that squared)" % (int((tmax-tmin)/dt))
    tarr   = np.arange(tmin,tmax,dt) + dt/2.
    tauarr = np.arange(tmin,tmax,dt) + dt/2.
    Mmixgrid = np.zeros((len(tarr),len(tauarr)))
    if numprocs == 1:
        for i,t in enumerate(tarr):
            for j,tau in enumerate(tauarr):
                Mmixgrid[i,j] = Mmix(t,tau,Vmixdot,rho)
    else:
        assert numprocs > 1
        pool = Pool(numprocs)
        for j,tau in enumerate(tauarr):
            Mmixer = functools.partial(Mmix,tau=tau,Vmixdot=Vmixdot,rho=rho)
            Mmixgrid[:,j] = pool.map(Mmixer,tarr)
        pool.close(); pool.join()
    return tarr,tauarr,Mmixgrid

def calc_Mmix_grid_fast(Vmix,rho,tmin=0,tmax=1000,dt=1.0,Vmixdotarr=None):
    """ Returns (tarr,tauarr,Mmixgrid) """
    if (tmax-tmin)/dt > 10**6:
        print "Warning: length of tarr is %i (Mmixgrid is that squared)" % (int((tmax-tmin)/dt))
    tarr = np.arange(tmin,tmax,dt) + dt/2.
    iarr = np.arange(len(tarr))
    Mmixgrid = np.zeros((len(tarr),len(tarr)))

    #Vmixdotarr = Vmixdot(tarr)
    if Vmixdotarr == None:
        Vmixarr = Vmix(np.concatenate((tarr,[tarr[-1]+dt])))
        Vmixdotarr = (Vmixarr[1:]-Vmixarr[:-1])/dt
    rhoarr = rho(tarr)
    for itau in iarr:
        for it in iarr[itau:]: #[(itau+1):]
            Mmixgrid[it,itau] = np.sum(Vmixdotarr[0:itau]*rhoarr[(it-itau):it])
    Mmixgrid = Mmixgrid * ((3.086 * 10**21)**3./(1.989*10**33)) * dt
    return tarr,tarr,Mmixgrid

def calc_Mmix_grid_fast_trap(Vmix,rho,tmin=0,tmax=1000,dt=1.0,Vmixdotarr=None):
    """ Returns (tarr,tauarr,Mmixgrid) """
    if (tmax-tmin)/dt > 10**6:
        print "Warning: length of tarr is %i (Mmixgrid is that squared)" % (int((tmax-tmin)/dt))
    #tarr = np.arange(tmin,tmax,dt) + dt/2.
    tarr = np.arange(tmin,tmax,dt) + dt/2.
    inttarr = np.arange(tmin,tmax+dt,dt)
    iarr = np.arange(len(tarr))
    Mmixgrid = np.zeros((len(tarr),len(tarr)))

    if Vmixdotarr == None:
        Vmixarr = Vmix(np.concatenate((inttarr,[inttarr[-1]+dt])))
        Vmixdotarr = (Vmixarr[1:]-Vmixarr[:-1])/dt
    rhoarr = rho(inttarr)
    for itau in iarr:
        for it in iarr[itau:]: #[(itau+1):]
            Mmixgrid[it,itau] = np.sum((Vmixdotarr[0:itau]*rhoarr[(it-itau):it]+Vmixdotarr[1:(itau+1)]*rhoarr[(it-itau+1):(it+1)])/2.0)
    Mmixgrid = Mmixgrid * ((3.086 * 10**21)**3./(1.989*10**33)) * dt
    return tarr,tarr,Mmixgrid

def get_one_DMlist(i,Mbins,tarr,tauarr,Mmixgrid,dt,dtau):
    Mmin = Mbins[i]; Mmax = Mbins[i+1]
    iix,iiy = np.where(np.logical_and(Mmixgrid > Mmin, Mmixgrid <= Mmax))
    these_t   = tarr[iix]
    these_tau = tauarr[iiy]
    these_dt  = np.zeros(len(iix)) + dt
    these_dtau= np.zeros(len(iiy)) + dtau
    return np.array([these_t,these_tau,these_dt,these_dtau]).transpose()

def get_DMlist(Mbins,tarr,tauarr,Mmixgrid,numprocs=1):
    nbins = len(Mbins)-1
    dt = tarr[1]-tarr[0]; dtau = dt # for now, same everywhere
    if numprocs==1:
        DMlist = []
        for i in xrange(nbins):
            DMlist.append(get_one_DMlist(i,Mbins,tarr,tauarr,Mmixgrid,dt,dtau))
    else:
        assert numprocs > 1
        pool = Pool(numprocs)
        this_func = functools.partial(get_one_DMlist,Mbins=Mbins,tarr=tarr,tauarr=tauarr,Mmixgrid=Mmixgrid,dt=dt,dtau=dtau)
        DMlist = pool.map(this_func,range(nbins))
        pool.close(); pool.join()
    return DMlist

def calc_fMk(k,Mbins,DMlist,Vmix,wISM,mufn,SFR,normalize=True,
             SFRlms=None):
    assert k > 0
    nbins = len(Mbins)-1
    assert nbins == len(DMlist)

    if SFRlms == None:
        SFRlms = SFR

    fMk = np.zeros(nbins)
    for i,DM in enumerate(DMlist):
        assert DM.shape[1] == 4
        t = DM[:,0]; tau = DM[:,1]; dt = DM[:,2]; dtau = DM[:,3]
        Vmixarr = Vmix(tau)
        muarr = mufn(t)
        wISMarr = wISM(k-1,muarr)
        sfr1 = SFR(t-tau)
        sfr2 = SFRlms(t)

        fMk[i] = np.sum(dt*dtau*Vmixarr*wISMarr*sfr1*sfr2)
    if normalize:
        fMk = fMk/np.sum(fMk)
    return fMk

def weight_chemgrid(kmax,chemgrid_in,paramfn,
                    elemnames=None,verbose=False):
    Nstar,numyields,kmax_tmp = chemgrid_in.shape
    assert kmax<=kmax_tmp
    chemarr = chemgrid_in[:,:,0:kmax]
    
    Mhalo,zvir,vturb,lturb,nSN,trecovery = paramfn()
    vturb *= 3.16/3.08 * .001 #km/s to kpc/yr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    if verbose:
        print "Mhalo %.1e vturb %.2e lturb %f" % (Mhalo,vturb,lturb)
        print "Dt %.2e uSN %.2e" % (Dt,uSN)
    th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = get_Vmixfn_K08(RHO,Dt=Dt)
    WISM = wISM_K05
    PSI = lambda t: uSN
    MUFN = get_mufn(VMIX,PSI)

    with warnings.catch_warnings():
        warnings.simplefilter('ignore') #ck gives a convergence warning which isn't important
        ck = calc_ck(kmax,WISM,MUFN,PSI)
    for k in range(kmax):
        chemarr[:,:,k] *= ck[k]
    chemarr = np.sum(chemarr,2)
    if elemnames != None:
        chemarr = convert_to_solar(elemnames,chemarr,verbose)
        return chemarr,ck

def convert_to_solar(elemnames,chemarr,verbose=False):
    chemarr = np.log10(chemarr)
    asplund09ZHsun = yields.map_elemnames_to_asplund09(elemnames)
    for z in range(len(elemnames)):
        chemarr[:,z] -= asplund09ZHsun[z]
        if verbose: 
            print "%s: min %3.2f max %3.2f" % (elemnames[z],np.min(chemarr[:,z]),np.max(chemarr[:,z]))
    return chemarr

def draw_from_distr(N,x,pdf,seed=None,eps=10**-10):
    if seed!=None:
        random.seed(seed)
    unifarr = random.rand(N)
    cdf = np.cumsum(pdf)
    assert np.abs(cdf[-1] - 1.0) < eps, 'cdf[-1] != 1.0 (to tolerance of %e): %f' % (eps,cdf[-1])
    cdf[-1]=1
    return np.array(x[np.searchsorted(cdf,unifarr)])

def calc_ck(kmax,wISM,mufn,SFR,tmin=0,tmax=1000):
    """ normalized array of weights of how many SN enrich a star """
    myarr =  [quad(lambda t: wISM(k,mufn(t))*SFR(t),tmin,tmax)[0] for k in range(1,kmax+1)]
    return np.array(myarr)/np.sum(myarr)
