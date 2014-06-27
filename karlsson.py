import numpy as np
from math import factorial
from scipy.integrate import quad,trapz
from scipy.interpolate import interp1d
import scipy.stats
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
    
def get_Dt_uSN(vturb,lturb,nSN,trecovery):
    vturb *= 3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    return Dt,uSN
def get_environment_fns(paramfn,logMdil=5,verbose=False):
    Mhalo,zvir,vturb,lturb,nSN,trecovery = paramfn()
    Dt,uSN = get_Dt_uSN(vturb,lturb,nSN,trecovery)
    if verbose:
        print "Mhalo %.1e zvir %f vturb %.2e lturb %f" % (Mhalo,zvir,vturb,lturb)
        print "Dt %.2e uSN %.2e" % (Dt,uSN)
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    RHO = th.get_rho_of_t_fn()
    VMIX = get_Vmixfn_K08(RHO,Dt=Dt,Mdil=10**logMdil)
    return th,RHO,VMIX

def envparams(envname):
    if envname=='atomiccoolinghalo':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=5
    elif envname=='atomiccoolinghalo4':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=4
    elif envname=='atomiccoolinghalo_lowmass':
        Mhalo = 10**7.4; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=5
    elif envname=='atomiccoolinghalo_lowmass4':
        Mhalo = 10**7.4; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=4
    elif envname=='minihalo':
        Mhalo = 10**6; zvir = 25
        nSN = 2; trecovery = 10
        logMdil=5
    elif envname=='minihalo4':
        Mhalo = 10**6; zvir = 25
        nSN = 2; trecovery = 10
        logMdil=4
    else:
        raise ValueError("invalid envname: "+envname)

    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil

def params_k08():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc
    Mhalo = 10**8; zvir = 10
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo,zvir,vturb,lturb
def params_minihalo():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**6; zvir = 25; nSN = 2; trecovery = 10
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_fatminihalo(): # TODO
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**7; zvir = 18; nSN = 2; trecovery = 10
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghaloearly():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 15; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghalo():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 10; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghalolate():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 7; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghalolate_lowvturb():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 7; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir/10.
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghalo_lowvturb():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**8; zvir = 10; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir/10.
    return Mhalo, zvir, vturb, lturb, nSN, trecovery
def params_atomiccoolinghalo_lowmass():
    #Mhalo/Msun, zvir, vturb/km/s, lturb/kpc, nSN, trecovery/Myr
    Mhalo = 10**7.4; zvir = 10; nSN = 10; trecovery = 300
    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    return Mhalo, zvir, vturb, lturb, nSN, trecovery

def get_mmixgrid_filename(envname,mmixgrid_foldername = "MMIXGRIDS"):
    return mmixgrid_foldername+'/'+envname+'_mmixgrid.hdf5'
def get_fMkkp_filename(envname,suffix,k,kp,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+envname+'_'+suffix+'_fM'+str(k)+'_'+str(kp)+'.npy'
def get_ckkp_filename(envname,suffix,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+envname+'_'+suffix+'_ckkp.hdf5'
def get_Mbinskkp_filename(envname,suffix,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+envname+'_'+suffix+'_Mbins.npy'
def get_Mplotkkp_filename(envname,suffix,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+envname+'_'+suffix+'_Mplot.npy'

def RHOP2f(filename,rhop2):
    if rhop2: return filename+'_rhop2'
    return filename
def get_fMk_filename(filename,k,rhop2=False,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+RHOP2f(filename,rhop2)+'_fM'+str(k)+'.npy'
def get_Mbins_filename(filename,rhop2=False,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+RHOP2f(filename,rhop2)+'_Mbins.npy'
def get_Mplot_filename(filename,rhop2=False,mmixdistr_foldername = "MMIXDISTR"):
    return mmixdistr_foldername+'/'+RHOP2f(filename,rhop2)+'_Mplot.npy'
def get_chemgrid_filename(filename,chemgrid_foldername = "CHEMGRIDS",
                          rhop2=False,postfix=''):
    return chemgrid_foldername+'/'+RHOP2f(filename,rhop2)+'_chemgrid'+postfix+'.hdf5'

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
    #return scipy.stats.poisson.pmf(k,mu) 

def wISM_III(k,mu):
    """ NOT normalized!"""
    return np.exp(-(mu-k)**2/4.) * np.exp(-mu) * mu**k/factorial(k) 
    #return np.exp(-(mu-k)**2/4.) * scipy.stats.poisson.pmf(k,mu)

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

def get_mufn(Vmix,uSN,tarr=None,tmin=0,tmax=2000,dt=1.0):
    if tarr == None:
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

def calc_fMk(k,Mbins,DMlist,Vmix,WISM,mufn,SFR,normalize=True,
             SFRlms=None, WISM2_0=None):
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
        WISMarr = WISM(k-1,muarr)
        if WISM2_0 != None: #ugh note WISM uses mu, this uses t
            WISMarr = WISMarr * WISM2_0(t)
        sfr1 = SFR(t-tau)
        sfr2 = SFRlms(t)

        fMk[i] = np.sum(dt*dtau*Vmixarr*WISMarr*sfr1*sfr2)
    if normalize:
        fMk = fMk/np.sum(fMk)
    return fMk

def calc_fMkkp(k,kp,Mbins,DMlist,VMIX,
               WISMII,MUII,UII,
               WISMIII,MUIII,UIII,
               normalize=True):
    assert k > 0 and kp > 0 and kp <= k
    nbins = len(Mbins)-1
    assert nbins == len(DMlist)

    fMkkp = np.zeros(nbins)
    for i,DM in enumerate(DMlist):
        assert DM.shape[1] == 4
        t = DM[:,0]; tau = DM[:,1]; dt = DM[:,2]; dtau = DM[:,3]
        Vmixarr = VMIX(tau) #assuming same Vmix for II and III
        muIIarr = MUII(t)
        muIIIarr= MUIII(t)

        if k-kp-1 >= 0:
            wprodII = WISMIII(kp,muIIIarr)*WISMII(k-kp-1,muIIarr)*UII(t-tau)
        else: wprodII = 0
        if kp-1 >= 0:
            wprodIII= WISMIII(kp-1,muIIIarr)*WISMII(k-kp,muIIarr)*UIII(t-tau)
        else: wprodIII = 0

        sfr = UII(t)

        fMkkp[i] = np.sum(dt*dtau*Vmixarr*(wprodII+wprodIII)*sfr)
    if normalize:
        fMkkp = fMkkp/np.sum(fMkkp)
    return fMkkp

def hist_chemgrid(chemarr,bins,elemindex=None,verbose=False):
    if elemindex==None:
        assert len(chemarr.shape)==2
        kmax = chemarr.shape[1]
    else:
        assert len(chemarr.shape)==3,'[SN,elem,k]'
        assert (elemindex >= 0) and (elemindex < chemarr.shape[2])
        kmax = chemarr.shape[2]
    minbin = bins[0]; maxbin = bins[-1]

    histlist = []
    for k in range(kmax):
        if elemindex==None:
            if verbose: warn_bins(chemarr[:,k],minbin,maxbin,prefix='k='+str(k))
            h,x = np.histogram(chemarr[:,k],bins=bins)
        else:
            if verbose: warn_bins(chemarr[:,elemindex,k],minbin,maxbin,prefix='k='+str(k))
            h,x = np.histogram(chemarr[:,elemindex,k],bins=bins)
        histlist.append(h)
    return x,histlist

def cumweight_list(mylist,ck,eps=10**-8):
    kmax = len(mylist)
    assert len(ck)==kmax
    myarr = np.array(mylist)
    #with warnings.catch_warnings():
    #    warnings.simplefilter('ignore') #ck gives a convergence warning which isn't important
    #ck = calc_ck(kmax,WISM,MUFN,PSI,tmin=tmin,tmax=tmax,wISM2_0=wISM2_0)
    weightedlist = []
    for thiskmax in range(kmax):
        thisck = ck[0:thiskmax]; thisck = thisck/np.sum(thisck)
        assert np.sum(thisck)-1 < eps
        thislist = np.transpose(myarr[0:thiskmax,:])
        weightedlist.append(np.dot(thislist,thisck))
    return weightedlist

def convert_to_solar(elemnames,chemarr,verbose=False):
    assert len(elemnames)==chemarr.shape[1]
    chemarr = np.log10(chemarr)
    asplund09ZHsun = yields.map_elemnames_to_asplund09(elemnames)
    if len(chemarr.shape)==2: #single chemarr
        for z in range(len(elemnames)):
            chemarr[:,z] -= asplund09ZHsun[z]
            if verbose: 
                print "%s: min %3.2f max %3.2f" % (elemnames[z],np.min(chemarr[:,z]),np.max(chemarr[:,z]))
    elif len(chemarr.shape)==3: #chemarr with multiple k
        for k in range(chemarr.shape[2]):
            for z in range(len(elemnames)):
                chemarr[:,z,k] -= asplund09ZHsun[z]
                if verbose: 
                    print "%i %s: min %3.2f max %3.2f" % (k,elemnames[z],np.min(chemarr[:,z]),np.max(chemarr[:,z]))
    return chemarr

def draw_from_distr(N,x,pdf,seed=None,eps=10**-10):
    if seed!=None:
        random.seed(seed)
    unifarr = random.rand(N)
    cdf = np.cumsum(pdf)
    assert np.abs(cdf[-1] - 1.0) < eps, 'cdf[-1] != 1.0 (to tolerance of %e): %f' % (eps,cdf[-1])
    cdf[-1]=1
    return np.array(x[np.searchsorted(cdf,unifarr)])

def calc_ck(kmax,wISM,mufn,SFR,tmin=0,tmax=1000,wISM2_0=None):
    """ normalized array of weights of how many SN enrich a star """
    if wISM2_0==None:
        myarr =  [quad(lambda t: wISM(k,mufn(t))*SFR(t),tmin,tmax)[0] for k in range(1,kmax+1)]
    else:
        myarr =  [quad(lambda t: wISM(k,mufn(t))*wISM2_0(t)*SFR(t),tmin,tmax)[0] for k in range(1,kmax+1)]
    return np.array(myarr)/np.sum(myarr)

def _calc_ckkp(tasknum,kmax,wISMII,mufnII,wISMIII,mufnIII,uII,tmin,tmax):
    k,kp = divmod(tasknum,kmax+1)
    k+=1 #k >= 1; k >= kp >= 0
    if kp>k: return 0.
    return quad(lambda t: wISMIII(kp,mufnIII(t))*wISMII(k-kp,mufnII(t))*uII(t),tmin,tmax)[0]

def calc_ckkp(kmax,wISMII,mufnII,wISMIII,mufnIII,uII,tmin=0,tmax=1000,numprocs=1):
    """ normalized array of weights of how many SN enrich a star """
    if numprocs==1:
        cgrid = np.zeros((kmax,kmax+1))
        for ik in range(kmax):
            k = ik+1
            start = time.time()
            for ikp in range(k+1):
                kp = ikp
                cgrid[ik,ikp] = quad(lambda t: wISMIII(kp,mufnIII(t))*wISMII(k-kp,mufnII(t))*uII(t),tmin,tmax)[0]
            print "k=%i: %4.2f" % (k,time.time()-start)
    else:
        assert numprocs > 1
        pool = Pool(numprocs)
        this_func = functools.partial(_calc_ckkp,kmax=kmax,
                                      wISMII=wISMII,mufnII=mufnII,
                                      wISMIII=wISMIII,mufnIII=mufnIII,uII=uII,
                                      tmin=tmin,tmax=tmax)
        cgrid = pool.map(this_func,range(kmax*(kmax+1)))
        pool.close(); pool.join()
        cgrid = np.array(cgrid).reshape((kmax,kmax+1))
    return cgrid/np.nansum(cgrid)

def compute_cfrac(C,Fe,CFe_crit=0.75,bins=None,
                  returnallinfo=False):
    Nstar=len(C)
    assert Nstar == len(Fe)
    CFe = C - Fe
    
    if bins==None:
        binwidth = 0.2
        binmin = np.ceil(np.min(Fe))
        while binmin > np.min(Fe):
            binmin -= binwidth
        binmax = np.floor(np.max(Fe))
        while binmin > np.min(Fe):
            binmax += binwidth
        bins = np.arange(binmin,binmax+binwidth/2.,binwidth)
    binplot = (bins[0:-1] + bins[1:])/2.0
    assert len(binplot) == len(bins)-1
    fracarr = np.zeros(len(binplot))
    NCarr = np.zeros(len(binplot))
    NTarr = np.zeros(len(binplot))

    Crich = CFe > CFe_crit
    ii = np.digitize(Fe,bins)
    for ibin in xrange(1,len(binplot)+1):
        starsinbin = (ii == ibin)
        NCarr[ibin-1] = np.sum(Crich[starsinbin])
        NTarr[ibin-1] = np.sum(starsinbin)
    fracarr = NCarr/NTarr
    errarr  = [np.abs(fracarr - err/NTarr) for err in scipy.stats.poisson.interval(0.95,NCarr)]

    if returnallinfo:
        return binplot,fracarr,errarr,NCarr,NTarr
    else:
        return binplot,fracarr,errarr

def warn_bins(arr,minbin,maxbin,prefix=''):
    N = float(len(arr))
    minarr = np.min(arr)
    maxarr = np.max(arr)
    if minarr < minbin:
        percentless = np.sum(arr < minbin)/N * 100.
        print "Warning: "+prefix+" %3.2f below lowest bin %3.2f" % (k, percentless,minbin)
    if maxarr > maxbin:
        percentmore = np.sum(arr > maxbin)/N * 100.
        print "Warning: "+prefix+" %3.2f above highest bin %3.2f" % (k, percentmore,maxbin)

def cfraclist_chemgrid(chemarr,FeHbins,
                       iC=0,iFe=5,CFe_crit=0.75,verbose=False):
    kmax = chemarr.shape[2]
    minbin = FeHbins[0]; maxbin = FeHbins[-1]

    cfraclist = []
    for k in range(kmax):
        if verbose: warn_bins(chemarr[:,iFe,k],minbin,maxbin,prefix='k='+str(k))
        FeH,Cfrac,CfracErr = compute_cfrac(chemarr[:,iC,k],chemarr[:,iFe,k],
                                           CFe_crit=CFe_crit,bins=FeHbins)
        cfraclist.append(Cfrac)
    return FeH,cfraclist
