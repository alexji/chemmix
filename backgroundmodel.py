import numpy as np
import sys
import time

from multiprocessing import Pool
import functools

import util
from alexinterp import interp1d
from tophat import TopHat
from scipy.interpolate import InterpolatedUnivariateSpline
from scipy.integrate import trapz

def get_Vmixfn(tarr,th,Dt,Vmax,getVmixdot=False,E51=1.):
    sigEfn = get_sigEfn(tarr,th,E51=E51)
    Vmixarr = Vmix_base(tarr,Dt,sigEfn,Vmax)
    Vmixfn = interp1d(tarr,Vmixarr)
    if getVmixdot:
        Vmixfnspline = InterpolatedUnivariateSpline(tarr,Vmixarr)
        Vmixdotspline = Vmixfnspline.derivative()
        Vmixdotarr = Vmixdotspline(tarr)
        Vmixdotfn = interp1d(tarr,Vmixdotarr)
        return Vmixfn,Vmixdotfn
    else:
        return Vmixfn
def get_Vmixarr(tarr,th,Dt,Vmax,getVmixdot=False,E51=1.):
    sigEfn = get_sigEfn(tarr,th,E51=E51)
    Vmixarr = Vmix_base(tarr,Dt,sigEfn,Vmax)
    if getVmixdot:
        Vmixfnspline = InterpolatedUnivariateSpline(tarr,Vmixarr)
        Vmixdotspline = Vmixfnspline.derivative()
        Vmixdotarr = Vmixdotspline(tarr)
        #dt = tarr[1]-tarr[0]
        #Vmixdotarr = (Vmixarr[1:]-Vmixarr[:-1])/dt
        #Vmixdotarr = np.concatenate((Vmixdotarr,[0]))
        return Vmixarr,Vmixdotarr
    else:
        return Vmixarr

def Vmix_base(t,Dt,sigEfn,Vmax):
    return 4*np.pi/3. * (6*Dt*t + sigEfn(t)**2)**(1.5) #kpc^3
    #return np.minimum(4*np.pi/3. * (6*Dt*t + sigEfn(t)**2)**(1.5),Vmax)
def get_sigEfn(tarr,th,E51=1,zm=10**-3,betaC06=1): #Cioffi et al 1988 Rmerge in kpc
    #nfn = th.get_n_of_t_fn()
    #narr = nfn(tarr)
    #sigEarr = .069 * E51**.316 * (narr**-.367) * zm**-.051 * betaC06**-.429
    sigE = .069 * E51**.316 * (th.nvir**-.367) * zm**-.051 * betaC06**-.429
    sigEarr = np.zeros(tarr.shape)+sigE
    return interp1d(tarr,sigEarr)

def _compute_mmix_grid(tasknum,dt,Vmixdotarr,rhoarr,Mmax):
    """ trapezoid rule on fixed timegrid """
    # The Mmax is needed because rho -> infinity at t=0; the SFR should prevent this situation
    it,itau = divmod(tasknum,len(Vmixdotarr))
    if itau>it: return 0
    if itau==0: return 0
    return min(Mmax,dt*(.5*Vmixdotarr[0]*rhoarr[it-itau] + np.sum(Vmixdotarr[1:itau]*rhoarr[(it-itau+1):it]) + .5*Vmixdotarr[itau]*rhoarr[it]))

def compute_mmix_grid(envname,numprocs=1,hires=False,timethis=True):
    tarr = gettarr(envname,hires=hires)
    th,rhofn,Vmixfn,Vmixdotfn = util.load_Vmix(envname,get_thrho=True,getVmixdot=True)
    rhoarr = rhofn(tarr)*1.477543 * 10**31 #g/cm^3 to Msun/kpc^3
    Vmixdotarr = Vmixdotfn(tarr) #kpc^3
    numtasks = len(tarr)**2

    Mhalo,zvir,vturb,lturb,nSN,trecovery,Dt,uSN,E51 = envparams(envname)
    Mmax = Mhalo*0.1551 #fb = .1551

    dt = float(tarr[1]-tarr[0])
    myfunc = functools.partial(_compute_mmix_grid,
                               dt=dt,Vmixdotarr=Vmixdotarr,
                               rhoarr=rhoarr,Mmax=Mmax)
    pool = Pool(numprocs)
    if timethis: start = time.time()
    out = pool.map(myfunc,xrange(numtasks))
    out = np.reshape(np.array(out),(len(tarr),len(tarr)))
    if timethis: print "compute_mmix_grid time: %f" % (time.time()-start)
    pool.close()
    return out

def _compute_ckk(tasknum,k3max,wII,muII,wIII,muIII,uII,tarr):
    if tasknum==0: return 0
    kII,kIII = divmod(tasknum,k3max+1)
    f = lambda t: wIII(kIII,muIII(t))*wII(kII,muII(t))*uII(t)
    y = f(tarr); y[np.isnan(y)] = 0
    return trapz(y,x=tarr)

def compute_ckk(envname,sfrname,
                 k2max,k3max,
                 numprocs=1):

    uII,uIII,muII,muIII,wII,wIII = util.load_sfr(envname,sfrname)
    tarr = gettarr(envname) #np.linspace(tmin,tmax+.01,.01)

    pool = Pool(numprocs)
    myfunc = functools.partial(_compute_ckk,k3max=k3max,
                               wII=wII,muII=muII,wIII=wIII,muIII=muIII,
                               uII=uII,tarr=tarr)
    cgrid = pool.map(myfunc,range((k2max+1)*(k3max+1)))
    cgrid = np.array(cgrid).reshape((k2max+1,k3max+1))
    pool.close()
    return cgrid/np.nansum(cgrid)

def _get_DMlist(i,Mbins,tarr,Mmixgrid):
    Mmin = Mbins[i]; Mmax = Mbins[i+1]
    iix,iiy = np.where(np.logical_and(Mmixgrid > Mmin, Mmixgrid <= Mmax))
    these_t   = tarr[iix]
    these_tau = tarr[iiy]
    return np.array([these_t,these_tau]).transpose()

def get_DMlist(Mbins,tarr,Mmixgrid,numprocs=1):
    nbins = len(Mbins)-1
    myfunc = functools.partial(_get_DMlist,
                               Mbins=Mbins,tarr=tarr,Mmixgrid=Mmixgrid)
    pool = Pool(numprocs)
    DMlist = pool.map(myfunc,range(nbins))
    pool.close()
    return DMlist

def _compute_fMkkII(DM,k2max,k3max,Vmix,uII,uIII,muII,muIII,wII,wIII,dt):
    outarr = np.zeros((k2max+1,k3max+1,1))
    t = DM[:,0]; tau = DM[:,1]
    Vmixarr = Vmix(tau)
    muIIarr = muII(t)
    muIIIarr= muIII(t)
    lmssfr = uII(t)
    uIIarr = uII(t-tau)
    preprod = dt*dt*Vmixarr*lmssfr*uIIarr
    for kII in range(k2max+1):
        if kII==0:
            outarr[kII,:,0] = 0
        else:
            for kIII in range(k3max+1):
                wprod = wIII(kIII,muIIIarr)*wII(kII-1,muIIarr)
                outarr[kII,kIII,0] = np.sum(preprod*wprod)
    return outarr
def _compute_fMkkIII(DM,k2max,k3max,Vmix,uII,uIII,muII,muIII,wII,wIII,dt):
    start = time.time()
    outarr = np.zeros((k2max+1,k3max+1,1))
    t = DM[:,0]; tau = DM[:,1]
    Vmixarr = Vmix(tau)
    muIIarr = muII(t)
    muIIIarr= muIII(t)
    lmssfr = uII(t)
    uIIIarr = uIII(t-tau)
    preprod = dt*dt*Vmixarr*lmssfr*uIIIarr
    for kII in range(k2max+1):
        for kIII in range(k3max+1):
            if kIII==0: 
                outarr[kII,kIII,0] = 0
            else:
                wprod = wIII(kIII-1,muIIIarr)*wII(kII,muIIarr)
                outarr[kII,kIII,0] = np.sum(preprod*wprod)
    print "lenDM: %7i time: %6.1f" % (len(DM),time.time()-start)
    sys.stdout.flush()
    return outarr

def _post_fMkk(fMkklist,k2max,k3max):
    fMkk = np.concatenate(fMkklist,axis=2)
    for kII in range(k2max+1):
        for kIII in range(k3max+1):
            fMkk[kII,kIII,:] = fMkk[kII,kIII,:]/np.sum(fMkk[kII,kIII,:])
    return fMkk

def compute_fMkk(envname,sfrname,
                  k2max,k3max,
                  logdM = .01,
                  numprocs=1,hires=False,timethis=True):
    allstart = time.time()
    tarr = gettarr(envname,hires=hires)
    dt = float(tarr[1]-tarr[0])
    logMmin = 2
    if 'minihalo' in envname: logMmax = 6
    if 'atomiccoolinghalo' in envname: logMmax = 8
    Mbins,Mplot = util.logMbinsMplot(logMmin,logMmax,logdM)
    np.save(util.fname_Mbinskk(envname,sfrname),Mbins)
    np.save(util.fname_Mplotkk(envname,sfrname),Mplot)

    Vmixfn = util.load_Vmix(envname)
    uII,uIII,muII,muIII,wII,wIII = util.load_sfr(envname,sfrname)

    mmixgrid = util.load_mmixgrid(envname)
    print "Getting DMlist..."; start = time.time()
    DMlist = get_DMlist(Mbins,tarr,mmixgrid,numprocs=numprocs)
    print "Time: %f" % (time.time()-start)

    pool = Pool(numprocs)
    myfuncII = functools.partial(_compute_fMkkII,
                               k2max=k2max,k3max=k3max,Vmix=Vmixfn,
                               uII=uII,uIII=uIII,muII=muII,muIII=muIII,
                               wII=wII,wIII=wIII,dt=dt)
    myfuncIII = functools.partial(_compute_fMkkIII,
                               k2max=k2max,k3max=k3max,Vmix=Vmixfn,
                               uII=uII,uIII=uIII,muII=muII,muIII=muIII,
                               wII=wII,wIII=wIII,dt=dt)
    print "Starting fII"; start = time.time()
    fMkkII  = pool.map(myfuncII,DMlist)
    fMkkII  = _post_fMkk(fMkkII,k2max,k3max)
    np.save(util.fname_fMkkII(envname,sfrname),fMkkII)
    print "Finished fII %f" % (time.time()-start)

    print "Starting fIII"; start = time.time()
    fMkkIII = pool.map(myfuncIII,DMlist)
    fMkkIII  = _post_fMkk(fMkkIII,k2max,k3max)
    np.save(util.fname_fMkkIII(envname,sfrname),fMkkIII)
    print "Finished fIII %f" % (time.time()-start)

    pool.close()
    print "Total time: %f" % (time.time()-allstart)

def gettarr(envname,hires=False):
    if hires: dt = .01
    else: dt = 0.1
    if 'atomiccoolinghalo' in envname:
        tmin = 0
        tmax = 1000
    elif 'minihalo' in envname:
        tmin = 0
        tmax = 200
    return np.arange(tmin,tmax+dt,dt)

def envparams(envname,gettarr=False):
    if envname=='atomiccoolinghalo':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        E51 = 1.0
    #elif envname=='atomiccoolinghalo4':
    #    Mhalo = 10**8; zvir = 10
    #    nSN = 10; trecovery = 300
    elif envname=='atomiccoolinghalo_lowDt':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        E51 = 1.0
    #elif envname=='atomiccoolinghalo_lowDt4':
    #    Mhalo = 10**8; zvir = 10
    #    nSN = 10; trecovery = 300
    elif envname=='atomiccoolinghalo_lowmass':
        Mhalo = 10**7.4; zvir = 10
        nSN = 10; trecovery = 300
        E51 = 1.0
    #elif envname=='atomiccoolinghalo_lowmass4':
    #    Mhalo = 10**7.4; zvir = 10
    #    nSN = 10; trecovery = 300
    elif envname=='minihalo':
        Mhalo = 10**6; zvir = 25
        nSN = 2; trecovery = 10
        E51 = 1.0
    #elif envname=='minihalo4':
    #    Mhalo = 10**6; zvir = 25
    #    nSN = 2; trecovery = 10
    elif envname=='minihalo_lowDt':
        Mhalo = 10**6; zvir = 25
        nSN = 2; trecovery = 10
        E51 = 1.0
    elif envname=='minihalo_lowE':
        Mhalo = 10**6; zvir = 25
        nSN = 2; trecovery = 10
        E51 = 0.01 #10^49 ergs
    else:
        raise ValueError("invalid envname: "+envname)

    th = TopHat(Mhalo=Mhalo,zvir=zvir)
    lturb = th.Rvir/10.; vturb = th.vvir
    if 'lowDt' in envname: vturb = vturb/10.

    vturbkpcMyr = vturb*3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturbkpcMyr * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3

    #rhovir = th.nvir * th.mu * 24705449.8 #Msun/kpc^3
    #Vmax = (th.fb*Mhalo)/(rhovir)
    return Mhalo,zvir,vturb,lturb,nSN,trecovery,Dt,uSN,E51
