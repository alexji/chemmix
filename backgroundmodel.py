import numpy as np
import util
from alexinterp import interp1d
from tophat import TopHat
from scipy.interpolate import InterpolatedUnivariateSpline

def get_Vmixfn(tarr,rhofn,Dt,logMdil,Vmax,getVmixdot=False):
    sigEfn = get_sigEfn(tarr,rhofn,logMdil)
    Vmixarr = Vmix_base(tarr,Dt,sigEfn,Vmax)
    Vmixfn = interp1d(tarr,Vmixarr)
    if getVmixdot:
        #Vmixfnspline = InterpolatedUnivariateSpline(tarr,Vmixarr)
        #Vmixdotspline = Vmixfnspline.derivative()
        #Vmixdotarr = Vmixdotspline(tarr)
        dt = tarr[1]-tarr[0]
        Vmixdotarr = (Vmixarr[1:]-Vmixarr[:-1])/dt
        Vmixdotarr = np.concatenate((Vmixdotarr,[0]))
        Vmixdotfn = interp1d(tarr,Vmixdotarr)
        return Vmixfn,Vmixdotfn
    else:
        return Vmixfn
def get_Vmixarr(tarr,rhofn,Dt,logMdil,Vmax,getVmixdot=False):
    sigEfn = get_sigEfn(tarr,rhofn,logMdil)
    Vmixarr = Vmix_base(tarr,Dt,sigEfn,Vmax)
    if getVmixdot:
        #Vmixfnspline = InterpolatedUnivariateSpline(tarr,Vmixarr)
        #Vmixdotspline = Vmixfnspline.derivative()
        #Vmixdotarr = Vmixdotspline(tarr)
        dt = tarr[1]-tarr[0]
        Vmixdotarr = (Vmixarr[1:]-Vmixarr[:-1])/dt
        Vmixdotarr = np.concatenate((Vmixdotarr,[0]))
        return Vmixarr,Vmixdotarr
    else:
        return Vmixarr

def Vmix_base(t,Dt,sigEfn,Vmax):
    return np.minimum(4*np.pi/3. * (6*Dt*t + sigEfn(t)**2)**(1.5),Vmax)
def get_sigEfn(tarr,rhofn,logMdil):
    Mdil = 10**logMdil
    rhoarr = rhofn(tarr)
    sigEarr = ((3*Mdil)/(4*np.pi*rhoarr))**(1./3.) * 4.07 * 10**-11 #kpc
    return interp1d(tarr,sigEarr)

def compute_mmix_grid(envname):
    tarr = gettarr(envname)
    th,rhofn,vmixfn = util.load_Vmix(envname,get_thrho=True)

    pass


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
        logMdil=5
    elif envname=='atomiccoolinghalo4':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=4
    elif envname=='atomiccoolinghalo_lowDt':
        Mhalo = 10**8; zvir = 10
        nSN = 10; trecovery = 300
        logMdil=5
    elif envname=='atomiccoolinghalo_lowDt4':
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
    if 'lowDt' in envname: vturb = vturb/10.

    vturbkpcMyr = vturb*3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturbkpcMyr * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3

    rhovir = th.nvir * th.mu * 24705449.8 #Msun/kpc^3
    Vmax = (th.fb*Mhalo)/(rhovir)
    return Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,Dt,uSN,Vmax
