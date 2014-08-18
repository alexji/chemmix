import numpy as np
from scipy import integrate
from alexinterp import interp1d
from tophat import TopHat
from optparse import OptionParser
import backgroundmodel as bgm
import util
import time
import functools
import sys
import scipy.special
import cosmology
from multiprocessing import Pool
#from numpy.fft import rfft,irfft
import sys

def compute_umufns(envname,sfrname,tarr,numprocs=1):
    th,Vmixfn = util.load_Vmix(envname,get_th=True)
    if 'Vmax' in sfrname:
        Mhalo = bgm.envparams(envname)[0]
        rhovir = th.nvir * th.mu * 24705449.8 #Msun/kpc^3
        rhomin = rhovir/8.
        Vmax = (th.fb*Mhalo)/(rhovir) #4*np.pi/3. * th.Rvir**3
        print "--compute_umufns: Using Vmax=%f" % (Vmax)
        Vmixfn = util.load_Vmix(envname,Vmax=Vmax)
    label = get_sfrlabel(envname,sfrname)
    u0II,u0III,ttII,ttIII,tstop = sfrparams(envname,sfrname)

    if 'burst' in sfrname:
        saveuIIIburst(label,tarr,Vmixfn,u0III,ttIII,50,numprocs=numprocs) #tburst=50Myr
        saveuIIdata(label,Vmixfn,th,u0II,ttII,tstop,numprocs=numprocs)
    elif 'h14' in sfrname:
        nSN = bgm.envparams(envname)[4]
        Rvir = th.Rvir
        saveuIIIh14(label,tarr,Vmixfn,nSN,Rvir,numprocs=numprocs)
        saveuIIdata(label,Vmixfn,th,u0II,ttII,tstop,numprocs=numprocs)
    elif 'betalogit' in sfrname:
        tvir = ttIII; trecovery = ttII-ttIII; tburst = trecovery
        saveuIIIbetalogit(label,tarr,Vmixfn,tvir,tburst,u0III,numprocs=numprocs)
        saveuIIbetalogit(label,tarr,Vmixfn,tvir,tburst,u0II,tstop,numprocs=numprocs)
    else:
        numsf = np.sum(tarr>ttIII); numskip = len(tarr)-numsf
        u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
        saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0III,ttIII,tstop,numprocs=numprocs)
        saveuIIdata(label,Vmixfn,th,u0II,ttII,tstop,numprocs=numprocs)

def plot_sfu(label,showplot=True):
    import pylab as plt
    fIII,tarr,uIII  = loaduIIIfn(label,retarrays=True)
    fII,tarr,uII    = loaduIIfn(label,retarrays=True)
    fIII,tarr,muIII = loadmuIIIfn(label,retarrays=True)
    fII,tarr,muII   = loadmuIIfn(label,retarrays=True)
    f,((ax1,ax2),(ax3,ax4)) = plt.subplots(2,2,sharex=True)
    ax1.plot(tarr,uIII); ax1.set_title(label+' uIII')
    ax2.plot(tarr,uII);  ax2.set_title(label+' uII')
    ax3.plot(tarr,muIII); ax3.set_title('muIII')
    ax4.plot(tarr,muII);  ax4.set_title('muII')
    plt.savefig('PLOTS/umugrid_'+label+'.png')
    if showplot: plt.show()

def saveuIIbetalogit(label,tarr,Vmixfn,tvir,tburst,u0II,tstop,numprocs=1):
    logitfn = get_logitfn(tvir,tburst)
    uarr = u0II*logitfn(tarr)
    uarr[tarr > tstop] = 0.
    ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def saveuIIdata(label,Vmixfn,th,u0,tthresh,tstop,numprocs=1):
    uIIIfn,tarr,uarr = loaduIIIfn(label,retarrays=True)
    muIIIfn = loadmuIIIfn(label)
    Qparr = Qp(muIIIfn(tarr))
    narr = (th.get_n_of_t_fn())(tarr)
    uarr = uII(Qparr,narr,u0=u0,nvir=th.nvir)
    uarr[tarr < tthresh] = 0.
    uarr[tarr > tstop]   = 0.

    ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = interp1d(tarr,muarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,muarr
    return fn

def saveuIIIburst(label,tarr,Vmixfn,u0III,ttIII,tburst,numprocs=1):
    uarr = np.zeros(len(tarr))
    ii = (tarr >= ttIII) & (tarr <= ttIII+tburst)
    uarr[ii] = u0III
    ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def saveuIIIh14(label,tarr,Vmixfn,nSN,Rvir,numprocs=1):
    h14sfr = get_h14sfrfn()
    h14norm = 2225.63075682507 #integrated with quad from 50 to 500
    uarr = h14sfr(tarr)/h14norm * nSN / (4*np.pi/3 * Rvir**3)
    ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def saveuIIIbetalogit(label,tarr,Vmixfn,tvir,tburst,u0III,numprocs=1):
    betafn = get_betafn(tvir,tburst,3)
    uarr = betafn(tarr)*u0III
    ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs,points=[tvir-tburst/2.,tvir+tburst/2.])
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0,tthresh,tstop):
    uarr,muarr = solveforuIII(tarr,u_init,Vmixfn,th,u0,tthresh,tstop,verbose=True)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = interp1d(tarr,muarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,muarr
    return fn

def get_h14sfrfn():
    import asciitable
    import scipy.interpolate as interpolate
    d = asciitable.read("SFRDATA/hirano14.csv")
    d.dtype.names = ('z','t','N')
    tint = np.concatenate(([0,50],d['t'][::-1],[500,1000]))
    Nint = np.concatenate(([0,0,],d['N'][::-1],[0,0]))
    f = interpolate.UnivariateSpline(tint,Nint)
    tmin = 55; tmax = 496 #found by trial and error to be close/above 0
    myfn = functools.partial(_h14sfrfn,f=f,tmin=tmin,tmax=tmax)
    return myfn
def _h14sfrfn(t,f,tmin,tmax):
    t = np.ravel(t)
    out = f(t)
    out[t<tmin] = 0
    out[t>tmax] = 0
    return out

def get_betafn(tvir,tburst,n):
    tmin = tvir-tburst/2.; tmax = tvir+tburst/2.
    norm = _betafnnorm(n,tmin,tmax)
    myfn = functools.partial(_betafn,n=n,tmin=tmin,tmax=tmax,norm=norm)
    return myfn
def _betafn(t,n,tmin,tmax,norm): #symmetrical beta function
    t = np.ravel(t)
    out = (t-tmin)**n * (tmax-t)**n
    out[t<tmin] = 0
    out[t>tmax] = 0
    return out/norm
def _betafnnorm(n,tmin,tmax): #found analytically in mathematica: n>0
    #return integrate.quad(lambda t: (t-tmin)**n * (tmax-t)**n,tmin,tmax)[0]
    return (2**(-1.-2*n) * (tmax-tmin)**(1.+2*n) * np.sqrt(np.pi) * scipy.special.gamma(1.+n))/scipy.special.gamma(1.5+n)
def get_logitfn(tvir,trec,trecfact=.25):
    return functools.partial(_logit,tvir=tvir,trec=trec,trecfact=trecfact)
def _logit(t,tvir,trec,trecfact):
    t = np.ravel(t)
    return 1./(1. + np.exp(-(t-(tvir+2*trec))/(trec*trecfact)))

def Qp(muIII):
    return np.exp(-muIII-(muIII**2)/4.)

def uII(Qparr,narr,u0=.4,nvir=.1):
    return u0 * (1-Qparr) * narr/nvir

def uIII(muarr,narr,u0=.04,nvir=.1):
    return u0 * Qp(muarr) * narr/nvir

def squarediff(a1,a2):
    return np.sqrt(np.sum((a1-a2)**2))

def solveforuIII(tarr,u_init,Vmixfn,th,u0,
                 tthresh,tstop,
                 itermax=30,threshold=10**-6,convergemin=1,
                 showplot=False,verbose=True,numprocs=1):
    assert len(tarr) == len(u_init)
    Vmixarr = Vmixfn(tarr)
    narr = (th.get_n_of_t_fn())(tarr)

    uarr = u_init
    convergecount = 0
    for i in xrange(itermax):
        if verbose: start = time.time()
        ufn = interp1d(tarr,uarr,bounds_error=False,fill_value=0)
        mufn = get_mufn(Vmixfn,ufn,tarr,numprocs=numprocs)
        muarr = mufn(tarr)
        unew = uIII(muarr,narr,u0=u0,nvir=th.nvir)
        unew[tarr<tthresh] = 0.
        unew[tarr>tstop]   = 0.

        #if verbose: print unew
        chi2 = squarediff(unew,uarr)/len(tarr)
        if verbose: print chi2,time.time()-start; sys.stdout.flush()
        uarr = unew
        if showplot:
            plotdebug(tarr,muarr,uarr,narr,Vmixarr,th)
        if chi2 < threshold: 
            if convergecount >= convergemin:
                if verbose:
                    print "solveforuIII converged in %i iterations" % (i+1)
                break
            else:
                convergecount += 1
        else:
            convergecount = 0
    if convergecount < convergemin:
        print "Warning: solveforuIII has not converged, stopping at %i iterations" % (itermax)

    muarr = mufn(tarr)
    return uarr,muarr

def plotdebug(tarr,muarr,uarr,narr,Vmixarr,th):
    import pylab as plt
    print uarr[99:110]
    f,axarr = plt.subplots(2,2)
    ((ax1,ax2),(ax3,ax4)) = axarr
    ax1.plot(tarr,muarr)
    ax1.set_ylabel('mu')
    ax2.plot(tarr,uarr,lw=2); ax2.plot(tarr,u(muarr,narr,nvir=th.nvir))
    ax2.set_ylabel('u'); ax2.set_yscale('log')
    ax3.plot(tarr,Qp(muarr))
    ax3.set_ylabel('Qp')#; ax3.set_yscale('log')
    #ax4.plot(tarr,narr)
    #ax4.set_ylabel('n')#; ax4.set_yscale('log')
    ax4.plot(tarr,Vmixarr)
    ax4.set_ylabel('Vmix')#; ax4.set_yscale('log')
    plt.show()

def get_sfrlabel(envname,sfrname):
    return envname+'__'+sfrname

def sfrparams(envname,sfrname,verbose=False):
    try:
        name,ts = sfrname.split('TS')
    except ValueError:
        raise ValueError("Invalid sfrname: "+sfrname+" (requires 'TS' to denote tstop)")
    ts = int(ts)

    Mhalo,zvir,vturb,lturb,nSN,trecovery,Dt,uSN,E51 = bgm.envparams(envname)
    u0III = uSN
    ttIII = 100
    if 'minihalo' in envname:
        ttII = ttIII + 10
    elif 'atomiccoolinghalo' in envname:
        ttII = ttIII + 300

    if 'burst' in sfrname:
        assert 'atomiccoolinghalo' in envname
        #ttIII = 100; ttII = 400; tburst = 50

    if 'h14' in sfrname:
        assert 'minihalo' in envname
        zvir = bgm.envparams(envname)[1]
        c = cosmology.cosmology(Ol=0.7,Om=0.3,h=0.7)
        u0III = -1; ttIII = -1 #not used
        ttII = c.t_age(zvir)
        u0II = 0.4
        #TODO I think this is the wrong way to proceed...
    if "fix" in sfrname:
        u0II = 0.4
    elif "10x" in sfrname:
        u0II = 10*u0III
    elif "mhsim" in sfrname:
        u0II = 2.0 #100Msun with TS150
    elif "achsim" in sfrname:
        u0II = 0.1 #10000Msun with TS500
    elif "achslo" in sfrname:
        u0II = 0.01 #1000Msun with TS500

    if 'u0IIMLT' in sfrname:
        pre,suf = sfrname.split('u0IIMLT')
        multfact = float(suf[0:2])
        u0II *= multfact
    if 'u0IIDIV' in sfrname:
        pre,suf = sfrname.split('u0IIDIV')
        multfact = float(suf[0:2])
        u0II /= multfact
    if 'u0IIIMLT' in sfrname:
        pre,suf = sfrname.split('u0IIIMLT')
        multfact = float(suf[0:2])
        u0III *= multfact
    if 'u0IIIDIV' in sfrname:
        pre,suf = sfrname.split('u0IIIDIV')
        multfact = float(suf[0:2])
        u0III /= multfact

    if 'betalogit' in sfrname:
        zvir = bgm.envparams(envname)[1]
        c = cosmology.cosmology(Ol=0.7,Om=0.3,h=0.7)
        ttIII = c.t_age(zvir)*1000. #tvir
        ttII = ttIII+trecovery
        # Normalize u0II based on total Pop II mass
        if 'minihalo' in envname: MpopII = 100.
        elif 'atomiccoolinghalo' in envname: MpopII = 10.**4
        numtomass = 184. #assuming salpeter slope:
        th,Vmixfn = util.load_Vmix(envname,get_th=True) #just for th
        Rvir = th.Rvir
        logitfn = get_logitfn(ttIII,trecovery)
        myint = integrate.quad(logitfn,0,ts)[0]
        u0II = MpopII/(myint*(4*np.pi/3.)*Rvir**3*numtomass)
        u0III = nSN/(4*np.pi/3 * Rvir**3)

    if verbose:
        print name
        print "u0II %3.2e u0III %3.2e" % (u0II,u0III)
        print "ttII %i ttIII %i" % (tII,ttIII)
        print "ts %i" % (ts)
    return u0II,u0III,ttII,ttIII,ts

def _get_mufn(t,Vmix,uSN,points):
    try:
        return integrate.quad(lambda tp: Vmix(t-tp)*uSN(tp),0,t,limit=1000,points=points)[0]
    except ValueError:
        return integrate.quad(lambda tp: Vmix(t-tp)*uSN(tp),0,t,limit=1000)[0]
def get_mufn(Vmix,uSN,tarr,numprocs=1,points=None):
    myfunc = functools.partial(_get_mufn,Vmix=Vmix,uSN=uSN,points=points)
    if numprocs==1:
        muarr = [myfunc(t) for t in tarr]
    else:
        pool = Pool(numprocs)
        muarr = pool.map(myfunc,tarr)
        pool.close()
    return interp1d(tarr,muarr)
#this doesn't work...? Need to read up on DFTs
#def get_mufn2(Vmix,uSN,tarr):
#    Vn = Vmix(tarr); un = uSN(tarr)
#    muarr = irfft(rfft(Vn)*rfft(un)) 
#    return interp1d(tarr,muarr)

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--plotmany', action='store_true',dest='plotmany', default=False)
    parser.add_option('--hires',action='store_true',dest='hires', default=False)
    options,args = parser.parse_args()

    envname = args[0]
    if options.plotmany:
        f = open('sfrnames.txt','r')
        sfrnames = [sfrname.strip() for sfrname in f.readlines()]
        f.close()
        print sfrnames
        for sfrname in sfrnames:
            plot_sfu(sfrlabel(envname,sfrname),showplot=False)
    else:
        sfrname = args[1]
        print "Running solveforu"
        print "envname:",envname
        print "sfrname:",sfrname
        
        tarr = bgm.gettarr(envname,hires=options.hires)
        start = time.time()
        compute_umufns(envname,sfrname,tarr)
        print "Finished:",time.time()-start
