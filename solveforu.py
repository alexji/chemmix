import numpy as np
from scipy import integrate
#from scipy.interpolate import interp1d
from tophat import TopHat
import karlsson
from optparse import OptionParser
import util
import time
import functools
    
def compute_umufns(envname,sfrname,tmin=.01,tmax=1000.,tres=.1):
    th,Vmixfn = util.load_Vmix(envname,get_th=True)
    label = util.get_sfrlabel(envname,sfrname)
    u0II,u0III,ttII,ttIII,tstop = util.get_sfrparams(sfrname,envname=envname)
    tarr = np.arange(tmin,tmax,tres)

    if 'flat' in sfrname:
        uIII = functools.partial(flattemplate,u0=u0III,tstart=ttIII,tstop=ttII)
        uII  = functools.partial(flattemplate,u0=u0II,tstart=ttII,tstop=tstop)
        print "Getting mufns for "+sfrname+"..."; start = time.time()
        muIII = karlsson.get_mufn(Vmixfn,uIII,tarr=tarr)
        muII  = karlsson.get_mufn(Vmixfn,uII,tarr=tarr)
        print "Took %f" % (time.time()-start)
        uIIarr  = uII(tarr);  muIIarr  = muII(tarr)
        uIIIarr = uIII(tarr); muIIIarr = muIII(tarr)
        np.savez('SFRDATA/'+label+'_uII.npz',tarr=tarr, uarr=uIIarr, muarr=muIIarr)
        np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uIIIarr,muarr=muIIIarr)
    else:
        numsf = np.sum(tarr>ttIII); numskip = len(tarr)-numsf
        u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
        saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0III,ttIII,tstop)
        saveuIIdata(label,Vmixfn,th,u0II,ttII,tstop)

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

def saveuIIdata(label,Vmixfn,th,u0,tthresh,tstop):
    uIIIfn,tarr,uarr = loaduIIIfn(label,retarrays=True)
    muIIIfn = loadmuIIIfn(label)#karlsson.get_mufn(Vmixfn,uIIIfn,tarr=tarr)
    Qparr = Qp(muIIIfn(tarr))
    narr = (th.get_n_of_t_fn())(tarr)
    uarr = uII(Qparr,narr,u0=u0,nvir=th.nvir)
    uarr[tarr < tthresh] = 0.
    uarr[tarr > tstop]   = 0.

    ufn = util.interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    mufn = karlsson.get_mufn(Vmixfn,ufn,tarr=tarr)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = util.interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = util.interp1d(tarr,muarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,muarr
    return fn

def saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0,tthresh,tstop):
    uarr,muarr = solveforuIII(tarr,u_init,Vmixfn,th,u0,tthresh,tstop,verbose=True)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = util.interp1d(tarr,uarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = util.interp1d(tarr,muarr,bounds_error=False,fill_value=0)
    if retarrays: return fn,tarr,muarr
    return fn

def flattemplate(t,u0,tstart,tstop):
    try:
        if t>=tstart and t<=tstop: return u0
        return 0.
    except:
        outarr = np.zeros(len(t))
        outarr[np.logical_and(t>=tstart, t<=tstop)] = u0
        return outarr

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
                 itermax=10,threshold=10**-6,convergemin=1,
                 showplot=False,verbose=True):
    assert len(tarr) == len(u_init)
    Vmixarr = Vmixfn(tarr)
    narr = (th.get_n_of_t_fn())(tarr)

    uarr = u_init
    convergecount = 0
    for i in xrange(itermax):
        if verbose: start = time.time()
        ufn = util.interp1d(tarr,uarr,bounds_error=False,fill_value=0)
        mufn = karlsson.get_mufn(Vmixfn,ufn,
                                 tarr=tarr)
        muarr = mufn(tarr)
        unew = uIII(muarr,narr,u0=u0,nvir=th.nvir)
        unew[tarr<tthresh] = 0.
        unew[tarr>tstop]   = 0.

        #if verbose: print unew
        chi2 = squarediff(unew,uarr)/len(tarr)
        if verbose: print chi2,time.time()-start
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

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--plotmany', action='store_true',dest='plotmany', default=False)
    parser.add_option('--tres',action='store',type='float',dest='tres',
                      default=.1)
    options,args = parser.parse_args()

    envname = args[0]
    if options.plotmany:
        f = open('sfrnames.txt','r')
        sfrnames = [sfrname.strip() for sfrname in f.readlines()]
        f.close()
        print sfrnames
        for sfrname in sfrnames:
            plot_sfu(util.get_sfrlabel(envname,sfrname),showplot=False)
    else:
        sfrname = args[1]
        print "Running solveforu"
        print "envname:",envname
        print "sfrname:",sfrname
        
        start = time.time()
        compute_umufns(envname,sfrname,tres=options.tres)
        print "Finished:",time.time()-start

#    label = 'atomiccoolinghalo'
#    if options.double: label = 'atomiccoolinghalodouble'
#    elif options.doubIII: label = 'atomiccoolinghalodoubIII'
#    elif options.ten: label = 'atomiccoolinghaloten'
#    elif options.tenIII: label = 'atomiccoolinghalotenIII'
#    elif options.tenth: label = 'atomiccoolinghalotenth'
#    elif options.tenthIII: label = 'atomiccoolinghalotenthIII'
#    if options.logMdil!=5: label += str(options.logMdil)
#
#    print label
#    if options.uIII:
#        tthresh=100
#        u0 = .04; th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
#        if options.double: u0 = u0*2
#        elif options.doubIII: u0 = u0*2
#        elif options.tenIII: u0 = u0*10
#        elif options.ten: u0 = u0*10
#        elif options.tenth: u0 = u0*0.1
#        elif options.tenthIII: u0 = u0*0.1
#        Vmixfn = karlsson.get_Vmixfn_K08(th.get_rho_of_t_fn(),Dt=Dt,Mdil=10**options.logMdil)
#        numsf = np.sum(tarr>tthresh); numskip=len(tarr)-numsf
#        u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
#        saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0,tthresh=tthresh)
#    
#    if options.uII:
#        tthresh=110
#        u0 = .4; th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
#        if options.double: u0 = u0*2
#        elif options.ten: u0 = u0*10
#        elif options.tenth: u0 = u0*0.1
#        Vmixfn = karlsson.get_Vmixfn_K08(th.get_rho_of_t_fn(),Dt=Dt,Mdil=10**options.logMdil)
#        saveuIIdata(label,Vmixfn,th,u0,tthresh=tthresh)
#    
#
#    if options.plotmany:
#        for suffix in ['','double','doubIII','ten','tenIII','tenth','tenthIII']:
#            plot_sfu('atomiccoolinghalo'+suffix,showplot=False)
