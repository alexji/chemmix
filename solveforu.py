import numpy as np
from scipy import integrate,interpolate
from tophat import TopHat
import karlsson
from optparse import OptionParser
import util
import time
    
def compute_umufns(envname,sfrname,tmin=.01,tmax=1000.,tres=.1):
    th,Vmixfn = util.load_Vmix(envname,get_th=True)
    label = util.get_sfrlabel(envname,sfrname)
    u0II,u0III,ttII,ttIII,tstop = util.get_sfrparams(sfrname)
    tarr = np.arange(tmin,tmax,tres)

    numsf = np.sum(tarr>ttIII); numskip = len(tarr)-numsf
    u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
    saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0III,tthresh=ttIII)
    saveuIIdata(label,Vmixfn,th,u0II,tthresh=ttII)

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

def saveuIIdata(label,Vmixfn,th,u0,tthresh=110):
    uIIIfn,tarr,uarr = loaduIIIfn(label,retarrays=True)
    muIIIfn = loadmuIIIfn(label)#karlsson.get_mufn(Vmixfn,uIIIfn,tarr=tarr)
    Qparr = Qp(muIIIfn(tarr))
    narr = (th.get_n_of_t_fn())(tarr)
    uarr = uII(Qparr,narr,u0=u0,nvir=th.nvir,tarr=tarr,tthresh=tthresh)

    ufn = interpolate.InterpolatedUnivariateSpline(tarr,uarr)
    mufn = karlsson.get_mufn(Vmixfn,ufn,tarr=tarr)
    muarr = mufn(tarr)
    np.savez('SFRDATA/'+label+'_uII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = interpolate.InterpolatedUnivariateSpline(tarr,uarr)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = interpolate.InterpolatedUnivariateSpline(tarr,muarr)
    if retarrays: return fn,tarr,muarr
    return fn

def saveuIIIdata(label,tarr,u_init,Vmixfn,th,u0,tthresh=100):
    uarr,muarr = solveforuIII(tarr,u_init,Vmixfn,th,u0,tthresh,verbose=True)
    np.savez('SFRDATA/'+label+'_uIII.npz',tarr=tarr,uarr=uarr,muarr=muarr)
def loaduIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    uarr = npz['uarr']; tarr = npz['tarr']
    fn = interpolate.InterpolatedUnivariateSpline(tarr,uarr)
    if retarrays: return fn,tarr,uarr
    return fn
def loadmuIIIfn(label,retarrays=False):
    npz = np.load('SFRDATA/'+label+'_uIII.npz')
    muarr = npz['muarr']; tarr = npz['tarr']
    fn = interpolate.InterpolatedUnivariateSpline(tarr,muarr)
    if retarrays: return fn,tarr,muarr
    return fn

def Qp(muIII):
    return np.exp(-muIII-(muIII**2)/4.)

def uII(Qparr,narr,u0=.4,nvir=.1,tarr=None,tthresh=100):
    uarr = u0 * (1-Qparr) * narr/nvir
    if tarr != None:
        uarr[tarr<tthresh] = 0
    return uarr

def uIII(muarr,narr,u0=.04,nvir=.1,tarr=None,tthresh=100):
    uarr = u0 * Qp(muarr) * narr/nvir
    if tarr != None:
        uarr[tarr<tthresh] = 0
    return uarr

def squarediff(a1,a2):
    return np.sqrt(np.sum((a1-a2)**2))

def solveforuIII(tarr,u_init,Vmixfn,th,u0,
                 tthresh,
                 itermax=10,threshold=10**-6,convergemin=1,
                 showplot=False,verbose=True):
    assert len(tarr) == len(u_init)
    Vmixarr = Vmixfn(tarr)
    narr = (th.get_n_of_t_fn())(tarr)

    uarr = u_init
    convergecount = 0
    for i in xrange(itermax):
        ufn = interpolate.InterpolatedUnivariateSpline(tarr,uarr)
        mufn = karlsson.get_mufn(Vmixfn,ufn,
                                 tarr=tarr)
        muarr = mufn(tarr)
        unew = uIII(muarr,narr,u0=u0,nvir=th.nvir,
                    tarr=tarr,tthresh=tthresh)
        print unew
        chi2 = squarediff(unew,uarr)/len(tarr)
        if verbose: print chi2
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

def test():
    import pylab as plt
    print "Testing solveforu"
    u0 = .04
    th = TopHat(Mhalo=10**8,zvir=10)
    Vmixfn = karlsson.get_Vmixfn_K08(th.get_rho_of_t_fn())
    print "nvir:",th.nvir

    tarr = np.arange(0.01,1000,1)
    numsf = np.sum(tarr>100); numskip=len(tarr)-numsf
    u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
    
    uarr,muarr = solveforuIII(tarr,u_init,Vmixfn,th,u0)
    plt.plot(tarr,uarr)
    plt.loglog(); plt.ylabel('u'); plt.xlabel('t')
    plt.xlim((100,1000)); plt.ylim((10**-8,10**0))
    plt.show()

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option('--plotmany', action='store_true',dest='plotmany', default=False)

#    parser.add_option('--uII', action='store_true',dest='uII', default=False)
#    parser.add_option('--uIII',action='store_true',dest='uIII',default=False)
#    parser.add_option('--double',action='store_true',dest='double',default=False)
#    parser.add_option('--doubIII',action='store_true',dest='doubIII',default=False)
#    parser.add_option('--ten',action='store_true',dest='ten',default=False)
#    parser.add_option('--tenIII',action='store_true',dest='tenIII',default=False)
#    parser.add_option('--tenth',action='store_true',dest='tenth',default=False)
#    parser.add_option('--tenthIII',action='store_true',dest='tenthIII',default=False)

    parser.add_option('--tres',action='store',type='float',dest='tres',
                      default=.1)
#    parser.add_option('--logmdil',action='store',type='int',dest='logMdil',
#                      default=5)
    options,args = parser.parse_args()

    print "Running solveforu"
    print "envname:",args[0]
    print "sfrname:",args[1]

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
