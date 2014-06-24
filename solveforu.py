import numpy as np
from scipy import integrate,interpolate
from tophat import TopHat
import karlsson
from optparse import OptionParser

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
    #test()
    parser = OptionParser()
    parser.add_option('--uII', action='store_true',dest='uII', default=False)
    parser.add_option('--uIII',action='store_true',dest='uIII',default=False)
    parser.add_option('--tres',action='store',type='float',dest='tres',
                      default=.1)
    parser.add_option('--logmdil',action='store',type='float',dest='logMdil',
                      default=5)
    options,args = parser.parse_args()
    
    Mhalo,zvir,vturb,lturb = karlsson.params_k08()
    Dt,uSN = karlsson.get_Dt_uSN(vturb,lturb,0,1)

    tarr = np.arange(.01,1000,options.tres)
    if options.uIII:
        tthresh=100
        u0 = .04; th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
        Vmixfn = karlsson.get_Vmixfn_K08(th.get_rho_of_t_fn(),Dt=Dt,Mdil=10**options.logMdil)
        numsf = np.sum(tarr>tthresh); numskip=len(tarr)-numsf
        u_init = np.concatenate((np.zeros(numskip),np.logspace(-2,-8,numsf)))
        if options.logMdil==5:
            saveuIIIdata('atomiccoolinghalo',tarr,u_init,Vmixfn,th,u0,tthresh=tthresh)
        else:
            saveuIIIdata('atomiccoolinghalo'+str(options.logMdil),tarr,u_init,Vmixfn,th,u0,tthresh=tthresh)
    
    if options.uII:
        tthresh=110
        u0 = .4; th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
        Vmixfn = karlsson.get_Vmixfn_K08(th.get_rho_of_t_fn(),Dt=Dt,Mdil=10**options.logMdil)
        if options.logMdil==5:
            saveuIIdata('atomiccoolinghalo',Vmixfn,th,u0,tthresh=tthresh)
        else:
            saveuIIdata('atomiccoolinghalo'+str(options.logMdil),Vmixfn,th,u0,tthresh=tthresh)
