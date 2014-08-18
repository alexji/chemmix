import numpy as np
import pylab as plt
import backgroundmodel as bgm
import util
import yields
from multiprocessing import Pool
import functools
import time
from plot_fMkk import plot_fMkk_grid, plot_fMkk_weighted
from scipy.integrate import trapz
import pickle
from montecarloyields import _montecarloyields
from plot_histlist import _plotone,calc_cfrac
from plot_util import plot1dhist

def _calc_ckk(tasknum,tarr,k49max,k51max,Vmix49,Vmix51,uII,u49,u51,muII,mu49,mu51,wII,w49,w51):
    if tasknum==0: return 0
    k49,k51 = divmod(tasknum,k51max+1)
    f = lambda t: w49(k49,mu49(t))*w51(k51,mu51(t))*wII(0,muII(t))*uII(t)
    y = f(tarr); y[np.isnan(y)] = 0
    return trapz(y,x=tarr)

def calc_ckk(numprocs,prefix,hires=False):
    print "Calculating ckk..."
    if hires: hiresstr="_hires"
    else: hiresstr = ""
    if 'Vm' in prefix: sfrname = 'betalogitVmaxTS200'
    else: sfrname = 'betalogitTS200'
    Vmix49 = util.load_Vmix('minihalo_lowE'+hiresstr)
    Vmix51 = util.load_Vmix('minihalo'+hiresstr)
    uII,u49,muII,mu49,wII,w49 = util.load_sfr('minihalo_lowE',sfrname)
    uII,u51,muII,mu51,wII,w51 = util.load_sfr('minihalo',sfrname)
    k49max=10;k51max=10
    tarr = bgm.gettarr('minihalo'+hiresstr)

    myfunc = functools.partial(_calc_ckk,tarr=tarr,k49max=k49max,k51max=k51max,
                               Vmix49=Vmix49,Vmix51=Vmix51,
                               uII=uII,u49=u49,u51=u51,
                               muII=muII,mu49=mu49,mu51=mu51,
                               wII=wII,w49=w49,w51=w51)
    pool = Pool(numprocs)
    cgrid = pool.map(myfunc,range((k49max+1)*(k51max+1)))
    pool.close()

    cgrid = np.array(cgrid).reshape((k49max+1,k51max+1))
    cgrid = cgrid/np.nansum(cgrid)
    np.save('MMIXDISTR/'+prefix+'_cgrid.npy',cgrid)
    print "Done!"

def plot_mh_ckk(prefix):
    k49max=10;k51max=10
    ckk = np.load('MMIXDISTR/'+prefix+'_cgrid.npy')
    plt.matshow(ckk)
    plt.colorbar()
    plt.xlabel(r'$k_{51}$')
    plt.ylabel(r'$k_{49}$')
    plt.savefig("PLOTS/"+prefix+'_ckk.png')

def _calc_fMkk49(DM,k49max,k51max,Vmix49,Vmix51,uII,u49,u51,muII,mu49,mu51,wII,w49,w51,dt):
    outarr = np.zeros((k49max+1,k51max+1,1))
    t = DM[:,0]; tau = DM[:,1]
    Vmix49arr = Vmix49(tau)
    mu49arr = mu49(t)
    mu51arr = mu51(t)
    lmssfr  = uII(t)
    u49arr  = u49(t-tau)

    muIIarr = muII(t)
    wIIarr  = wII(0,muIIarr)

    preprod49 = dt*dt*Vmix49arr*lmssfr*u49arr*wIIarr
    for k49 in range(k49max+1):
        if k49==0: outarr[k49,:,0]=0
        else: 
            for k51 in range(k51max+1):
                wprod49 = w49(k49-1,mu49arr)*w51(k51,mu51arr)
                outarr[k49,k51,0] = np.sum(preprod49*wprod49)
    return outarr
def _calc_fMkk51(DM,k49max,k51max,Vmix49,Vmix51,uII,u49,u51,muII,mu49,mu51,wII,w49,w51,dt):
    outarr = np.zeros((k49max+1,k51max+1,1))
    t = DM[:,0]; tau = DM[:,1]
    Vmix51arr = Vmix51(tau)
    mu49arr = mu49(t)
    mu51arr = mu51(t)
    lmssfr  = uII(t)
    u51arr  = u51(t-tau)

    muIIarr = muII(t)
    wIIarr  = wII(0,muIIarr)

    preprod51 = dt*dt*Vmix51arr*lmssfr*u51arr*wIIarr
    for k51 in range(k51max+1):
        if k51==0: outarr[:,k51,0]=0
        else:
            for k49 in range(k49max+1):
                wprod51 = w49(k49,mu49arr)*w51(k51-1,mu51arr)
                outarr[k49,k51,0] = np.sum(preprod51*wprod51)
    return outarr
def calc_fMkk(numprocs,prefix,hires=False):
    if hires: hiresstr="_hires"
    else: hiresstr = ""
    sfrname = 'betalogitTS200'
    #fMkkE51 = np.load(util.fname_fMkkIII('minihalo',sfrname))
    #fMkkE49 = np.load(util.fname_fMkkIII('minihalo_lowE',sfrname))
    Vmix49 = util.load_Vmix('minihalo_lowE'+hiresstr)
    Vmix51 = util.load_Vmix('minihalo'+hiresstr)
    uII,u49,muII,mu49,wII,w49 = util.load_sfr('minihalo_lowE',sfrname)
    uII,u51,muII,mu51,wII,w51 = util.load_sfr('minihalo',sfrname)
    
    k49max=10;k51max=10
    Mbins,Mplot = util.logMbinsMplot(2,6,.01)
    
    tarr = bgm.gettarr('minihalo'+hiresstr)
    dt = float(tarr[1]-tarr[0])

    pool = Pool(numprocs)
    print "Running E49"
    print "Getting DMlist..."; start = time.time()
    if hires:
        start2 = time.time()
        mmixgrid = util.load_mmixgrid_smem('minihalo_lowE'+hiresstr)
        DMlist = bgm.get_DMlist_smem(Mbins,tarr,mmixgrid,numprocs=numprocs)
        print "Time to load shared mmixgrid: %f" % (time.time()-start2)
    else:        
        mmixgrid = util.load_mmixgrid('minihalo_lowE'+hiresstr)
        DMlist = bgm.get_DMlist(Mbins,tarr,mmixgrid,numprocs=numprocs)
    print "Time: %f" % (time.time()-start)
    print "Calculating fMkk49..."; start = time.time()
    myfunc = functools.partial(_calc_fMkk49,
                               k49max=k49max,k51max=k51max,
                               Vmix49=Vmix49,Vmix51=Vmix51,
                               uII=uII,u49=u49,u51=u51,
                               muII=muII,mu49=mu49,mu51=mu51,
                               wII=wII,w49=w49,w51=w51,
                               dt=dt)
    outarrlist = pool.map(myfunc,DMlist)
    fMkk49 = np.concatenate(outarrlist,axis=2)
    print "Time: %f" % (time.time()-start)

    print "Running E51"
    print "Getting DMlist..."; start = time.time()
    if hires:
        start2 = time.time()
        mmixgrid = util.load_mmixgrid_smem('minihalo'+hiresstr)
        DMlist = bgm.get_DMlist_smem(Mbins,tarr,mmixgrid,numprocs=numprocs)
        print "Time to load shared mmixgrid: %f" % (time.time()-start2)
    else:        
        mmixgrid = util.load_mmixgrid('minihalo'+hiresstr)
        DMlist = bgm.get_DMlist(Mbins,tarr,mmixgrid,numprocs=numprocs)
    print "Time: %f" % (time.time()-start)
    print "Calculating fMkk51..."; start = time.time()
    myfunc = functools.partial(_calc_fMkk51,
                               k49max=k49max,k51max=k51max,
                               Vmix49=Vmix49,Vmix51=Vmix51,
                               uII=uII,u49=u49,u51=u51,
                               muII=muII,mu49=mu49,mu51=mu51,
                               wII=wII,w49=w49,w51=w51,
                               dt=dt)
    outarrlist = pool.map(myfunc,DMlist)
    fMkk51 = np.concatenate(outarrlist,axis=2)
    print "Time: %f" % (time.time()-start)


    print "normalizing distributions..."
    for k49 in range(k49max+1):
        for k51 in range(k51max+1):
            fMkk49[k49,k51,:] = fMkk49[k49,k51,:]/np.sum(fMkk49[k49,k51,:])
            fMkk51[k49,k51,:] = fMkk51[k49,k51,:]/np.sum(fMkk51[k49,k51,:])
    fMkk49[np.isnan(fMkk49)] = 0.
    fMkk51[np.isnan(fMkk51)] = 0.

    np.save('MMIXDISTR/'+prefix+'_fMkk49.npy',fMkk49)
    np.save('MMIXDISTR/'+prefix+'_fMkk51.npy',fMkk51)
    print "Done!"

def plot_mh_fMkk(prefix):
    k49max=10;k51max=10
    Mbins,Mplot = util.logMbinsMplot(2,6,.01)
    fMkk49 = np.load('MMIXDISTR/mh_fMkk49.npy')
    fMkk51 = np.load('MMIXDISTR/mh_fMkk51.npy')
    ckk = np.load('MMIXDISTR/'+prefix+'_cgrid.npy')
    Mmax = 10**6 * .1551
    plot_fMkk_grid('PLOTS/'+prefix+'_fMkk4951.png',
                   k49max,k51max,Mplot,fMkk49,fMkk51,Mmax,
                   (15,15),'Mmix','k_{49}','k_{51}',(10**4.2,10**5.5),(0,.1))
    plot_fMkk_weighted('PLOTS/'+prefix+'_wfMkk4951.png',
                       k49max,k51max,Mplot,fMkk49,fMkk51,Mmax,ckk,
                       'Mmix','k_{49}','k_{51}',(10**4.2,10**5.5),(0,.1))

def mh_montecarloyields(numprocs,N,prefix):
    print "Running monte carlo yields"
    fname = "CHEMGRIDS/"+prefix+'.list'
    binwidth = .1
    binlist = [np.arange(-10,0+binwidth,binwidth) for i in range(6)]
    
    y49 = yields.ryI05N06(1.0) #only carbon rich
    y51 = yields.ryI05N06(0.0) #only carbon normal
    XH=0.75

    k49max=10;k51max=10
    Mbins,Mplot = util.logMbinsMplot(2,6,.01)
    fMkk49 = np.load("MMIXDISTR/"+prefix+"_fMkk49.npy")
    fMkk51 = np.load("MMIXDISTR/"+prefix+"_fMkk51.npy")
    masstonum = 1.0/(y49.elemmass * XH)

    start = time.time()
    pool = Pool(numprocs)
    myfunc = functools.partial(_montecarloyields,k2max=k49max,k3max=k51max,
                               Mplot=Mplot,fMkkII=fMkk49,fMkkIII=fMkk51,
                               yII=y49,yIII=y51,masstonum=masstonum,Nstars=N,
                               binlist=binlist)
    dictlist = pool.map(myfunc,range((k49max+1)*(k51max)+1))
    pool.close()
    print "Finished! Time %f" % (time.time()-start)
    f = open(fname,'w')
    pickle.dump(dictlist,f)
    f.close()

def plot_mh_mcy(prefix):
    k49max=10;k51max=10
    fname = "CHEMGRIDS/"+prefix+'.list'
    f = open(fname,'r')
    histlist = pickle.load(f)
    f.close()
    _plotone("PLOTS/"+prefix+"_cfracgrid.png",k49max,k51max,histlist,0.5,
             kIIlabel='k_{49}',kIIIlabel='k_{51}')
    _plotone("PLOTS/"+prefix+"_cfracgridzoom.png",k49max,k51max,histlist,0.5,
             kIIlabel='k_{49}',kIIIlabel='k_{51}',xlim=(-7,-1))

    ckk = np.load('MMIXDISTR/'+prefix+'_cgrid.npy')
    h,x = histlist[1][(5,5)]
    wCfrac = np.zeros(len(h))
    wFe    = np.zeros(len(h))

    for tasknum,histdict in enumerate(histlist):
        k49,k51 = divmod(tasknum,k51max+1)
        if k49==0 and k51==0: continue
        FeH,Cfrac = calc_cfrac(histdict,CFecrit=0.75)
        h,x = histdict[(5,5)]
        wCfrac = wCfrac + ckk[k49,k51]*Cfrac
        wFe = wFe + ckk[k49,k51]*h
    wFe = wFe/np.sum(wFe)
    fig,ax = plt.subplots()
    plot1dhist(wFe,x,ax=ax,color='black')
    ax.plot(FeH,wCfrac)
    ax.plot([-4,-3],[0.75,0.3],'ko-')
    ax.set_xlim((-7,-1))
    ax.set_xlabel('[Fe/H]')
    plt.savefig('PLOTS/'+prefix+'_wcfrac.png',bbox_inches='tight')

if __name__=="__main__":
    prefix='mh'
    numprocs=8
    #calc_ckk(numprocs,prefix,hires=False)
    #plot_mh_ckk(prefix)
    #calc_fMkk(numprocs,prefix,hires=False)
    #plot_mh_fMkk(prefix)
    #mh_montecarloyields(numprocs,10**5,prefix)
    plot_mh_mcy(prefix)
