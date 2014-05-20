import numpy as np
import karlsson
import time
import pylab as plt
import h5py

from tophat import TopHat

def run_compute_fMk(filename,logMmin,logMmax,logdM,Mhalo,vturb,lturb,nSN,trecovery,kmin=1,kmax=20,saveplot=True,numprocs=1):
    """
    Assume a TopHat density function
    Mhalo: Msun
    vturb: km/s
    lturb: kpc
    nSN: number
    trecovery: Myr
    """
    vturb *= 3.16/3.08 * .001 #km/s to kpc/yr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    print filename
    print "Mhalo",Mhalo,"vturb",vturb,"lturb",lturb
    print "Dt",Dt,"uSN",uSN

    th = TopHat(Mhalo=Mhalo,nvir=0.1,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO,Dt=Dt)
    MMIXGRID_FILENAME = karlsson.get_mmixgrid_filename(filename)

    WISM = karlsson.wISM_K05
    PSI = lambda t: uSN
    #PSI = karlsson.uSN_G08_const
    MUFN = karlsson.get_mufn(VMIX,PSI)
    ## IMPORTANT: Mbins and Mplot
    Mbins,Mplot = karlsson.logMbinsMplot(logMmin,logMmax,logdM)
    #Mbins = 10**np.arange(2,8.01,.01)
    #Mplot = 10**(np.arange(2,8.0,.01)+.005)
    #print "len Mplot:",len(Mplot)
    print Mbins.shape,Mplot.shape
    
    # Retrieve DMList
    f = h5py.File(MMIXGRID_FILENAME)
    tarr = np.array(f['tarr']); tauarr = tarr
    Mmixgrid = np.array(f['Mmixgrid'])
    print "mmixgrid shape:",Mmixgrid.shape
    f.close()
    print "Starting get_DMlist with",numprocs,"processors"
    start = time.time()
    DMlist = karlsson.get_DMlist(Mbins,tarr,tauarr,Mmixgrid,numprocs=numprocs)
    print "get_DMlist time:",time.time()-start

    # Calculate and save fMk
    if saveplot:
        plt.figure()
    for k in range(kmin,kmax+1):
        filename_k = karlsson.get_fMk_filename(filename,k)
        start = time.time()
        fMk = karlsson.calc_fMk(k,Mbins,DMlist,VMIX,WISM,MUFN,PSI)
        print "calc_fM",k,time.time()-start,len(fMk)
        np.save(filename_k,fMk)
        #fMk = np.load(FILENAME+'_fM'+str(k)+'.npy')
        #print np.sum(fMk)
        if saveplot:
            plt.plot(Mplot,fMk,label=str(k))
    np.save(karlsson.get_Mbins_filename(filename),Mbins)
    np.save(karlsson.get_Mplot_filename(filename),Mplot)
    if saveplot:
        plt.savefig(filename+'_fMk.png',bbox_inches='tight')
        plt.show()

if __name__=="__main__":
    ## minihalo
    filename='minihalo'
    Mhalo,vturb,lturb,nSN,trecovery = karlsson.params_minihalo()
    run_compute_fMk(filename,0,7,.01,
                    Mhalo,vturb,lturb,nSN,trecovery,numprocs=1)

    ## atomic cooling halo
    filename='atomiccoolinghalo'
    Mhalo,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo()
    run_compute_fMk(filename,0,8,.01,
                    Mhalo,vturb,lturb,nSN,trecovery,numprocs=1)

#    fMk_foldername = "MMIXDISTR/"

    #kmin = 1
    #kmax = 20
    #NUMPROCS = 1
    ##tmin = 0.0; tmax = 1000.0; dt = 1.0
    #tmin = 0.0; tmax = 1000.0; dt = 0.03
    ##FILENAME = "TESTb_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    ##FILENAME = "TESTc_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    ##VMIX = karlsson.Vmix_K05
    ##RHO = karlsson.rho_K05A
    ##WISM = karlsson.wISM_K05
    ##MUFN = karlsson.get_mufn_K05A()
    ##PSI = karlsson.uSN_K05A
    ##Mbins = np.linspace(0,3,301)*10**6
    ##Mplot = (Mbins[:-1]+Mbins[1:])/2.0
    #FILENAME = "Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    #th = TopHat()
    #RHO = th.get_rho_of_t_fn()
    #VMIX = karlsson.get_Vmixfn_K08(RHO)
    #WISM = karlsson.wISM_K05
    #PSI = karlsson.uSN_G08_const
    #MUFN = karlsson.get_mufn(VMIX,PSI)
    #Mbins = 10**np.arange(2,8.01,.01)
    #Mplot = 10**(np.arange(2,8.0,.01)+.005)
    #print "len Mplot:",len(Mplot)
    #
    ## Retrieve DMList
    ##saved_files = np.load(FILENAME+".npz")
    #f = h5py.File(FILENAME+'.hdf5')
    #tarr = np.array(f['tarr']); tauarr = tarr
    #Mmixgrid = np.array(f['Mmixgrid'])
    #print Mmixgrid.shape
    #f.close()
    #start = time.time()
    #DMlist = karlsson.get_DMlist(Mbins,tarr,tauarr,Mmixgrid,numprocs=NUMPROCS)
    #print "get_DMlist time:",time.time()-start
#
#    # Calculate and save fMk
#    plt.figure()
#    for k in range(kmin,kmax+1):
#        start = time.time()
#        fMk = karlsson.calc_fMk(k,Mbins,DMlist,VMIX,WISM,MUFN,PSI)
#        print "calc_fM",k,time.time()-start,len(fMk)
#        np.save(FILENAME+'_fM'+str(k),fMk)
#        #fMk = np.load(FILENAME+'_fM'+str(k)+'.npy')
#        #print np.sum(fMk)
#        plt.plot(Mplot,fMk,label=str(k))
#    #leg = plt.legend(loc='best')
#    #for label in leg.get_texts():
#    #    label.set_fontsize('xx-small')
#    plt.savefig(FILENAME+'_fMk.png',bbox_inches='tight')
#    plt.show()

