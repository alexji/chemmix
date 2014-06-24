import numpy as np
import karlsson
import time
import pylab as plt
import h5py
from optparse import OptionParser

from tophat import TopHat

def run_compute_fMk(filename,logMmin,logMmax,logdM,Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,rhop2=False,kmin=1,kmax=20,saveplot=False,numprocs=1,k08=False):
    """
    Assume a TopHat density function
    Mhalo: Msun
    vturb: km/s
    lturb: kpc
    nSN: number
    trecovery: Myr
    """
    vturb *= 3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    print filename
    if (rhop2): print "psiLMS ~ rho"
    print "Mhalo",Mhalo,"zvir",zvir
    print "vturb",vturb,"lturb",lturb
    print "Dt",Dt,"uSN",uSN

    th = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
    RHO = th.get_rho_of_t_fn()
    VMIX = karlsson.get_Vmixfn_K08(RHO,Dt=Dt,Mdil=10**logMdil)
    MMIXGRID_FILENAME = karlsson.get_mmixgrid_filename(filename)
    
    WISM = karlsson.wISM_K05
    PSI = lambda t: uSN
    MUFN = karlsson.get_mufn(VMIX,PSI)
    Mbins,Mplot = karlsson.logMbinsMplot(logMmin,logMmax,logdM)
    #print Mbins.shape,Mplot.shape
    
    if (k08): 
        import solveforu as sfu
        print "Using Karlsson+08 SFR only looking at stars w/o Pop II enrichment"
        print "(ignoring uSN here)"
        print MMIXGRID_FILENAME
        filename = 'k08'
        if logMdil != 5: filename += str(logMdil)
        WISM = karlsson.wISM_III
        PSI    = sfu.loaduIIIfn('atomiccoolinghalo')
        PSIlms = sfu.loaduIIfn('atomiccoolinghalo')
        MUFN   = sfu.loadmuIIIfn('atomiccoolinghalo')
        MUFNlms= sfu.loadmuIIfn('atomiccoolinghalo')
        #WISMlms_0= lambda mu: karlsson.wISM_K05(0,mu) #e^-mu
        #MUFN    = karlsson.get_mufn(VMIX,PSI)
        #MUFNlms = karlsson.get_mufn(VMIX,PSIlms)

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
        filename_k = karlsson.get_fMk_filename(filename,k,rhop2)
        start = time.time()
        if rhop2:
            fMk = karlsson.calc_fMk(k,Mbins,DMlist,VMIX,WISM,MUFN,PSI,
                                    SFRlms=RHO)
        elif k08:
            fMk = karlsson.calc_fMk(k,Mbins,DMlist,VMIX,WISM,MUFN,PSI,
                                    SFRlms=PSIlms)#,WISMlms_0=WISMlms_0)
        else:
            fMk = karlsson.calc_fMk(k,Mbins,DMlist,VMIX,WISM,MUFN,PSI)
        print "calc_fM",k,time.time()-start,len(fMk)
        np.save(filename_k,fMk)
        #fMk = np.load(FILENAME+'_fM'+str(k)+'.npy')
        #print np.sum(fMk)
        if saveplot:
            plt.plot(Mplot,fMk,label=str(k))
    np.save(karlsson.get_Mbins_filename(filename,rhop2=rhop2),Mbins)
    np.save(karlsson.get_Mplot_filename(filename,rhop2=rhop2),Mplot)
    if saveplot:
        plt.savefig('PLOTS/'+filename+'_fMk.png',bbox_inches='tight')
        plt.show()

if __name__=="__main__":
    parser = OptionParser()
    parser.add_option("--k08",action='store_true',dest='k08',default=False)
    parser.add_option("--minihalo",action='store_true',dest='minihalo',default=False)
    parser.add_option("--atomiccoolinghalo",action='store_true',dest='atomiccoolinghalo',default=False)
    parser.add_option("--atomiccoolinghaloearly",action='store_true',dest='atomiccoolinghaloearly',default=False)
    parser.add_option("--atomiccoolinghalolate", action='store_true',dest='atomiccoolinghalolate', default=False)
    parser.add_option("--atomiccoolinghalolatelowvturb",action='store_true',dest='atomiccoolinghalolate_lowvturb',default=False)
    parser.add_option("--atomiccoolinghalolowmass",action='store_true',dest='atomiccoolinghalo_lowmass',default=False)
    parser.add_option("--logmdil",action='store',type='int',dest='logMdil',default=5)
    parser.add_option("--rhop2",action="store_true",dest='rhop2',default=False)
    parser.add_option("--save",action="store_true",dest='saveplot',default=False)
    parser.add_option("--numprocs","-j",action='store',type='int',dest='numprocs',default=1)
    options,args = parser.parse_args()
    logMdil = options.logMdil

    ## Compute Karlsson
    if options.k08:
        filename='lores_atomiccoolinghalo' #uses same mmixgrid; fMk will be different
        Mhalo,zvir,vturb,lturb = karlsson.params_k08()
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,-1,-1,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs,
                        k08=True)

    ## minihalo
    if options.minihalo:
        filename='minihalo'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_minihalo()
        run_compute_fMk(filename,0,7,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)

    ## atomic cooling halo
    if options.atomiccoolinghalo:
        filename='atomiccoolinghalo'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)
    if options.atomiccoolinghaloearly:
        filename='atomiccoolinghaloearly'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghaloearly()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)
    if options.atomiccoolinghalolate:
        filename='atomiccoolinghalolate'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalolate()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)
    if options.atomiccoolinghalolate_lowvturb:
        filename='atomiccoolinghalolate_lowvturb'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalolate_lowvturb()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)
    if options.atomiccoolinghalo_lowmass:
        filename='atomiccoolinghalo_lowmass'
        if logMdil != 5: filename += str(logMdil)
        Mhalo,zvir,vturb,lturb,nSN,trecovery = karlsson.params_atomiccoolinghalo_lowmass()
        run_compute_fMk(filename,0,8,.01,
                        Mhalo,zvir,vturb,lturb,nSN,trecovery,logMdil,
                        options.rhop2,saveplot=options.saveplot,numprocs=options.numprocs)
