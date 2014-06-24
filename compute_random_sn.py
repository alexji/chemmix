import numpy as np
import h5py

import time
from multiprocessing import Pool
import functools
import numpy.random as random

import karlsson
import yields

def relative_imf(marr,alpha,norm=False):
    Mmax = float(marr[-1])
    imfpdf = np.array([(M/Mmax)**(-alpha) for M in marr])
    if norm:
        imfpdf = imfpdf/np.sum(imfpdf)
    return imfpdf

def run_one_star(j,k,Mplot,fMk,sntypearr,sntypepdf,yieldobj,masstonum,gaussianprofile=False):
                 #,carbonenhanceprob=0.0,carbonenhancefactor=100.):
    # draw k mixing masses, SN types
    mixing_masses = karlsson.draw_from_distr(k,Mplot,fMk)
    sn_types = karlsson.draw_from_distr(k,sntypearr,sntypepdf)
    # get the k * numyields element production array
    yields = yieldobj(sn_types) #np.array([yieldobj(sn_type) for sn_type in sn_types])
    # with a certain probability, add carbon enhancement
    #if carbonenhanceprob > 0:
    #    yields[random.rand(k) < carbonenhanceprob, 0] *= carbonenhancefactor
    # weighted sum to get 1 * numyields, which goes into output
    if gaussianprofile:
        pass #TODO!!!
    else:
        weightedyields = masstonum*np.sum((yields.transpose()/mixing_masses),1).transpose()
    return weightedyields

def run_compute_random_sn(filename,sntypepdf,yieldobj,kmax,
                          headernotes='',postfix='',
                          Mmax=None,rhop2=False,
                          XH=0.75,Nstars=10**5,numprocs=10):
    Mplot = np.load(karlsson.get_Mplot_filename(filename,rhop2=rhop2))
    masstonum = 1.0/(yieldobj.elemmass * XH)
    output = np.zeros((Nstars,yieldobj.numyields,kmax))
    sntypearr = np.arange(1,yieldobj.numtypes+1)
    assert len(sntypepdf)==len(sntypearr) #consistency check
    if Mmax != None:
        keepii = Mplot <= Mmax
        badii  = Mplot > Mmax
        Mplot = np.concatenate((Mplot[keepii],[Mmax]))

    pool = Pool(numprocs)
    print "Using",numprocs,"processors to run",filename,"with N =",Nstars
    for k in xrange(1,kmax+1):
        print "Starting k =",k
        start = time.time()
        fMk = np.load(karlsson.get_fMk_filename(filename,k,rhop2=rhop2))
        if Mmax != None:
            badprob = np.sum(fMk[badii])
            fMk = np.concatenate((fMk[keepii],[badprob]))
        starmaker = functools.partial(run_one_star,k=k,Mplot=Mplot,fMk=fMk,
                                      sntypearr=sntypearr,sntypepdf=sntypepdf,
                                      yieldobj=yieldobj,masstonum=masstonum)
        this_output = pool.map(starmaker,xrange(Nstars))
        output[:,:,k-1] = np.array(this_output)
        print "  time:",time.time()-start
    pool.close()#; pool.join()

    postfix = '_'+yieldobj.shortname+postfix
    outfile = karlsson.get_chemgrid_filename(filename,postfix=postfix,rhop2=rhop2)
    print "Finished! Saving file to",outfile
    f = h5py.File(outfile,"w")
    f['chemgrid']=output
    f.attrs['yields'] = yieldobj.name
    f.attrs['Nstars'] = Nstars
    f.attrs['XH'] = XH
    f.attrs['kmax'] = kmax
    f.attrs['notes'] = headernotes
    f.close()
    #np.save(outfile,output)
    #np.save(FILENAME+"_chemgrid_N06_Cenhanced_p0.1_N"+str(Nstars),output)

if __name__=="__main__":
    reload(yields)
    RHOP2 = False
    NUMPROCS = 15
    fb = 0.1551

    i05n06 = yields.I05N06yields()
    p = 0.1; sntypepdf = np.array([p,1-p])
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,
    #                      postfix='_p0.1',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate_lowvturb',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo_lowmass',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)
    run_compute_random_sn('k08',sntypepdf,i05n06,15,Mmax=fb*10**8,
                          postfix='_p0.1Mmax',Nstars=10**6,numprocs=NUMPROCS)

    p = 0.5; sntypepdf = np.array([p,1-p])
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,
    #                      postfix='_p0.5',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate_lowvturb',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo_lowmass',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)
    run_compute_random_sn('k08',sntypepdf,i05n06,15,Mmax=fb*10**8,
                          postfix='__p0.5Mmax',Nstars=10**6,numprocs=NUMPROCS)

    p = 0.9; sntypepdf = np.array([p,1-p])
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('minihalo4',sntypepdf,i05n06,10,Mmax=fb*10**6,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo4',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,
    #                      postfix='_p0.9',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate_lowvturb',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo_lowmass',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    run_compute_random_sn('k08',sntypepdf,i05n06,15,Mmax=fb*10**8,
                          postfix='_p0.9Mmax',Nstars=10**6,numprocs=NUMPROCS)
    
    p = 0.95; sntypepdf = np.array([p,1-p])
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate_lowvturb',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo_lowmass',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)
    run_compute_random_sn('k08',sntypepdf,i05n06,15,Mmax=fb*10**8,
                          postfix='_p0.95Mmax',Nstars=10**6,numprocs=NUMPROCS)

    p = 0.99; sntypepdf = np.array([p,1-p])
    #run_compute_random_sn('atomiccoolinghalo',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghaloearly',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalolate_lowvturb',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo_lowmass',sntypepdf,i05n06,15,Mmax=fb*10**8,
    #                      postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)
    run_compute_random_sn('k08',sntypepdf,i05n06,15,Mmax=fb*10**8,
                          postfix='_p0.99Mmax',Nstars=10**6,numprocs=NUMPROCS)

    #n06y = yields.nomoto06interpyields()
    #Nsntypepdf = relative_imf(n06y.massarr,2.35,norm=True) #salpeter imf
    #
    #n06p1 = yields.nomoto06interpyields_Cenhance(p=0.1,f=100)
    #Csntypepdf1 = n06p1.modify_sntypepdf(Nsntypepdf)
    #
    #n06p2 = yields.nomoto06interpyields_Cadd(p=0.1,MC=0.8)
    #Csntypepdf2 = n06p2.modify_sntypepdf(Nsntypepdf)
    #
    #hw1 = yields.hw10yields(1.2,'S4',0)
    #Hsntypepdf = relative_imf(hw1.massarr,2.35,norm=True)
    #H135sntypepdf = relative_imf(hw1.massarr,1.35,norm=True)
    #Hflatsntypepdf = relative_imf(hw1.massarr,0.00,norm=True)
    #
    #nhw1 = yields.mixN06HW10yields(0.1,1.2,'S4',0)
    #NHWsntypepdf1 = nhw1.modify_sntypepdf(Nsntypepdf)
    #nhw2 = yields.mixN06HW10yields(0.5,1.2,'S4',0)
    #NHWsntypepdf2 = nhw2.modify_sntypepdf(Nsntypepdf)
    #nhw3 = yields.mixN06HW10yields(0.9,1.2,'S4',0)
    #NHWsntypepdf3 = nhw3.modify_sntypepdf(Nsntypepdf)

    #run_compute_random_sn('minihalo',Nsntypepdf,n06y,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',Nsntypepdf,n06y,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #run_compute_random_sn('minihalo',Csntypepdf1,n06p1,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',Csntypepdf1,n06p1,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #run_compute_random_sn('minihalo',Csntypepdf2,n06p2,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',Csntypepdf2,n06p2,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #run_compute_random_sn('minihalo',Hsntypepdf,hw1,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',Hsntypepdf,hw1,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #run_compute_random_sn('minihalo',H135sntypepdf,hw1,10,
    #                      headernotes='alpha = 1.35 imf',postfix='_a1.35',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',H135sntypepdf,hw1,15,
    #                      headernotes='alpha = 1.35 imf',postfix='_a1.35',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #
    #run_compute_random_sn('minihalo',Hflatsntypepdf,hw1,10,
    #                      headernotes='flat imf',postfix='_flat',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',Hflatsntypepdf,hw1,15,
    #                      headernotes='flat imf',postfix='_flat',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #run_compute_random_sn('minihalo',NHWsntypepdf1,nhw1,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',NHWsntypepdf1,nhw1,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #
    #run_compute_random_sn('minihalo',NHWsntypepdf2,nhw2,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',NHWsntypepdf2,nhw2,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #
    #run_compute_random_sn('minihalo',NHWsntypepdf3,nhw3,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',NHWsntypepdf3,nhw3,15,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
    #n = yields.nomoto06yields()
    #nimf = np.array([0,0,0,0,1,0,0])
    #run_compute_random_sn('minihalo',nimf,n,10,
    #                      headernotes='M=25',postfix='_M25',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    #run_compute_random_sn('atomiccoolinghalo',nimf,n,15,
    #                      headernotes='M=25',postfix='_M25',
    #                      Nstars=10**6,numprocs=NUMPROCS)
    
