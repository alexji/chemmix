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
                          headernotes='',
                          XH=0.75,Nstars=10**5,numprocs=10):
    Mplot = np.load(karlsson.get_Mplot_filename(filename))
    masstonum = 1.0/(yieldobj.elemmass * XH)
    output = np.zeros((Nstars,yieldobj.numyields,kmax))
    sntypearr = np.arange(1,yieldobj.numtypes+1)

    pool = Pool(numprocs)
    print "Using",numprocs,"processors to run",filename,"with N =",Nstars
    for k in xrange(1,kmax+1):
        print "Starting k =",k
        start = time.time()
        fMk = np.load(karlsson.get_fMk_filename(filename,k))
        starmaker = functools.partial(run_one_star,k=k,Mplot=Mplot,fMk=fMk,
                                      sntypearr=sntypearr,sntypepdf=sntypepdf,
                                      yieldobj=yieldobj,masstonum=masstonum)
        this_output = pool.map(starmaker,xrange(Nstars))
        output[:,:,k-1] = np.array(this_output)
        print "  time:",time.time()-start
    pool.close()#; pool.join()

    postfix = '_'+yieldobj.shortname
    outfile = karlsson.get_chemgrid_filename(filename,postfix=postfix)
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
    n06y = yields.nomoto06interpyields()
    sntypepdf = relative_imf(n06y.massarr,2.35,norm=True) #salpeter imf
    #run_compute_random_sn('minihalo',sntypepdf,n06y,10,
    #                      headernotes='salpeter imf',
    #                      Nstars=10**6,numprocs=20)

    run_compute_random_sn('atomiccoolinghalo',sntypepdf,n06y,15,
                          headernotes='salpeter imf',
                          Nstars=10**6,numprocs=20)
