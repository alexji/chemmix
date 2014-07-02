import numpy as np
import numpy.random as random
import time
from multiprocessing import Pool
from optparse import OptionParser
import functools

import karlsson
import util

def _run_compute_chemgrid(tasknum,kmax,Mplot,fMkkp,yieldII,yieldIII,masstonum,Nstars,
                          timethis=True):
    if timethis: start = time.time()
    #assert yieldII.numyields == yieldIII.numyields
    numyields = yieldII.numyields

    ik,ikp = divmod(tasknum,kmax+1); k=ik+1; kp=ikp
    if kp>k: return np.zeros((Nstars,numyields))
    this_fMkkp = fMkkp[ik,ikp,:]

    yieldmassratioIII = None; yieldmassratioII = None
    if kp > 0:
        mixing_massesIII= karlsson.draw_from_distr(kp*Nstars,Mplot,this_fMkkp)
        yieldarrIII= yieldIII.draw_yields(kp*Nstars)
        yieldmassratioIII= np.reshape((yieldarrIII.T/ mixing_massesIII).T,(Nstars,(kp),numyields))
    if k-kp > 0:
        mixing_massesII = karlsson.draw_from_distr((k-kp)*Nstars,Mplot,this_fMkkp)
        yieldarrII = yieldII.draw_yields((k-kp)*Nstars)
        yieldmassratioII = np.reshape((yieldarrII.T / mixing_massesII).T, (Nstars,(k-kp),numyields))

    #yieldmassratio has shape (Nstars,k,numyields)
    if yieldmassratioII != None and yieldmassratioIII != None:
        yieldmassratio = np.concatenate((yieldmassratioII,yieldmassratioIII),axis=1)
    elif yieldmassratioII == None:
        yieldmassratio = yieldmassratioIII
    elif yieldmassratioIII == None:
        yieldmassratio = yieldmassratioII
    else:
        raise ValueError("Somehow invalid k,kp: k=%i, kp=%i" % (k,kp))

    #yieldmassratio has shape (Nstars,numyields)
    yieldmassratio = np.sum(yieldmassratio,axis=1)
    yieldnumratio  = yieldmassratio*masstonum
    if timethis: print "task %i (k=%i, kp=%i): %f" % (tasknum,k,kp,time.time()-start)
    return yieldnumratio

def run_compute_chemgrid(envname,sfrname,yII,yIII,
                         Mmax,
                         Nstars=10**6,numprocs=1,
                         postfix='',XH=0.75):
    assert yII.elemnames == yIII.elemnames
    kmax,ckkp = util.load_ckkp(envname,sfrname)
    masstonum = 1.0/(yII.elemmass * XH)

    Mplot = np.load(karlsson.get_Mplotkkp_filename(envname,sfrname))
    keepii = Mplot <= Mmax
    badii  = Mplot > Mmax
    Mplot = np.concatenate((Mplot[keepii],[Mmax]))
    fMkkp = np.load(karlsson.get_fMkkp_filename(envname,sfrname))
    badprob = np.sum(fMkkp[:,:,badii],axis=2).reshape((fMkkp.shape[0],fMkkp.shape[1],1))
    fMkkp = np.concatenate((fMkkp[:,:,keepii],badprob),axis=2)

    
    chemgrid_filename = karlsson.get_chemgridkkp_filename(envname,sfrname,postfix=postfix)
    pool = Pool(numprocs)
    print "Using",numprocs,"processors to run",chemgrid_filename,"with N =",Nstars
    print "yII:",str(yII)
    print "yIII:",str(yIII)
    numyields = yII.numyields
    assert numyields == yIII.numyields

    start = time.time()
    this_func = functools.partial(_run_compute_chemgrid,kmax=kmax,Mplot=Mplot,fMkkp=fMkkp,
                                  yieldII=yII,yieldIII=yIII,masstonum=masstonum,
                                  Nstars=Nstars,timethis=False)
    numtasks = kmax*(kmax+1)
    chemlist = pool.map(this_func,range(numtasks))
    pool.close()
    print "Finished! Time: %f" % (time.time()-start)

    #copy into an array for saving
    start = time.time()
    output = np.zeros((Nstars,numyields,kmax,kmax+1))
    for i,arr in enumerate(chemlist):
        ik,ikp = divmod(i,kmax+1)
        output[:,:,ik,ikp] = arr
    np.save(chemgrid_filename,output)
    print "Time to copy/save grid: %f" % (time.time()-start)

def relative_imf(marr,alpha,norm=False):
    Mmax = float(marr[-1])
    imfpdf = np.array([(M/Mmax)**(-alpha) for M in marr])
    if norm:
        imfpdf = imfpdf/np.sum(imfpdf)
    return imfpdf

if __name__=="__main__":
    import yields
    parser = OptionParser()
    parser.add_option('-j','--numprocs',action='store',type='int',dest='numprocs',default=1)
    options,args = parser.parse_args()

    #envname = args[0]; sfrname = args[1]
    envname = 'atomiccoolinghalo'
    sfrname = 'fidTS400'
    postfix = 'test'
    fb = .1551; Mhalo = karlsson.envparams(envname)[0]
    Mmax = fb * Mhalo

    salimf  = relative_imf(np.array([13,15,18,20,25,30,40.]),2.35,norm=True) #salpeter imf
    popIIy  = yields.ryN06(salimf)
    
    popIIIy_arr = [yields.ryI05N06(p) for p in [.5]] #[.1,.25,.5,.75,.9,.95]]
    
    for popIIIy in popIIIy_arr:
        run_compute_chemgrid(envname,sfrname,popIIy,popIIIy,Mmax,Nstars=10**4,
                             numprocs=options.numprocs,postfix=postfix)
