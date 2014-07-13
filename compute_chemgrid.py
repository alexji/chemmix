import numpy as np
import numpy.random as random
import time
from multiprocessing import Pool
from optparse import OptionParser
import functools
import sys

import karlsson
import util
import gc
import pickle
import subprocess
import resource

def _get_tmpdictname(tasknum,tmplabel):
    return 'CHEMGRIDS/'+tmplabel+'.'+str(tasknum)+'.tmp'

def _run_compute_chemgrid(tasknum,kmax,kpmax,Mplot,fMkkp,yieldII,yieldIII,masstonum,Nstars,
                          timethis=True,
                          weighthist=False,binlist=None,tmplabel=''):
    gc.collect()
    if tmplabel != '': savetmpfiles=True
    if timethis: start = time.time()
    numyields = yieldII.numyields

    ik,ikp = divmod(tasknum,kpmax+1); k=ik+1; kp=ikp
    if kp>k: return np.zeros((Nstars,numyields))
    this_fMkkp = fMkkp[ik,ikp,:]

    yieldratioIII = None; yieldratioII = None
    if kp > 0:
        mixing_massesIII= karlsson.draw_from_distr(kp*Nstars,Mplot,this_fMkkp)
        yieldarrIII= yieldIII.draw_yields(kp*Nstars)
        yieldratioIII= np.reshape((yieldarrIII.T/ mixing_massesIII).T,(Nstars,(kp),numyields))
    if k-kp > 0:
        mixing_massesII = karlsson.draw_from_distr((k-kp)*Nstars,Mplot,this_fMkkp)
        yieldarrII = yieldII.draw_yields((k-kp)*Nstars)
        yieldratioII = np.reshape((yieldarrII.T / mixing_massesII).T, (Nstars,(k-kp),numyields))

    #yieldratio has shape (Nstars,k,numyields)
    if yieldratioII != None and yieldratioIII != None:
        yieldratio = np.concatenate((yieldratioII,yieldratioIII),axis=1)
    elif yieldratioII == None:
        yieldratio = yieldratioIII
    elif yieldratioIII == None:
        yieldratio = yieldratioII
    else:
        raise ValueError("Somehow invalid k,kp: k=%i, kp=%i" % (k,kp))
    yieldratioII = None; yieldratioIII = None
    gc.collect()

    #yieldratio has shape (Nstars,numyields)
    yieldratio = np.sum(yieldratio,axis=1) #mass ratios
    yieldratio  = yieldratio*masstonum     #number ratios
    if timethis: print "task %i (k=%i, kp=%i): %f" % (tasknum,k,kp,time.time()-start); sys.stdout.flush()

    #compute the histogram if weighthist is specified (to save memory)
    if weighthist:
        if timethis: start = time.time()
        elemnames = yieldII.elemnames
        yieldratio = karlsson.convert_to_solar(elemnames,yieldratio,verbose=False)        
        outdict = karlsson.hist_chemgridkkp_one(yieldratio,binlist,elemnames)
        if timethis: print "task %i histogram time: %f" % (tasknum,time.time()-start); sys.stdout.flush()
        if timethis: print "task %i memory: %f MB" % (tasknum, resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1000)
        if savetmpfiles:
            f = open(_get_tmpdictname(tasknum,tmplabel),'w')
            pickle.dump(outdict,f)
            f.close()
            return -1
        else:
            return outdict
    else:
        return yieldratio

def run_compute_chemgrid(envname,sfrname,yII,yIII,
                         Mmax,savegrid=False,
                         Nstars=10**6,numprocs=1,
                         postfix='',XH=0.75,
                         weighthist=False,binlist=None,
                         savetmpfiles=False):
    assert yII.elemnames == yIII.elemnames
    kmax,kpmax,ckkp = util.load_ckkp(envname,sfrname)
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
    if not savegrid: print "(will not save chemgrid to save space)"
    print "yII:",str(yII)
    print "yIII:",str(yIII)
    numyields = yII.numyields
    assert numyields == yIII.numyields

    start = time.time()
    if savetmpfiles:
        tmplabel = '_'+envname+'_'+sfrname
        this_func = functools.partial(_run_compute_chemgrid,kmax=kmax,kpmax=kpmax,
                                      Mplot=Mplot,fMkkp=fMkkp,
                                      yieldII=yII,yieldIII=yIII,masstonum=masstonum,
                                      Nstars=Nstars,
                                      timethis=True,
                                      weighthist=weighthist,binlist=binlist,
                                      tmplabel=tmplabel)
    else:
        this_func = functools.partial(_run_compute_chemgrid,kmax=kmax,kpmax=kpmax,
                                      Mplot=Mplot,fMkkp=fMkkp,
                                      yieldII=yII,yieldIII=yIII,masstonum=masstonum,
                                      Nstars=Nstars,
                                      timethis=True,
                                      weighthist=weighthist,binlist=binlist)
    numtasks = kmax*(kpmax+1)
    chemlist = pool.map(this_func,xrange(numtasks))
    pool.close()
    print "Finished! Time: %f" % (time.time()-start)

    if not weighthist:
        #copy into an array for saving/output
        start = time.time()
        output = np.zeros((Nstars,numyields,kmax,kpmax+1))
        for i,arr in enumerate(chemlist):
            ik,ikp = divmod(i,kpmax+1)
            output[:,:,ik,ikp] = arr
        output = karlsson.convert_to_solar(yII.elemnames,output,verbose=False)
        if savegrid: np.save(chemgrid_filename,output)
        print "Time to copy/convert[/save] grid: %f" % (time.time()-start)
    else:
        #already have a list of histograms; weight and combine them
        outdict = {}
        start = time.time()
        for irow,erow in enumerate(elemnames):
            for icol,ecol in enumerate(elemnames):
                if icol > irow: continue
                for i,thisdict in enumerate(chemlist):
                    if savetmpfiles:
                        f = open(_get_tmpdictname(i,tmplabel),'r')
                        thisdict = pickle.load(f)
                        f.close()
                    ik,ikp = divmod(i,kpmax+1); k = ik+1; kp = ikp
                    if kp > k: continue
                    this_c = ckkp[ik,ikp]
                    if irow==icol:
                        h,x = thisdict[(irow,icol)]
                        if i==0:
                            outdict[(irow,icol)] = h*this_c,x
                        else:
                            th,tx = thisdict[(irow,icol)]
                            oh,ox = outdict[(irow,icol)]
                            assert (tx==ox).all()
                            outdict[(irow,icol)] = (oh + th*this_c,ox)
                    else:
                        h,x,y = thisdict[(irow,icol)]
                        if i==0:
                            outdict[(irow,icol)] = h*this_c,x,y
                        else:
                            th,tx,ty = thisdict[(irow,icol)]
                            oh,ox,oy = outdict[(irow,icol)]
                            assert (tx==ox).all() and (ty==oy).all()
                            outdict[(irow,icol)] = (oh + th*this_c,ox,oy)
                    if savetmpfiles:
                        subprocess.call(['rm '+_get_tmpdictname(i,tmplabel)],shell=True)
        print "Time to weight/combine hists: %f" % (time.time()-start)
    return outdict

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
    parser.add_option('-w','--weight',action='store_true',dest='weight',default=False)
    parser.add_option('--savegrid',action='store_true',dest='savegrid',default=False)
    parser.add_option('-s','--savetmpfiles',action='store_true',dest='savetmpfiles',default=False)
    options,args = parser.parse_args()
    if options.weight: print "-w: Calculating histograms on the fly and weighting afterwards"

    envname = args[0]; sfrname = args[1]
    fb = .1551; Mhalo = karlsson.envparams(envname)[0]
    Mmax = fb * Mhalo

    salimf  = relative_imf(np.array([13,15,18,20,25,30,40.]),2.35,norm=True) #salpeter imf
    popIIy  = yields.ryN06(salimf)
    elemnames = popIIy.elemnames
    
    popIIIy_arr = [yields.ryI05N06(p) for p in [0,.1,.3,.5,.7,.9,1.]]
    #popIIIy_arr = [yields.ryI05N06(p) for p in [.3,.5,.7,.9,1.]]
    #popIIIy_arr = [yields.ryI05N06(p) for p in [1.]]
    
    binwidth = .1
    binlist = [np.arange(-9,1+binwidth,binwidth) for i in range(6)]

    for popIIIy in popIIIy_arr:
        postfix = 'p%0.2f' % (popIIIy.p)
        kmax,kpmax,ckkp = util.load_ckkp(envname,sfrname)
        chemgrid = run_compute_chemgrid(envname,sfrname,popIIy,popIIIy,Mmax,
                                        savegrid=options.savegrid,Nstars=10**6,
                                        numprocs=options.numprocs,postfix=postfix,
                                        weighthist=options.weight,binlist=binlist,
                                        savetmpfiles=options.savetmpfiles)
        if not options.weight:
            start = time.time()
            outdict = karlsson.hist_chemgridkkp(ckkp,chemgrid,binlist,elemnames)
            print "Histogramming done: %f" % (time.time()-start)
        else:
            outdict = chemgrid
        outdict['yII'] = popIIy
        outdict['yIII'] = popIIIy
        outdict['binlist'] = binlist

        histname = karlsson.get_chemgridhist_filename(envname,sfrname,postfix=postfix)
        f = open(histname,'w')
        pickle.dump(outdict,f)
        f.close()
