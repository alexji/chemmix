import numpy as np
import pylab as plt
import time

from karlsson import *

if __name__=="__main__":
    label = "10000_dt0.1"
    def get_filename(k):
        return "fM"+str(k)+"_K05_"+label

    Mbins = np.linspace(0,3,301)*10**6
    Mplot = (Mbins[:-1]+Mbins[1:])/2.0

    mufn = get_mufn_K05A()
    aLMS = 0.835
    aSN  = 0.002
    psi_K05 = lambda t: uSN_K05A(t)/aLMS
    
# Retrieve DMList
    NUMPROCS = 20
    filename = "Mmixgrid_"+label
    saved_files = np.load(filename+".npz")
    tarr = saved_files['tarr']; tauarr = saved_files['tauarr']
    Mmixgrid = saved_files['Mmixgrid']
    start = time.time()
    DMlist = get_DMlist(Mbins,tarr,tauarr,Mmixgrid,numprocs=NUMPROCS)
    print "get_DMlist time:",time.time()-start

# Calculate and save fMk
    plt.figure()
    for k in range(1,20):
        start = time.time()
        fMk = calc_fMk(k,Mbins,DMlist,Vmix_K05,wISM_K05,mufn,psi_K05)
        print "calc_fM",k,time.time()-start
        np.save(get_filename(k),fMk)
        plt.plot(Mplot,fMk,label=str(k))
    plt.savefig('run_k05.png',bbox_inches='tight')

## Test calc_ck
#    ck = calc_ck(19,wISM_K05,mufn,uSN_K05A)
#    Mmixweighted = np.zeros(len(Mplot))
#    for ck,k in zip(ck,range(1,20)):
#        fMk = np.load(get_filename(k)+'.npy')
#        Mmixweighted += fMk*ck
#    plt.plot(Mplot,Mmixweighted)
#    plt.show()

## Test draw_from_distr
#    fM1 = np.load(get_filename(1)+'.npy')
#    Mrand = draw_from_distr(10000,Mplot,fM1)
#    plt.hist(Mrand,bins=Mbins)
#    plt.show()

    
