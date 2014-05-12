import numpy as np
import pylab as plt

if __name__=="__main__":
    tmin = 0.0; tmax = 1000.0; dt = 0.1
    #FILENAME = "TESTb_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    #FILENAME = "TESTc_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
    #Mbins = np.linspace(0,3,301)*10**6
    FILENAME = "Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
    Mbins = np.linspace(0,5,501)*10**6
    data = np.load(FILENAME+'.npz')
    tarr = data['tarr']
    mgrid = data['Mmixgrid']
    gtrzero = np.where(mgrid > 0)
    mlist = mgrid[gtrzero].reshape((-1,1))
    
    plt.figure()
    plt.imshow(np.log10(mgrid.transpose()),
               extent=[min(tarr),max(tarr),min(tarr),max(tarr)],
               origin='lower')
    plt.colorbar()
    plt.savefig(FILENAME+'_ttaugrid.png')

    plt.figure()
    plt.hist(mlist,bins=Mbins)
    plt.savefig(FILENAME+'_Mhist.png')
    
    plt.show()
