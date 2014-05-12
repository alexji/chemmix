import numpy as np
import pylab as plt
from karlsson import *


if __name__=="__main__":
    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    asplund09solar = np.array([8.43,8.69,7.60,7.51,6.34,7.50]) - 12 #log10(Z/H) solar
    tmin = 0.0; tmax = 1000.0; dt = 0.03
    Nstars = 10**5
    print "tmax,dt:",tmax,dt
    print "Nstars:",Nstars
    fileprefix="Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)

    mufn = get_mufn_K05A()
    psi_K05 = uSN_K05A
    numyields = 6
    bins = np.arange(-6.5,0.05,.05)
    parr = [0.00,.001,.01,.25]
    plt.figure()
    maxkmax = 6
    plt.subplots_adjust(hspace=0.0,wspace=0.0)
    for kmax in range(1,maxkmax+1):
        allchemarr = [np.load(fileprefix+'_chemgrid_N06_N'+str(Nstars)+'.npy'),
                      np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.001_N'+str(Nstars)+'.npy'),
                      np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.01_N'+str(Nstars)+'.npy'),
                      np.load(fileprefix+'_chemgrid_N06_Cenhanced_p0.25_N'+str(Nstars)+'.npy')]
        for i,chemarr in enumerate(allchemarr):
            allchemarr[i] = chemarr[:,:,0:kmax]

        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            ck = calc_ck(kmax,wISM_K05,mufn,psi_K05)
        for k in range(kmax):
            for i in range(len(allchemarr)):
                allchemarr[i][:,:,k] *= ck[k]
        for i in range(len(allchemarr)):
            allchemarr[i] = np.log10(np.sum(allchemarr[i],2))
            for z in range(numyields):
                allchemarr[i][:,z] -= asplund09solar[z]
                #print np.min(chemarr[:,z]), np.max(chemarr[:,z])
        #for chemarr in allchemarr:
        #    for z in range(numyields):
        #        print np.min(chemarr[:,z]), np.max(chemarr[:,z])

        #plt.title('[C/H]')
        for i,chemarr in enumerate(allchemarr):
            plt.subplot(maxkmax,len(allchemarr),1+((kmax-1)*4)+i)
            plt.hist(chemarr[:,0],bins=bins,normed=True)
            plt.yticks([])
            if kmax==1:
                plt.title('p=%4.3f' % parr[i])
            if kmax==maxkmax:
                plt.xticks([-6,-5,-4,-3,-2,-1,0])
                plt.xlabel('[C/H]')
            else:
                plt.xticks([])
            plt.xlim([-7,0])
            plt.ylim([0,2.5])
            if i==0:
                plt.ylabel(r'$k_{max}=$'+str(kmax))
    plt.savefig('carbdistr.png',bbox_inches='tight')
    #plt.show()
