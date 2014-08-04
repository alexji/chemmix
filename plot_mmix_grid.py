import numpy as np
import pylab as plt
import util
import backgroundmodel as bgm
from optparse import OptionParser

if __name__=="__main__":
    parser = OptionParser()
    options,args = parser.parse_args()
    envname=args[0]

    tarr = bgm.gettarr(envname); tmax = np.max(tarr)
    mmg = util.load_mmixgrid(envname)
    mmg = np.log10(mmg) #plot log
    print np.nanmin(mmg[np.isfinite(mmg)]),np.nanmax(mmg)
    if 'minihalo' in envname:
        levels = np.arange(1,5.5,.2)
    if 'atomiccoolinghalo' in envname:
        levels = np.arange(1,7.5,.2)

    plt.figure()
    CS = plt.contour(mmg,levels,extent=(0,tmax,0,tmax),aspect=1)
    plt.clabel(CS,levels,inline=1,fontsize=8,fmt="%1.1f")
    plt.ylabel('t'); plt.xlabel('tau')
    plt.savefig('PLOTS/mmixgrid_'+envname+'.png',bbox_inches='tight')

    #plt.figure()
    #plt.hist(mmg,bins=levels)
    #plt.savefig('PLOTS/mmixgridhist_'+envname+'.png',bbox_inches='tight')


#    tmin = 0.0; tmax = 1000.0; dt = 0.1
#    #FILENAME = "TESTb_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
#    #FILENAME = "TESTc_Mmixgrid_K05_tmax"+str(tmax)+"_dt"+str(dt)
#    #Mbins = np.linspace(0,3,301)*10**6
#    FILENAME = "Mmixgrid_K08_tmax"+str(tmax)+"_dt"+str(dt)
#    Mbins = np.linspace(0,5,501)*10**6
#    data = np.load(FILENAME+'.npz')
#    tarr = data['tarr']
#    mgrid = data['Mmixgrid']
#    gtrzero = np.where(mgrid > 0)
#    mlist = mgrid[gtrzero].reshape((-1,1))
#    
#    plt.figure()
#    plt.imshow(np.log10(mgrid.transpose()),
#               extent=[min(tarr),max(tarr),min(tarr),max(tarr)],
#               origin='lower')
#    plt.colorbar()
#    plt.savefig(FILENAME+'_ttaugrid.png')
#
#    plt.figure()
#    plt.hist(mlist,bins=Mbins)
#    plt.savefig(FILENAME+'_Mhist.png')
#    
#    plt.show()
