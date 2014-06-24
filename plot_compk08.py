import numpy as np
import pylab as plt
import h5py

from optparse import OptionParser

import karlsson
from solveforu import loaduIIfn,loaduIIIfn,loadmuIIfn,loadmuIIIfn
from tophat import TopHat

def labelaxes(ax,xlabel,ylabel,k,irow,icol,numrows,numcols,columnnames):
    if irow == numrows - 1: ax.set_xlabel(xlabel)
    #else: ax.xaxis.set_ticklabels([])
    ax.set_ylabel(r'$'+ylabel+'='+str(k)+'$')
    #if irow==0: ax.set_title(columnnames[icol])

if __name__=="__main__":
    ADDWIIFLAG = True

    label = 'I05N06_p0.5Mmax'
    print "Using",label
    file1 = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_'+label+'.hdf5'
    file2 = 'CHEMGRIDS/k08_chemgrid_'+label+'.hdf5'
    paramfn1 = karlsson.params_atomiccoolinghalo
    paramfn2 = karlsson.params_atomiccoolinghalo
    column1 = 'ACH'; column2 = 'ACH with K08'
    columnnames = [column1,column2]
    logMdil=5
    label += '_compk08'
    numrows = 7; numcols = 2
    klist = np.arange(14)+1 #code assumes contiguous klist starting at 1
    kmax = klist[-1]
    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    FeHbins = np.arange(-10.5,0.5+0.01,0.25)

    Mhalo,zvir,vturb,lturb,nSN,trecovery = paramfn1()
    vturb *= 3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    th1 = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
    RHO1 = th1.get_rho_of_t_fn()
    VMIX1 = karlsson.get_Vmixfn_K08(RHO1,Dt=Dt,Mdil=10**logMdil)
    WISM1 = karlsson.wISM_K05
    UFN1  = lambda t: uSN if type(t)==float else uSN+np.zeros(len(t)); print uSN
    MUFN1 = karlsson.get_mufn(VMIX1,UFN1)
    ck1 = karlsson.calc_ck(kmax,WISM1,MUFN1,UFN1,tmin=.01,tmax=999.9)
    print "Finished",column1

    Mhalo,zvir,vturb,lturb,nSN,trecovery = paramfn2()
    vturb *= 3.16/3.08 * .001 #km/s to kpc/Myr
    Dt =  vturb * lturb / 3.0 #kpc^2/Myr
    uSN = nSN/(trecovery * (4*np.pi/3.) * (10*lturb)**3) #SN/Myr/kpc^3
    th2 = TopHat(Mhalo=Mhalo,zvir=zvir,fb=0.1551,mu=1.4)
    RHO2 = th2.get_rho_of_t_fn()
    VMIX2 = karlsson.get_Vmixfn_K08(RHO1,Dt=Dt,Mdil=10**logMdil)
    uIIfn  = loaduIIfn('atomiccoolinghalo'); uIIIfn = loaduIIIfn('atomiccoolinghalo')
    WISM2  = karlsson.wISM_III
    UFN2   = uIIfn
    MUFN2  = loadmuIIIfn('atomiccoolinghalo')#karlsson.get_mufn(VMIX2,uIIIfn)
    MUFNII_2= loadmuIIfn('atomiccoolinghalo')#karlsson.get_mufn(VMIX2,uIIfn)
    wISM2_0 = lambda t: karlsson.wISM_K05(0,MUFNII_2(t))
    ck2 = karlsson.calc_ck(kmax,WISM2,MUFN2,UFN2,tmin=.01,tmax=999.9)
    print "Finished",column2

    plotprefix = 'compk08_'
    if ADDWIIFLAG:
        ck1 = karlsson.calc_ck(kmax,WISM2,MUFN2,UFN2,tmin=.01,tmax=999.9,wISM2_0=wISM2_0)
        ck2 = karlsson.calc_ck(kmax,WISM1,MUFN1,UFN1,tmin=.01,tmax=999.9,wISM2_0=wISM2_0)
        plotprefix = 'compk08wII_'
    plotprefix += label

    f = plt.figure()
    ax = f.gca()
    ax.bar(klist-.5,ck1,1,edgecolor='blue',color='blue',alpha=.2,label=column1)
    ax.bar(klist-.5,ck2,1,edgecolor='red',color='red',alpha=.2,label=column2)
    ax.set_xlabel('k'); ax.set_ylabel('ck')
    ax.legend(fontsize='small',loc='best')
    plt.savefig('PLOTS/'+plotprefix+'ck.png')

    f,((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3,sharex=True)
    dummy1,tarr,dummy2 = loadmuIIfn('atomiccoolinghalo',retarrays=True)
    for k in klist:
        ax1.plot(tarr,WISM1(k,MUFN1(tarr)),label=str(k))
        ax4.plot(tarr,WISM2(k,MUFN2(tarr)))
    ax1.legend(fontsize='xx-small',loc='best')
    ax1.set_title('wISM '+column1); ax4.set_title('wISM '+column2)
    ax1.set_ylim((0,.4)); ax4.set_ylim((0,.4))

    ax2.plot(tarr,UFN1(tarr)); ax2.set_title('uII (arbitrary flat) '+column1)
    ax5.plot(tarr,UFN2(tarr)); ax5.plot(tarr,10*uIIIfn(tarr)); ax5.set_title('uII/III '+column2)
    ax2.set_yticks([]); ax5.set_ylim((0,.5))

    ax3.plot(tarr,MUFN1(tarr)); ax3.set_title('muIII '+column1)
    ax6.plot(tarr,MUFN2(tarr)); ax6.set_title('muIII '+column2)

    ax4.set_xlabel('t (Myr)'); ax5.set_xlabel('t (Myr)'); ax6.set_xlabel('t (Myr)')
    plt.savefig('PLOTS/'+plotprefix+'wISM.png')

    f,axarr = plt.subplots(7,2,sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    f1 = h5py.File(file1,'r')
    chemarr1 = np.array(f1['chemgrid'])
    f1.close()
    f2 = h5py.File(file2,'r')
    chemarr2 = np.array(f2['chemgrid'])
    f2.close()
    histlist1 = []; histlist2 = []
    for icol in range(numcols):
        for irow in range(numrows):
            ax = axarr[irow,icol]
            k = icol*numrows + irow + np.min(klist)
            chemarr1k = karlsson.convert_to_solar(elemnames,chemarr1[:,:,k-1],verbose=False)
            chemarr2k = karlsson.convert_to_solar(elemnames,chemarr2[:,:,k-1],verbose=False)
            h1,x1 = np.histogram(chemarr1k[:,5],bins=FeHbins)
            h2,x2 = np.histogram(chemarr2k[:,5],bins=FeHbins)
            histlist1.append(h1); histlist2.append(h2)

            ax.bar(x1[:-1],h1,np.diff(x1),
                   edgecolor='blue',color='blue',alpha=.2)
            ax.bar(x2[:-1],h2,np.diff(x2),
                   edgecolor='red',color='red',alpha=.2)
            labelaxes(ax,'[Fe/H]','k',k,irow,icol,numrows,numcols,columnnames)
            ax.set_ylim((0,6*10**5)); ax.set_xlim((np.min(FeHbins),np.max(FeHbins)))
    plt.savefig('PLOTS/'+plotprefix+'k.png')
    
    f,axarr = plt.subplots(7,2,sharex=True,sharey=True)
    f.subplots_adjust(hspace=0)
    histlist1 = np.array(histlist1)
    histlist2 = np.array(histlist2)

    for icol in range(numcols):
        for irow in range(numrows):
            ax = axarr[irow,icol]
            kmax = icol*numrows + irow + np.min(klist)
            thisck1 = ck1[0:kmax]; thisck1 = thisck1/np.sum(thisck1)
            thisck2 = ck2[0:kmax]; thisck2 = thisck2/np.sum(thisck2)

            thishistlist1 = np.transpose(histlist1[0:kmax,:])
            thishistlist1 = thishistlist1 * thisck1
            thishist1 = np.sum(thishistlist1,axis=1)

            thishistlist2 = np.transpose(histlist2[0:kmax,:])
            thishistlist2 = thishistlist2 * thisck2
            thishist2 = np.sum(thishistlist2,axis=1)

            ax.bar(x1[:-1],thishist1,np.diff(x1),
                   edgecolor='blue',color='blue',alpha=.2)
            ax.bar(x2[:-1],thishist2,np.diff(x2),
                   edgecolor='red',color='red',alpha=.2)
            labelaxes(ax,'[Fe/H]','k_{max}',kmax,irow,icol,numrows,numcols,columnnames)
            ax.set_ylim((0,3*10**5)); ax.set_xlim((np.min(FeHbins),np.max(FeHbins)))
    plt.savefig('PLOTS/'+plotprefix+'kmax.png')

    plt.show()
    
