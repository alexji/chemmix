import numpy as np
import pylab as plt
import h5py
from optparse import OptionParser
import scipy.stats

import karlsson
from tophat import TopHat

def get_subplot_num(row,col,numrows,numcols):
    assert row < numrows; assert col < numcols
    return 1+row*numcols+col

if __name__=="__main__":
    numrows = 5; numcols = 2
    plotname = "cfrac_feh.png"
    minihalofiles = [#'CHEMGRIDS/minihalo_chemgrid_N06i.hdf5',
                     #'CHEMGRIDS/minihalo_chemgrid_N06ip0.100f100.0.hdf5',
                     'CHEMGRIDS/minihalo_chemgrid_N06_M25.hdf5',
                     'CHEMGRIDS/minihalo_chemgrid_HW10E1.2S4m0.hdf5',
                     #'CHEMGRIDS/minihalo_chemgrid_HW10E1.2S4m0_a1.35.hdf5',
                     'CHEMGRIDS/minihalo_chemgrid_HW10E1.2S4m0_flat.hdf5',
                     #'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.1.hdf5',
                     'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.5.hdf5',
                     'CHEMGRIDS/minihalo_chemgrid_mixN06HW10p0.9.hdf5']
    atomiccoolinghalofiles = [#'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06i.hdf5',
                              #'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06ip0.100f100.0.hdf5',
                              'CHEMGRIDS/atomiccoolinghalo_chemgrid_N06_M25.hdf5',
                              'CHEMGRIDS/atomiccoolinghalo_chemgrid_HW10E1.2S4m0.hdf5',
                              #'CHEMGRIDS/atomiccoolinghalo_chemgrid_HW10E1.2S4m0_a1.35.hdf5',
                              'CHEMGRIDS/atomiccoolinghalo_chemgrid_HW10E1.2S4m0_flat.hdf5',
                              #'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.1.hdf5',
                              'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.5.hdf5',
                              'CHEMGRIDS/atomiccoolinghalo_chemgrid_mixN06HW10p0.9.hdf5']
    yieldlabels = [#'N06',
                   #'N06 Cx100 (p=.1)',
                   'N06 M25',
                   'HW10 a2.35', 
                   #'HW10 a1.35', 
                   'HW10 flat', 
                   #'NHW (p=.1)', 
                   'NHW (p=.5)',
                   'NHW (p=.9)']

    allfiles = [minihalofiles,atomiccoolinghalofiles]

    assert numcols == len(allfiles)
    for filelist in allfiles: assert numrows == len(filelist)

    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    bins = np.arange(-6.0,-1.0+0.01,0.1)

    plt.figure(figsize=(8,10))
    # Minihalo column
    for col,filelist in enumerate(allfiles):
        if col == 0: 
            kmax = 2
            paramfn = karlsson.params_minihalo
            plt.subplot(numrows,numcols,1); 
            plt.title(r'minihalo $k_{\rm max}='+str(kmax)+'$')
        if col == 1: 
            kmax = 5
            paramfn = karlsson.params_atomiccoolinghalo
            plt.subplot(numrows,numcols,2); 
            plt.title(r'atomic cooling halo $k_{\rm max}='+str(kmax)+'$')
        for row,filename in enumerate(filelist):
            plt.subplot(numrows,numcols,get_subplot_num(row,col,numrows,numcols))
            f = h5py.File(filename,'r')
            chemarr = np.array(f['chemgrid'])
            f.close()
            chemarr,ck = karlsson.weight_chemgrid(kmax,chemarr,paramfn,
                                                  elemnames=elemnames,verbose=True)
            FeH,Cfrac,CfracErr = karlsson.compute_cfrac(chemarr[:,0],chemarr[:,5],
                                                        bins=bins)
            h,x = np.histogram(chemarr[:,5],bins=bins)
            h = h/float(np.max(h))
            plt.bar(x[:-1],h,np.diff(x),
                    edgecolor='#D3D3D3',color='#D3D3D3')
            plt.errorbar(FeH,Cfrac,yerr=CfracErr,
                         color='black',marker='.',drawstyle='steps-mid')
            plt.ylim((0,1.1)); plt.xlim((np.min(bins),np.max(bins)))
            #plt.ylabel(r'$N$, $f_{Crich}$')

            plt.plot([-4,-3],[0.75,0.3],'bo-')
            if col == 0: plt.ylabel(yieldlabels[row])
            else: plt.gca().yaxis.set_ticklabels([])
            if row == numrows-1: plt.xlabel('[Ca/H]') #plt.xlabel('[Fe/H]')
            else: plt.gca().xaxis.set_ticklabels([])
    plt.subplots_adjust(hspace=0,wspace=0)
    #plt.savefig("PLOTS/"+plotname)
    plt.show()
