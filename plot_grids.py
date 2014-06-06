import numpy as np
import pylab as plt
import h5py
from optparse import OptionParser
import scipy.stats
from matplotlib import cm, colors
from matplotlib.font_manager import FontProperties

import karlsson
from tophat import TopHat
from plot_cfrac import compute_cfrac

def get_subplot_num(row,col,numrows,numcols):
    assert row < numrows; assert col < numcols
    return 1+row*numcols+col

class ImageFollower(object):
    'update image in response to changes in clim or cmap on another image'
    def __init__(self, follower):
        self.follower = follower
    def __call__(self, leader):
        self.follower.set_cmap(leader.get_cmap())
        self.follower.set_clim(leader.get_clim())

def labelaxes(xlabel,kmax,KMAXFLAG,ax,irow,icol,numrows,numcols):
    if icol == 0: 
        if KMAXFLAG: ax.set_ylabel(r'$k_{\rm max}='+str(kmax)+'$')
        else: ax.set_ylabel(r'$k='+str(kmax)+'$')
    else: ax.yaxis.set_ticklabels([])
    if irow == numrows-1: ax.set_xlabel(xlabel)
    else: ax.xaxis.set_ticklabels([])

if __name__=="__main__":
    KMAXFLAG = False
    MDIL=5

    numrows = 10; numcols = 2
    #label = 'HW10E1.2S4m0'
    #label = 'HW10E1.2S4m0_a1.35'
    #label = 'HW10E1.2S4m0_flat'
    #label = 'N06i'
    #label = 'mixN06HW10p0.5'
    #label = 'N06ip0.100f100.0'
    #label = 'N06ip0.100MC0.800000'
    #label = 'I05T07_p0.5'
    #label = 'I05T07_p0.5Mmax'
    #label = 'I05T07_p0.1'
    #label = 'I05T07_p0.1Mmax'
    #label = 'I05T07_p0.9'
    #label = 'I05T07_p0.9Mmax'

    #label = 'I05N06_p0.1'
    #label = 'I05N06_p0.5'
    #label = 'I05N06_p0.9'
    #label = 'I05N06_p0.1Mmax'
    #label = 'I05N06_p0.5Mmax'
    #label = 'I05N06_p0.9Mmax'
    #label = 'I05N06_p0.95Mmax'
    label = 'I05N06_p0.99Mmax'
    print "Using",label
    ## TODO
    #if MDIL==5:
    #    file1 = 'CHEMGRIDS/minihalo_chemgrid_'+label+'.hdf5'
    #    file2 = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_'+label+'.hdf5'
    #    paramfn1 = karlsson.params_minihalo
    #    paramfn2 = karlsson.params_atomiccoolinghalo
    #    column1 = 'minihalo'
    #    column2 = 'atomiccoolinghalo'
    #else:
    #    file1 = 'CHEMGRIDS/minihalo'+str(MDIL)+'_chemgrid_'+label+'.hdf5'
    #    file2 = 'CHEMGRIDS/atomiccoolinghalo'+str(MDIL)+'_chemgrid_'+label+'.hdf5'
    #    paramfn1 = karlsson.params_minihalo
    #    paramfn2 = karlsson.params_atomiccoolinghalo
    #    column1 = 'minihalo'
    #    column2 = 'atomiccoolinghalo'
    #    label += '_Mdil'+str(MDIL)
    #file1 = 'CHEMGRIDS/atomiccoolinghaloearly_chemgrid_'+label+'.hdf5'
    ##file2 = 'CHEMGRIDS/atomiccoolinghalo_chemgrid_'+label+'.hdf5'
    #file2 = 'CHEMGRIDS/atomiccoolinghalolate_chemgrid_'+label+'.hdf5'
    #paramfn1 = karlsson.params_atomiccoolinghaloearly
    ##paramfn2 = karlsson.params_atomiccoolinghalo
    #paramfn2 = karlsson.params_atomiccoolinghalolate
    #column1 = 'ACH early'
    ##column2 = 'ACH'
    #column2 = 'ACH late'
    #label += '_achearly-late'
    ##klist = np.arange(1,11)
    #klist = np.arange(6,16)
    
    file1 = 'CHEMGRIDS/atomiccoolinghalolate_lowvturb_chemgrid_'+label+'.hdf5'
    file2 = 'CHEMGRIDS/atomiccoolinghalo_lowmass_chemgrid_'+label+'.hdf5'
    paramfn1 = karlsson.params_atomiccoolinghalolate_lowvturb
    paramfn2 = karlsson.params_atomiccoolinghalo_lowmass
    column1 = 'ACH late lowvturb (1/10)'
    column2 = 'ACH lowmass (7.4)'
    label += '_achlowDt'
    klist = np.arange(6,16)

    elemnames = ['C', 'O', 'Mg', 'Si', 'Ca', 'Fe']
    #FeHbins = np.arange(-8.0,-1.0+0.01,0.25)
    #CFebins = np.arange(-1.0,5.0+0.01,0.25)
    #CHbins  = np.arange(-8.0,-1.0+0.01,0.25)
    FeHbins = np.arange(-9.5,0.5+0.01,0.25)
    CFebins = np.arange(-1.0,5.0+0.01,0.25)
    CHbins  = np.arange(-6.0,2.5+0.01,0.25)

    cmap = cm.Greys; images = []
    vmin = 1e40; vmax = -1e40

    fig1 = plt.figure(figsize=(8,10))
    fig2 = plt.figure(figsize=(8,10))
    fig3 = plt.figure(figsize=(8,10))
    fig1.text(0.5, 0.95, label,horizontalalignment='center',fontproperties=FontProperties(size=16))
    fig2.text(0.5, 0.95, label,horizontalalignment='center',fontproperties=FontProperties(size=16))
    fig3.text(0.5, 0.95, label,horizontalalignment='center',fontproperties=FontProperties(size=16))
    fig1.subplots_adjust(hspace=0,wspace=0)
    fig2.subplots_adjust(hspace=0,wspace=0)
    fig3.subplots_adjust(hspace=0,wspace=0)

    for icol in range(numcols):
        if icol==0: 
            print file1
            f = h5py.File(file1,'r')
            paramfn = paramfn1
        if icol==1: 
            print file2
            f = h5py.File(file2,'r')
            paramfn = paramfn2
        if icol==2: 
            print file3
            f = h5py.File(file3,'r')
            paramfn = paramfn3
        chemarr = np.array(f['chemgrid'])
        f.close()
        for irow in range(numrows):
            kmax = klist[irow]
            if KMAXFLAG:
                chemarrk,ck = karlsson.weight_chemgrid(kmax,chemarr,paramfn,elemnames=elemnames,verbose=False)
            else: #don't weight the kmax with each other
                chemarrk = karlsson.convert_to_solar(elemnames,chemarr[:,:,kmax-1],verbose=False) 
            FeH,Cfrac,CfracErr = compute_cfrac(chemarrk[:,0],chemarrk[:,5],
                                               bins=FeHbins)
            ax1 = fig1.add_subplot(numrows,numcols,get_subplot_num(irow,icol,numrows,numcols))
            ax2 = fig2.add_subplot(numrows,numcols,get_subplot_num(irow,icol,numrows,numcols))
            ax3 = fig3.add_subplot(numrows,numcols,get_subplot_num(irow,icol,numrows,numcols))
            ## Fig 1
            h,x = np.histogram(chemarrk[:,5],bins=FeHbins)
            h = h/float(np.max(h))
            ax1.bar(x[:-1],h,np.diff(x),
                    edgecolor='#D3D3D3',color='#D3D3D3')
            ax1.errorbar(FeH,Cfrac,yerr=CfracErr,
                         color='black',marker='.',drawstyle='steps-mid')
            ax1.set_ylim((0,1.1)); ax1.set_xlim((np.min(FeHbins),np.max(FeHbins)))
            ax1.plot([-4,-3],[0.75,0.3],'bo-')
            labelaxes('[Fe/H]',kmax,KMAXFLAG,ax1,irow,icol,numrows,numcols)
            ## Fig 2
            xmin = np.min(FeHbins); xmax = np.max(FeHbins)
            ymin = np.min(CFebins); ymax = np.max(CFebins)
            H,xedges,yedges = np.histogram2d(chemarrk[:,0]-chemarrk[:,5],chemarrk[:,5],bins=(CFebins,FeHbins))
            vmin = min(vmin, np.amin(H)); vmax = max(vmax, np.amax(H))
            images.append(ax2.imshow(H,extent=[xmin,xmax,ymin,ymax],
                       interpolation='nearest',origin='lower',
                       aspect='auto',cmap=cmap))
            #X, Y = 0.5*(xedges[1:]+xedges[:-1]), 0.5*(yedges[1:]+yedges[:-1])
            #ax2.contour(X,Y,H.T, origin='lower') #doesn't work I think
            ax2.plot([xmin,xmax],[.75,.75],'r:')
            ax2.set_xlim((xmin,xmax)); ax2.set_ylim((ymin,ymax))
            labelaxes('[Fe/H]',kmax,KMAXFLAG,ax2,irow,icol,numrows,numcols)
            ## Fig 3
            h,x = np.histogram(chemarrk[:,0],bins=CHbins)
            h = h/float(np.max(h))
            ax3.bar(x[:-1],h,np.diff(x),
                    edgecolor='#D3D3D3',color='#D3D3D3')
            ax3.set_ylim((0,1.1)); ax3.set_xlim((np.min(CHbins),np.max(CHbins)))
            labelaxes('[C/H]',kmax,KMAXFLAG,ax3,irow,icol,numrows,numcols)
    fig1.axes[0].set_title(column1)
    fig2.axes[0].set_title(column1)
    fig3.axes[0].set_title(column1)
    fig1.axes[numrows].set_title(column2)
    fig2.axes[numrows].set_title(column2)
    fig3.axes[numrows].set_title(column2)
    #fig1.axes[2*numrows].set_title(column3)
    #fig2.axes[2*numrows].set_title(column3)
    #fig3.axes[2*numrows].set_title(column3)
    
    norm = colors.Normalize(vmin=vmin,vmax=vmax)
    for i,im in enumerate(images):
        im.set_norm(norm)
        if i>0:
            images[0].callbacksSM.connect('changed', ImageFollower(im))
    cax = fig2.add_axes([0.2, 0.02, 0.6, 0.01])
    fig2.colorbar(images[0], cax, orientation='horizontal')

    if KMAXFLAG:
        fig1.savefig('PLOTS/cfrac_feh_kmax_'+label+'.png')
        fig2.savefig('PLOTS/cfe_feh_kmax_'+label+'.png')
        fig3.savefig('PLOTS/ch_kmax_'+label+'.png')
    else:
        fig1.savefig('PLOTS/cfrac_feh_k_'+label+'.png')
        fig2.savefig('PLOTS/cfe_feh_k_'+label+'.png')
        fig3.savefig('PLOTS/ch_k_'+label+'.png')
    plt.show()
