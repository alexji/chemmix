import numpy as np
import pylab as plt
import matplotlib.cm as cm
import util

def plot_ckkp(ckkp,envname=None,sfrname=None,
              ax=None,maxx=None,maxy=None,
              writelabelx=True,writelabely=True,
              title=None):
    if envname!=None and sfrname != None:
        kmax,ckkp = util.load_ckkp(envname,sfrname,full_grid=True)
    kmax = ckkp.shape[0] #kmax * (kmax+1)

    logminbin = -8; logmaxbin = 0
    if ax == None: ax=plt.gca()

    im = ax.matshow(np.log10(ckkp),vmin=logminbin,vmax=logmaxbin,
                    extent=(0,kmax,kmax,1))
    ax.xaxis.tick_bottom()
    if maxx!=None:
        ax.set_xlim(0,maxx)
    if maxy!=None:
        ax.set_ylim(maxy,1)
    if (title==None): 
        ax.set_title(util.default_filename(envname,sfrname))
    else: ax.set_title(title)
    if writelabelx: ax.set_xlabel(r'$k^\prime$ (Pop III SN)')
    else: ax.set_xticklabels([])
    if writelabely: ax.set_ylabel(r'$k$ (total SN)')
    else: ax.set_yticklabels([])
    return im

if __name__=="__main__":
    fig,axarr = plt.subplots(2,4)#,sharex=True,sharey=True)

    envname = 'atomiccoolinghalo'
    sfrnames = ['fidTS'+str(tcut) for tcut in (300,400,500,600,700,800,900,1000)]
    for iplot,sfrname in enumerate(sfrnames):
        irow,icol = divmod(iplot,axarr.shape[1])
        im = plot_ckkp(0,envname=envname,sfrname=sfrname,
                       ax=axarr[irow,icol],maxx=60,maxy=110,title=sfrname,
                       writelabelx=(irow==(axarr.shape[0]-1)),writelabely=(icol==0))
    cax = fig.add_axes([.9,.1,.03,.8])
    fig.colorbar(im,cax=cax)
    fig.subplots_adjust(wspace=0)
    plt.savefig('PLOTS/fid_ckkp.png',bbox_inches='tight')
    plt.show()
