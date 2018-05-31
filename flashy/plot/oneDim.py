import flashy.datahaul.hdf5yt as reader
import flashy.utils as ut
from .globals import *
import matplotlib.pyplot as plt
from matplotlib import gridspec
from ..datahaul.plainText import dataMatrix
from ..nuclear import sortNuclides, elemSplit

def fileProfile(fname, species=[], thresh=1e-6, sprange=[0.0, 0.0], 
                filetag='prof', show=False, byM=True, 
                geom='cartesian', direction=[]):
    """plots a checkpoint file ray through the domain.
    
    Args:
        fname (str): filename of checkpoint.
        species (list of str): list of species names to plot, defaults to all.
        thresh (float): ymin for species fraction plot.
        sprange (list of float): if set, change the range of the abundance plot.
        filetag (str): change prefix of png output files.
        show (bool): return figure instead of saving it to file.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
        geom (str): geometry (['cartesian'], 'spherical').
        direction (list of float): angles of lineout. (polar, azimuth)
            also sets dimensionality (empty is 1D).
    
    """
    if 'chk' not in fname:
        print("Not a checkpoint, skipping abundance plots")
        plotsp = False
    else:
        plotsp = True
    ad, allsp = reader.getLineout(fname, geom=geom, direction=direction)
    time, _, _, _, paths = reader.getMeta(fname)
    fig = plt.figure(figsize=(13, 8))

    if plotsp:
        layout = (3, 2)
        ax4 = plt.subplot2grid(layout, (0, 1), aspect="auto", adjustable='box-forced', rowspan=2)
        ax4.set_prop_cycle(cc)
    else:
        layout = (3,1)
    ax1 = plt.subplot2grid(layout, (0, 0), aspect='auto')
    ax2 = plt.subplot2grid(layout, (1, 0), aspect="auto", sharex=ax1, adjustable='box-forced')
    ax3 = plt.subplot2grid(layout, (2, 0), aspect="auto", sharex=ax1, adjustable='box-forced')

    ax1.loglog(ad[0], ad[1], color='black')
    ax1.set_ylabel('$\\frac{g}{cm^3}$', rotation=0, labelpad=15)
    ax1.annotate("{:.5f} s".format(time),
                 xy=(0.0, 0.0), xytext=(0.05, 0.15), size=12,
                 textcoords='axes fraction', xycoords='axes fraction')
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))

    ax2.loglog(ad[0], ad[2], color='red')
    ax2.set_ylabel('$K$', rotation=0, labelpad=15)
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    ax3.loglog(ad[0], ad[3], color='blue')
    ax3.set_ylabel('$\\frac{dyn}{cm^2}$', rotation=0, labelpad=15)
    ax3.set_xlabel('Radius ($cm$)')
    ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    if plotsp:
        if not species:
            species = allsp
        spoffset = len(ad)-len(allsp)
        allsp = sortNuclides(allsp)
        if byM:
            xs = ut.byMass(ad[0], ad[1])
            ax4.set_xlabel('Mass ($M_{\odot}$)')
        else:
            xs = ad[0]
            ax4.set_xlabel('Radius ($cm$)')
        for s in range(len(allsp)):
            if allsp[s] not in species:
                continue
            tag = '$^{{{}}}{}$'.format(*elemSplit(allsp[s], invert=True))
            ax4.semilogy(xs, ad[s+spoffset], label=tag, alpha=0.9)
        lgd = ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
                         columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                         numpoints=3, handletextpad=0.0)
        if sum(sprange)!=0.0:
            ax4.set_xlim(sprange)
        ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
        ax4.set_ylim(thresh, 2.0)
        ax4.set_ylabel('$X_{i}$', rotation=0, labelpad=10)
    if plotsp and len(species)<30:
        plt.tight_layout(pad=1.0, h_pad=0.05, w_pad=0.5, rect=(0, 0, 1.0,1.0))
    else:
        plt.tight_layout(pad=1.0, h_pad=0.05, w_pad=0.5, rect=(0, 0, 0.65,1.0))
    plt.subplots_adjust(hspace=0.001)
    
    if show:
        return
    else:
        num = paths[1][-5:]
        otpf, _ = os.path.split(paths[0])
        tag = filetag
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp,  bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)
        print("Wrote: {}".format(savp))

        
def flashProfile(fname, species=[], thresh=1e-6, xrange=[0.0, 0.0], 
                 filetag='prof', show=False, byM=True, 
                 geom='cartesian', direction=[]):
    
    ad, allsp = reader.getLineout(fname, geom=geom, direction=direction, 
                                  srcnames=False)
    time, _, _, _, paths = reader.getMeta(fname)
    keys = ['radius', 'density', 'temperature', 'pressure']
    keys +=allsp
    prof = dataMatrix([keys, ad.transpose()])
    fig = plotDMat(prof, byM=byM, species=species, 
                   thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.84, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='g'))
    if show:
        return fig
    else:
        num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(paths[0], filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))
    
        
def plainTprofile(fname, species=[], thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots main properties of a plain text file.
    
    Args:
        fname (str): filename of checkpoint.
        species (list of str): list of species names to plot, defaults to all.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """
    prof = dataMatrix(fname)
    return plotDMat(prof, species=species, thresh=thresh, xrange=xrange, byM=byM)
    

def plotDMat(prof, species=[], thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots main properties of a profile object.
    
    Args:
        prof (dataMatrix): dataMatrix obj.
        species (list of str): list of species names to plot, defaults to all.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """  
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    # count species to format plot
    if species:
        keys = sortNuclides(species) 
    else:
        keys = sortNuclides(prof.species)
        
    if len(keys)>30:
        ncol = 6
        pad = 0.0
    else:
        ncol = 2
        pad = -1.0

    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False  
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (len(plotp)+2, 2)
    ax1 = plt.subplot2grid(layout, (0, 0))
    for s in range(len(keys)):
        tag = '$^{{{}}}{}$'.format(*elemSplit(keys[s], invert=True))
        ax1.semilogy(xs, getattr(prof, keys[s]), label=tag, alpha=0.9)
    lgd = ax1.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.05), 
                     columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                     numpoints=3, handletextpad=0.0)
    ax1.axhline(1e0, linewidth=1, linestyle=':', color='black')
    ax1.set_ylim(thresh, 2.0)
    ax1.set_ylabel('$X_{i}$', rotation=0, labelpad=10)
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax1.tick_params(labelbottom=False) 
    if log:
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        
    ax2 = plt.subplot2grid(layout, (1, 0), sharex=ax1)
    ax2.semilogy(xs, prof.density)
    ax2.set_ylabel('$\\frac{g}{cm^3}$', rotation=0, labelpad=15)
#     ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.tick_params(labelbottom=False) 
    
    for i, p in enumerate(plotp):
#         print(i, p, plotp)
        ax3 = plt.subplot2grid(layout, (i+2, 0), sharex=ax1)
        for j in range(i+1):
            ax3.plot([], [])
        ax3.semilogy(xs, getattr(prof, p))#, color='k')
        ax3.set_ylabel(p.capitalize(), rotation=90, labelpad=15)
        ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
        ax3.tick_params(labelbottom=False)
#     ax2.set_xlabel('Radius ($cm$)')
    ax3.set_xlabel(xlab)
    ax3.tick_params(labelbottom=True)
    if sum(xrange)!=0.0:
        ax1.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=pad)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotSpecies(ax, dmatr, byMass=True, thresh=1e-4):
    """draws species from a profile object.
    
    Args:
        ax(mpl.axes): axes instance to draw on.
        dmatr(dataMatrix): profile object to extract data.
        byMass(bool): plot by Mass or radius (toggle).
        thresh(float): lower bound for plot.
    
    """
    if byMass:
        absc = 'mass'
    else:
        absc = 'radius'
    for s in dmatr.species:
        line = simplePlot(ax, dmatr, absc, s, log=True)
        line[0].set_label(s)
    l = ax.set_ylabel('$X_i$')
    l = ax.set_ylim([thresh, 2.0])
    

def simplePlot(ax, dmatr, absc, attr, log=True, **kwargs):
    """draws a pair of properties from a profile object.
    
    Args:
        ax(mpl.axes): axes instance to draw on.
        dmatr(dataMatrix): profile object to extract data.
        absc(str): abscissa for the plot. 
        attr(str): attribute to plot.
        
    """
    if log:
        line = plt.loglog(getattr(dmatr, absc), getattr(dmatr, attr), **kwargs)
    else:
        line = plt.plot(getattr(dmatr, absc), getattr(dmatr, attr), **kwargs)
    l = plt.xlabel(absc)
    l = plt.ylabel(attr)
    return line
