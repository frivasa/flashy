import flashy.datahaul.hdf5yt as reader
import flashy.utils as ut
from flashy.plot._globals import *
import matplotlib.pyplot as plt
from flashy.datahaul.plainText import dataMatrix
from flashy.nuclear import sortNuclides, elemSplit

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
        print "Not a checkpoint, skipping abundance plots"
        plotsp = False
    else:
        plotsp = True
    ad, allsp = reader.getLineout(fname, geom=geom, direction=direction)
    time, _, _, _, paths = reader.getMeta(fname)
    fig = plt.figure(figsize=(13, 8))

    if plotsp:
        layout = (3, 2)
        ax4 = plt.subplot2grid(layout, (0, 1), aspect="auto", adjustable='box-forced', rowspan=2)
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
    #ax3.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))

    ax2.loglog(ad[0], ad[2], color='red')
    ax2.set_ylabel('$K$', rotation=0, labelpad=15)
    #ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    ax3.loglog(ad[0], ad[3], color='blue')
    ax3.set_ylabel('$\\frac{dyn}{cm^2}$', rotation=0, labelpad=15)
    ax3.set_xlabel('Radius ($cm$)')
    #ax3.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
    
    if plotsp:
        if not species:
            species = allsp
        spoffset = len(ad)-len(allsp)
        allsp = sortNuclides(allsp)
        styleIter = colIter()
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
            c, ls = styleIter.next()
            ax4.semilogy(xs, ad[s+spoffset], label=tag, color=c, linestyle=ls, alpha=0.9)
        lgd = ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
          columnspacing=0.0, labelspacing=0.0, markerfirst=True, 
          numpoints=2)
        if sum(sprange)!=0.0:
            ax4.set_xlim(sprange)
        ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
        ax4.set_ylim(thresh, 2.0)
        ax4.set_ylabel('$X_{i}$', rotation=0, labelpad=10)
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
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp,  bbox_extra_artists=(lgd,), bbox_inches='tight')
        plt.close(fig)
        print "Wrote: {}".format(savp)


def plainTprofile(fname, species=[], thresh=1e-6, sprange=[0.0, 0.0], 
                  filetag='prof', byM=True):
    """plots main properties of a plain text file.
    
    Args:
        fname (str): filename of checkpoint.
        species (list of str): list of species names to plot, defaults to all.
        thresh (float): ymin for species fraction plot.
        sprange (list of float): if set, change the range of the abundance plot.
        filetag (str): change prefix of png output files.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """
    prof = dataMatrix(fname)
    fig = plt.figure(figsize=(12, 9))
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    layout = (len(plotp)+1, 2)
    
    ax1 = plt.subplot2grid(layout, (0, 0), aspect='auto')
    ax1.loglog(prof.radius, prof.density, color='black')
    ax1.set_ylabel('$\\frac{g}{cm^3}$', rotation=0, labelpad=15)
    ax1.set_title('Total M: {} $M_\odot$'.format(prof.meta['mass']))
    #ax3.yaxis.set_major_formatter(StrMethodFormatter('{x:.2e}'))
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))
    colorIt = colIter2()
    for i, p in enumerate(plotp):
        print i, p, plotp
        ax2 = plt.subplot2grid(layout, (i+1, 0), aspect="auto", sharex=ax1, adjustable='box-forced')
        ax2.loglog(prof.radius, getattr(prof, p), color=colorIt.next())
        ax2.set_ylabel(p.capitalize(), rotation=90, labelpad=15)
        ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.set_xlabel('Radius ($cm$)')

    ax4 = plt.subplot2grid(layout, (0, 1), aspect="auto", adjustable='box-forced', rowspan=2)
    if species:
        keys = sortNuclides(species)
    else:
        keys = sortNuclides(prof.species)
    styleIter = colIter()
    if byM:
        xs = prof.masses
        ax4.set_xlabel('Mass ($M_{\odot}$)')
    else:
        xs = prof.radius
        ax4.set_xlabel('Radius ($cm$)')
    for s in range(len(keys)):
        tag = '$^{{{}}}{}$'.format(*elemSplit(keys[s], invert=True))
        c, ls = styleIter.next()
        ax4.semilogy(xs, getattr(prof, keys[s]), label=tag, color=c, linestyle=ls, alpha=0.9)
    lgd = ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
      columnspacing=0.0, labelspacing=0.0, markerfirst=True, 
      numpoints=2)
    if sum(sprange)!=0.0:
        ax4.set_xlim(sprange)
    ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
    ax4.set_ylim(thresh, 2.0)
    ax4.set_ylabel('$X_{i}$', rotation=0, labelpad=10)
    plt.tight_layout(pad=1.0, h_pad=0.05, w_pad=0.5, rect=(0, 0, 0.8,1.0))
    plt.subplots_adjust(hspace=0.001)
    return
