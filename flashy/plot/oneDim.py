import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, StrMethodFormatter
import flashy.datahaul.hdf5yt as reader
import flashy.utils as ut
from flashy.plot._globals import *
from flashy.datahaul.plainText import getProfData

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
    ad, species = reader.getLineout(fname, geom=geom, direction=direction)
    time, _, _, _, paths = reader.getMeta(fname)
    fig = plt.figure(figsize=(12, 8))

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
        spoffset = len(ad)-len(species)-1
        species = sortNuclides(species)
        styleIter = colIter()
        # don't skip any plot to ensure colors stick to species, and legend doesn't 
        # shapeshift.
        if byM:
            xs = ut.byMass(ad[0], ad[1])
            ax4.set_xlabel('Mass ($M_{\odot}$)')
        else:
            xs = ad[0]
            ax4.set_xlabel('Radius ($cm$)')
        for s in range(len(species)):
            tag = '$^{{{}}}{}$'.format(*elemSplit(species[s], invert=True))
            c, ls = styleIter.next()
            ax4.semilogy(xs, ad[s+spoffset], label=tag, color=c, linestyle=ls, alpha=0.7)
        lgd = ax4.legend(ncol=5, loc='upper left', bbox_to_anchor=(1.0, 1.0), 
          columnspacing=0.5, labelspacing=0.5, markerfirst=False, 
          numpoints=4)
        if sum(sprange)!=0.0:
            ax4.set_xlim(sprange)
        ax4.axhline(1e0, linewidth=1, linestyle=':', color='black')
        ax4.set_ylim(thresh, 2.0)
        ax4.set_ylabel('$X_{i}$', rotation=0, labelpad=15)
    plt.tight_layout(pad=1.0, h_pad=0.0, w_pad=0.5, rect=(0,0,0.67,1))
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


def plotter(ax, filename, key, xkey='Radius', label='', log=1, **kwargs):
    """plots 'key' values from a data dictionary"""
    dataDict = getProfData(filename)
    ax.set_ylabel(key)
    ax.set_xlabel(xkey)
    if log:
        ax.loglog(dataDict[xkey], dataDict[key], label=label, **kwargs)
    else:
        ax.plot(dataDict[xkey], dataDict[key], label=label, **kwargs)


def percentDiff(ax, file1, file2, diffkey, xkey='Radius', label='', **kw):
    """plots the percentage difference for diffkey for two filenames."""
    d1 = getProfData(file1)
    d2 = getProfData(file2)
    jake = np.interp(d1[xkey], d2[xkey], d2[diffkey])
    norm = np.max(jake)
    ax.plot(d1[xkey], 100*(d1[diffkey]-jake)/norm, label=label, **kw)