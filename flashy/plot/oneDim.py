from flashy.datahaul.hdf5yt import getLineout
from flashy.datahaul.hdfdirect import directMeta
import flashy.utils as ut
from .globals import *
import matplotlib.pyplot as plt
from matplotlib import gridspec
from flashy.datahaul.plainText import dataMatrix
from flashy.nuclear import sortNuclides, elemSplit
from flashy.simulation import simulation
from scipy.integrate import trapz


def PARsimProfile(fname, simfolder='', thresh=1e-6, xrange=[0.0, 0.0], 
                  filetag='prof', batch=False, byM=False, 
                  geom='cartesian', direction=[]):
    """overloaded method for IPP. 
    WARN: slowest pos will not match for non x-axis directions.
    """
    ad, allsp = getLineout(fname, geom=geom, direction=direction, srcnames=False)
    time, _, _, _, paths = directMeta(fname)
    
    sim = simulation(simfolder)
    n, t = sim.time2step(time)
    step = sim.steps[int(n)]
    xm, ym, zm = step.slowx, step.slowy, step.slowz
    rm = np.sqrt(xm**2+ym**2+zm**2)
    keys = ['radius', 'density', 'temperature', 'pressure']
    keys +=allsp
    prof = dataMatrix([keys, ad.transpose()])
    fig = plotDMat(prof, byM=byM, 
                   thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.84, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    axlist = fig.axes
    for ax in axlist[1:]:
        ax.axvline(rm, ls=':', color='red', alpha=0.7)
    if not batch:
        return fig
    else:
        num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(os.path.dirname(paths[0]), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def flashProfile(fname, thresh=1e-6, xrange=[0.0, 0.0], 
                 filetag='prof', batch=False, byM=True, 
                 geom='cartesian', direction=[]):
    """Plot bulk properties and species in a chekpoint through a ray.
    
    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction.
        xrange(float list): abscissa range.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        geom(str): geometry of the file.
        direction(float list): list of spherical angles (alt, azimuth), empty for 1D.
    
    Returns:
        (mpl figure) or (None)
    
    """
    ad, allsp = getLineout(fname, geom=geom, direction=direction, 
                                  srcnames=False)
    time, _, _, _, paths = directMeta(fname)
    keys = ['radius', 'density', 'temperature', 'pressure']
    keys +=allsp
    prof = dataMatrix([keys, ad.transpose()])
    fig = plotDMat(prof, byM=byM, 
                   thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.8, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    if not batch:
        return fig
    else:
        num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(os.path.dirname(paths[0]), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def flashSpecies(fname, thresh=1e-6, filetag='spec', batch=False, byM=True, 
                 geom='cartesian', direction=[], plotall=False):
    """Plot species and aggregated masses in a chekpoint through a ray.
    
    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction and mass yields.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        geom(str): geometry of the file.
        direction(float list): list of spherical angles (alt, azimuth), empty for 1D.
        plotall(bool): force plotting every species found.
    
    Returns:
        (mpl figure) or (None)
    
    """
    fields = ['density', 'temperature', 'pressure', 'velx']
    ad, allsp = getLineout(fname, geom=geom, direction=direction, 
                           fields=fields, srcnames=False)
    time, _, _, _, paths = directMeta(fname)
    keys = ['radius'] + fields
    keys +=allsp
    prof = dataMatrix([keys, ad.transpose()])

    fig = plt.figure(figsize=(10, 8))
    layout = (3,3)
    # plot species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(ax1, prof, byMass=byM, thresh=thresh, plotall=plotall)
    ax1.set_xlabel('Mass ($M_{\odot}$)')
    # timestamp and legend
    a = ax1.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.65, 0.1), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    lgd = ax1.legend(ncol=4, loc='upper left', bbox_to_anchor=(1.05, 1.0), 
                    columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                    numpoints=3, handletextpad=0.0)
    lgd.get_frame().set_edgecolor('k')
    
    # plot masses
    ax2 = plt.subplot2grid(layout, (1,0), colspan=2)
    for i, sp in enumerate(prof.species):
        props = next(ax2._get_lines.prop_cycler)
        if i in skip:
            continue
        a = elemSplit(sp)[1]
        spmass = trapz(getattr(prof, sp), x=prof.masses)
        ax2.scatter(a, spmass, label=sp, color=props['color'])
    ax2.set_ylabel('Total Mass ($M_\odot$)', rotation=90, labelpad=5)
    ax2.set_xlabel('Mass Number (A)')
    ax2.set_ylim([thresh, 1.5])
    ax2.set_yscale("log", nonposy='clip')
    #ax2.yaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax2.xaxis.set_minor_formatter(StrMethodFormatter(''))
    
    # plot mass vs speed
    ax3 = plt.subplot2grid(layout, (2, 0), colspan=2)
    nbins = 60
    for i, sp in enumerate(prof.species):
        props = next(ax2._get_lines.prop_cycler)
        if i in skip:
            continue
        counts, bins = np.histogram(abs(prof.velx), bins=nbins)
        joined = sorted(zip(abs(prof.velx), getattr(prof, sp)))
        speeds, masses = zip(*joined)
        start = 0
        weights = []
        for c in counts:
            massf = sum(masses[start:start+c])
            weights.append(massf/c)
            start+=c
        weights = np.array(weights)
        bins = bins[1:]
        mpln, mplbins, patches = ax3.hist(bins, bins=len(bins), weights=weights, histtype='step', log=True, label=sp)
        ax3.set_ylim([1e-6, 2])
    ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.xaxis.set_major_formatter(customFormatter(10))
    ax3.xaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.set_xlabel('Speed ($10^{10}$ cm/s)')
    fig.subplots_adjust(hspace=0.4)
    
    if not batch:
        return fig
    else:
        num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(os.path.dirname(paths[0]), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def plainTprofile(fname, thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots main properties of a plain text file.
    
    Args:
        fname (str): filename of checkpoint.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """
    prof = dataMatrix(fname)
    return plotDMat(prof, thresh=thresh, xrange=xrange, byM=byM)
    

def plotDMat(prof, thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots main properties of a profile object.
    
    Args:
        prof (dataMatrix): dataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """  
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    keys = sortNuclides(prof.species)
    ncol = 5
    labelspace = -0.15

    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False  
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (len(plotp)+2, 3)
    # Species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(ax1, prof, byMass=byM, thresh=thresh, plotall=False)
    # remove last(lowest) yticklabel to avoid overlap
    ax1.get_yticklabels()[1].set_visible(False)
    lgd = ax1.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.05), 
                     columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                     numpoints=3, handletextpad=0.0)
    ax1.axhline(1e0, linewidth=1, linestyle=':', color='black')
    ax1.set_ylim(thresh, 2.0)
    ax1.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    ax1.yaxis.set_label_coords(labelspace, 0.5)
    ax1.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax1.tick_params(labelbottom=False) 
    if log:
        ax1.set_yscale('log')
        ax1.set_xscale('log')
        
    ax2 = plt.subplot2grid(layout, (1, 0), sharex=ax1, colspan=2)
    ax2.semilogy(xs, prof.density)
    ax2.set_ylabel('Density($\\frac{g}{cm^3}$)', size=13, rotation=90, labelpad=0)
    ax2.yaxis.set_label_coords(labelspace, 0.5)
    ax2.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax2.tick_params(labelbottom=False) 
    
    for i, p in enumerate(plotp):
        ax3 = plt.subplot2grid(layout, (i+2, 0), sharex=ax1, colspan=2)
        for j in range(i+1):
            ax3.plot([], [])
        ax3.semilogy(xs, getattr(prof, p))#, color='k')
        u = ut.getUnit(p)
        spacer = labelspace if '\\' in u else labelspace - 0.02
        ax3.set_ylabel('{}({})'.format(p.capitalize(), u), rotation=90, size=13, labelpad=0)
        ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
        ax3.yaxis.set_label_coords(spacer, 0.5)
        ax3.tick_params(labelbottom=False)
    ax3.set_xlabel(xlab)
    ax3.tick_params(labelbottom=True)
    if sum(xrange)!=0.0:
        ax1.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.0)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotSpecies(ax, dmatr, byMass=True, thresh=1e-4, plotall=False):
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
    skip = []
    if plotall:
        for s in dmatr.species:
            line = simplePlot(ax, dmatr, absc, s, log=True)
            line[0].set_label(s)
    else:
        for i, s in enumerate(dmatr.species):
            props = next(ax._get_lines.prop_cycler)
            if np.max(getattr(dmatr, s))< thresh:
                skip.append(i)
                continue
            else:
                tag = '$^{{{}}}{}$'.format(*elemSplit(s, invert=True))
                line = simplePlot(ax, dmatr, absc, s, log=True, 
                                  color=props['color'], ls=props['linestyle'])
                line[0].set_label(tag)
    l = ax.set_ylabel('$X_i$')
    l = ax.set_ylim([thresh, 2.0])
    return skip
    

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
    l = plt.xlabel(absc.capitalize())
    l = plt.ylabel(attr.capitalize())
    return line
