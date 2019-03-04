from ..datahaul.hdf5yt import getLineout
from ..datahaul.hdfdirect import directMeta
import flashy.utils as ut
from .globals import *
import matplotlib.pyplot as plt
from matplotlib import gridspec
from ..datahaul.plainText import dataMatrix
from ..nuclear import sortNuclides, elemSplit, decayYield, getMus
from ..simulation import simulation
from ..post import nonRelFermi, extRelFermi, speedHisto
from scipy.integrate import trapz


def plotSpeedHistogram(fname, geom='cartesian',
                       dimension=2, ref='x', antipode=False):
    """Returns figure with speed vs mass fraction 
    histogram for an equatorial slice (this is only a probe since it 
    does not take the whole hemisphere).
    
    Args:
        fname(str): file name.
        
    Returns:
        (mpl.figure)
    
    """
    f, ax = plt.subplots()
    he, cnos, imes, iges, bins = speedHisto(fname, geom=geom, dimension=dimension, 
                                            ref=ref, antipode=antipode)
    tmass = he + cnos + imes + iges
    tmass[tmass==0.0] = 1.0
    mpln, mplbins, patches = ax.hist(bins, bins=len(bins), weights=np.nan_to_num(iges/tmass), 
                                     histtype='step', log=True, label='IGE')
    mpln, mplbins, patches = ax.hist(bins, bins=len(bins), weights=np.nan_to_num(imes/tmass), 
                                     histtype='step', log=True, label='IME')
    mpln, mplbins, patches = ax.hist(bins, bins=len(bins), weights=np.nan_to_num(cnos/tmass), 
                                     histtype='step', log=True, label='CNO')
    mpln, mplbins, patches = ax.hist(bins, bins=len(bins), weights=np.nan_to_num(he/tmass), 
                                     histtype='step', log=True, label='He')
    ax.set_ylim([1e-6,1.5])
    ax.set_xlim([-0.5e9, 7e9])
    ax.axhline(1, ls='--', alpha=0.6, c='k')
    lgd = ax.legend()
    return f


def fetchData(fname, direction, fields=[]):
    """builds a profile dataMatrix and retrieves metadata from a file.
    
    Args:
        fname(str): filename.
        direction(float list): spherical angles for lineout (polar angle and azimuth).
        fields(str list): specify field to extract (optional).
        
    Returns:
        (float): timestamp for the file.
        (dict): parameter dictionary for the checkpoint.
        (str list): filesystem paths for the file.
        (dataMatrix): block of data for plotting.

    """
    time, pars, _, _, paths = directMeta(fname)
    if len(direction)>(pars['dimensionality']-1):
        print("WARN: Direction doesn't match dimensionality: {}".format(pars['dimensionality']))
    keys = ['radius']
    if fields:
        ad, allsp = getLineout(fname, geom=pars['geometry'], fields=fields,
                               direction=direction, srcnames=False)
        keys += fields
    else:
        ad, allsp = getLineout(fname, geom=pars['geometry'],
                               direction=direction, srcnames=False)
        keys += ['density', 'temperature', 'pressure']
    keys += allsp
    return time, pars, paths, dataMatrix([keys, ad.transpose()])


# def writeFig(fig, paths, filetag):
#     """writes figure to file according to folders in path.
    
#     Args:
#         fig(mpl.figure): matplotlib object to store.
#         paths(str list): output paths.
#         filetag(str): preffix for output file.
        
#     Returns:
#         (str): destination path of the file.
#         (str): file suffix number.
    
#     """
#     num = paths[1][-5:]  # checkpoint number 'flash_hdf5_chk_0001'
#     dest = os.path.join(os.path.dirname(paths[0]), filetag)
#     name = os.path.join(dest, '{}{}.png'.format(filetag, num))
#     os.makedirs(dest, exist_ok=True)  # bless you, p3
#     plt.savefig(name, format='png')
#     plt.close(fig)
#     print("Wrote: {}".format(name))
#     return dest, num


def shockFollow(fname, simfolder='', thresh=1e-4, batch=False, byM=False, 
                direction=[], wakesize=1e8, inward=False):
    """plot shock wake for a checkpoint. 
    WARN: slowest pos will not match for non x-axis directions.
    
    Args:
        fname(str): path to checkpoint file.
        simfolder(str): simulation run metadata folder.
        thresh(float): species threshold.
        batch(bool): write file instead of returning figure.
        byM(bool): plot by mass coordinate instead of radius.
        direction(str list): polar angle and azimuth of profile ray.
        wakesize(float): extent of plot in fraction of shock position.
        inward(bool): shock to follow.
    
    Return:
        (mpl figure) or (None)
    
    """
    fig, meta = PARsimProfile(fname, simfolder=simfolder, thresh=thresh, 
                              filetag='', batch=False, byM=byM, 
                              direction=direction, meta=True)
    # readjust the plot to center the shock
    ax = plt.gca()
    outsh, inwsh, slowx, paths = meta
    if inward:
        left, right = inwsh*0.95, inwsh*(1.0 + wakesize)
        filetag = 'inward_shock'
    else:
        left, right = outsh*(1.0 - wakesize), outsh*1.05
        filetag = 'outward_shock'
    ax.set_xlim([left, right])
    #I/O
    if not batch:
        return fig
    else:
        writeFig(fig, paths, filetag)


def PARsimProfile(fname, simfolder='', thresh=1e-4, xrange=[0.0, 0.0], 
                  filetag='metaprof', batch=False, byM=False, direction=[], meta=False):
    """bloated method for IPP.
    WARN: slowest pos will not match for non x-axis directions.
    
    Args:
        fname(str): path to checkpoint file.
        simfolder(str): simulation run metadata folder.
        thresh(float): species threshold.
        xrange(float list): abscissa range for all plots.
        filetag(str): prefix for output files (if any).
        batch(bool): write file instead of returning figure.
        byM(bool): plot by mass coordinate instead of radius.
        direction(str list): polar angle and azimuth of profile ray.
    
    Return:
        (mpl figure) or (None)
    
    """
    
    time, pars, paths, prof = fetchData(fname, direction)
    sim = simulation(simfolder)
    n, t = sim.time2step(time)
    step = sim.steps[int(n)]
    xm, ym, zm = step.slowx, step.slowy, step.slowz
    rm = np.sqrt(xm**2+ym**2+zm**2)
    fig = plotDMatMerged(prof, byM=byM, 
                         thresh=thresh, xrange=xrange)
    # add extra marks from simulation obj
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.84, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    axlist = fig.axes
    for ax in axlist[1:]:
        p = rm if not byM else prof.masses[np.where(prof.radius>rm)[0][0]]
        ax.axvline(p, ls=':', color='brown', alpha=0.5, label='slowest')
    if sim.CJ:
        times, xins, cjins, xouts, cjouts = np.genfromtxt(sim.CJ[0], unpack=True)
        if byM:
            xo, xi = xouts[pars['checkpointfilenumber']], xins[pars['checkpointfilenumber']]
            no = np.where(prof.radius>xo)[0][0]
            ni = np.where(prof.radius>xi)[0][0]
            po, pi = prof.masses[no], prof.masses[ni]
        else:
            po, pi = xouts[pars['checkpointfilenumber']], xins[pars['checkpointfilenumber']]
        ax.axvline(po, ls=':', color='green', alpha=0.5, label='shock')
        ax.axvline(pi, ls=':', color='red', alpha=0.5, label='shock')
    lgd = ax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.00, 0.50), 
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    if meta:
        markings = [po, pi, p, paths]
        return fig, markings
    #I/O
    if not batch:
        return fig
    else:
        writeFig(fig, paths, filetag)


def flashProfile(fname, thresh=1e-6, xrange=[0.0, 0.0], yrange=[0.0, 0.0],
                 filetag='prof', batch=False, byM=True, direction=[], grid=False):
    """Plot bulk properties and species in a chekpoint through a ray.
    
    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction.
        xrange(float list): abscissa range.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        direction(float list): list of spherical angles (alt, azimuth), empty for 1D.
        grid(bool): add radial positions as a line grid.
    
    Returns:
        (mpl figure) or (None)
    
    """
    time, pars, paths, prof = fetchData(fname, direction)
    fig = plotDMatMerged(prof, byM=byM, thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.82, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    if direction:
        b = ax.annotate("Angle(s): {}".format(','.join([str(i) for i in direction])),
                        xy=(0.0, 0.0), xytext=(0.82, 0.05), size=12,
                        textcoords='figure fraction', xycoords='figure fraction', 
                        bbox=dict(boxstyle='round', fc='w', ec='k'))
    if sum(yrange) != 0.0:
        ax.set_ylim(yrange)
    if grid:
        for ax in fig.axes:
            drawGrid(ax, prof.radius)
    if not batch:
        return fig
    else:
        writeFig(fig, paths, filetag)


def flashDegeneracy(fname, thresh=1e-6, filetag='deg', batch=False, 
                    byM=True, direction=[],  xrange=[0.0, 0.0]):
    time, pars, paths, prof = fetchData(fname, direction)
    fig = plotDegen(prof, byM=byM, 
                   thresh=thresh, xrange=xrange)
    ax = plt.gca()
    a = ax.annotate("{:.5f} s".format(time),
                    xy=(0.0, 0.0), xytext=(0.88, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    if direction:
        b = ax.annotate("Angle(s): {}".format(','.join([str(i) for i in direction])),
                        xy=(0.0, 0.0), xytext=(0.82, 0.05), size=12,
                        textcoords='figure fraction', xycoords='figure fraction', 
                        bbox=dict(boxstyle='round', fc='w', ec='k'))
    if not batch:
        return fig
    else:
        writeFig(fig, paths, filetag)
    

def flashSpecies(fname, thresh=1e-6, filetag='spec', batch=False, xrange=[0.0, 0.0],
                 byM=True, direction=[], vmax=4e9, plotall=False):
    """Plot species and aggregated masses in a chekpoint through a ray.
    Masses are a trapezoid integrated ray for each species.
    
    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction and mass yields.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
        direction(float list): list of spherical angles (alt, azimuth), empty for 1D.
        vmax(floar): limit for species velocities.
        plotall(bool): force plotting every species found.
    
    Returns:
        (mpl figure) or (None)
    
    """
    fields = ['density', 'temperature', 'pressure', 'velx']
    time, pars, paths, prof = fetchData(fname, direction, fields)
    fig = plt.figure(figsize=(10, 8))
    layout = (3,3)
    # plot species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(ax1, prof, byMass=byM, thresh=thresh, plotall=plotall)
    if byM:
        ax1.set_xlabel('Mass ($M_{\odot}$)')
    else:
        ax1.set_xlabel('Radius (cm)')
    if sum(xrange) != 0.0:
        ax1.set_xlim(xrange)
    # timestamp and legend
    a = ax1.annotate("{:.5f} s".format(time),
                     xy=(0.0, 0.0), xytext=(0.65, 0.1), size=12,
                     textcoords='figure fraction', xycoords='figure fraction', 
                     bbox=dict(boxstyle='round', fc='w', ec='k'))
    lgd = ax1.legend(ncol=4, loc='upper left', bbox_to_anchor=(1.05, 1.0), 
                     columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                     numpoints=3, handletextpad=0.0)
    lgd.get_frame().set_edgecolor('k')
    
    # write otp files for yields
    massfiletuples = []
    
    # plot masses
    ax2 = plt.subplot2grid(layout, (1,0), colspan=2)
    for i, sp in enumerate(prof.species):
        props = next(ax2._get_lines.prop_cycler)
        spmass = trapz(getattr(prof, sp), x=prof.masses)
        massfiletuples.append((sp, spmass))
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
    
    # decay yield
    names, decmasses = decayYield(*zip(*massfiletuples))

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
        mpln, mplbins, patches = ax3.hist(bins, bins=len(bins), weights=weights, 
                                          histtype='step', log=True, label=sp)
    ax3.set_ylim([thresh, 2])
    ax3.set_xlim([0.0, vmax])
    ax3.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.xaxis.set_major_formatter(customFormatter(9))
    ax3.xaxis.set_minor_formatter(StrMethodFormatter(''))
    ax3.set_xlabel('Speed ($10^{9}$ cm/s)')
    fig.subplots_adjust(hspace=0.4)
    
    if not batch:
        return fig
    else:
        dest, num = writeFig(fig, paths, filetag)
        yieldfile = os.path.join(dest, '{}{}.yield'.format(filetag, num))
        with open(yieldfile, 'w') as f:
            f.write('# {}\n'.format(paths[1]))
            f.write('# time: {}\n'.format(time))
            f.write('\n'.join(['{} {}'.format(*a) for a in massfiletuples]))
            f.write('\n')
        print("Wrote: {}".format(yieldfile))
        decayfile = os.path.join(dest, '{}{}.decayed'.format(filetag, num))
        with open(decayfile, 'w') as f:
            f.write('# {}\n'.format(paths[1]))
            f.write('# time: {}\n'.format(time))
            f.write('\n'.join(['{} {}'.format(*a) for a in zip(names, decmasses)]))
            f.write('\n')
        print("Wrote: {}".format(decayfile))


def plainTprofile(fname, thresh=1e-4, xrange=[0.0, 0.0], byM=True, merged=False):
    """plots main properties of a plain text file.
    
    Args:
        fname (str): filename of checkpoint.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """
    prof = dataMatrix(fname)
    if merged:
        return plotDMatMerged(prof, thresh=thresh, xrange=xrange, byM=byM)
    else:
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
    ncol = 4
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


def plotDMatMerged(prof, thresh=1e-4, xrange=[0.0, 0.0], byM=True, alpha=1.0):
    """plots main properties of a profile in only two axes, merging 
    thermodynamic properties.
    
    Args:
        prof (dataMatrix): dataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    """  
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    # DEBUG
#     print(len(prof.species))
#     print(prof.t3)
    keys = sortNuclides(prof.species)
    ncol = 4
    labelspace = -0.1
    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False  
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (2, 2)
    # Species
    spax = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(spax, prof, byMass=byM, thresh=thresh, plotall=False, alpha=alpha)
    # remove last(lowest) yticklabel to avoid overlap
    spax.get_yticklabels()[1].set_visible(False)
    lgd = spax.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.02), 
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    spax.axhline(1e0, linewidth=1, linestyle=':', color='black')
    spax.set_ylim(thresh, 2.0)
    spax.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    spax.yaxis.set_label_coords(labelspace, 0.5)
    spax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    spax.tick_params(labelbottom=False)
    spax.set_xlabel('')
    if log:
        spax.set_yscale('log')
        spax.set_xscale('log')
    # Thermodynamical variables (reference is, as always, density)
    tdax = plt.subplot2grid(layout, (1, 0), colspan=2, sharex=spax)
    tdax.semilogy(xs, prof.density, label='Density')
    ylabels = ['$\\frac{g}{cm^3}$']
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.tick_params(labelbottom=False)
    
    for i, p in enumerate(plotp):
        u = ut.getUnit(p)
        if p=='pressure':
            tdax.plot(xs, getattr(prof, p)/1e16, label=p.capitalize())#, color='k')
            ylabels.append(u+'$\\times10^{-16}$')
        else:
            tdax.plot(xs, getattr(prof, p), label=p.capitalize())#, color='k')
            ylabels.append(u)
        #spacer = labelspace if '\\' in u else labelspace - 0.02
    #tdax.set_ylabel('{}({})'.format(p.capitalize(), u), rotation=90, size=13, labelpad=0)
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.tick_params(labelbottom=False)
    tdax.set_xlabel(xlab)
    tdax.tick_params(labelbottom=True)
    tdax.set_ylabel('({})'.format(','.join(ylabels)), size=13, rotation=90, labelpad=0)
    tdax.yaxis.set_label_coords(labelspace, 0.5)
    lgd = tdax.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.00, 0.50), 
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    if sum(xrange)!=0.0:
        spax.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.01)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotDegen(prof, thresh=1e-4, xrange=[0.0, 0.0], byM=True):
    """plots composition and degeneracy parameters for a lineout,
    namely Y_e and \eta = T / T_fermi.
    \eta near zero implies degeneracy, 
    while \eta >> 1 or negative implies non-degenerate matter.
    
    Args:
        prof (dataMatrix): dataMatrix obj.
        thresh (float): ymin for species fraction plot.
        xrange (list of float): if set, change the range of the plots.
        byM (bool): abundance plot xaxis (by Mass or by Radius).
    
    Returns:
        (mpl figure) or (None)
    
    """  
    fig = plt.figure()
    skip = ['radius', 'masses', 'density']
    plotp = [x for x in prof.bulkprops if x not in skip]
    keys = sortNuclides(prof.species)
    ncol = 4
    labelspace = -0.1
    if byM:
        xs = prof.masses
        xlab = 'Mass ($M_{\odot}$)'
        log = False  
    else:
        xs = prof.radius
        xlab = 'Radius ($cm$)'
        log = True

    layout = (2, 2)
    # Species
    spax = plt.subplot2grid(layout, (0, 0), colspan=2)
    skip = plotSpecies(spax, prof, byMass=byM, thresh=thresh, plotall=False)
    # remove last(lowest) yticklabel to avoid overlap
    spax.get_yticklabels()[1].set_visible(False)
    lgd = spax.legend(ncol=ncol, loc='upper left', bbox_to_anchor=(1.00, 1.02), 
                      columnspacing=0.0, labelspacing=0.0, markerfirst=False, 
                      numpoints=3, handletextpad=0.0, edgecolor='k')
    spax.axhline(1e0, linewidth=1, linestyle=':', color='black')
    spax.set_ylim(thresh, 2.0)
    spax.set_ylabel('Mass Frac.($X_{i}$)', size=13, rotation=90, labelpad=0)
    spax.yaxis.set_label_coords(labelspace, 0.5)
    spax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    spax.tick_params(labelbottom=False) 
    if log:
        spax.set_yscale('log')
        spax.set_xscale('log')
    # Ye and Yion
    yes, yis = [], []
    npnts = len(getattr(prof, prof.species[0]))
    for i in range(npnts):
        xis = []
        for s in prof.species:
            xis.append(getattr(prof, s)[i])
        invyi, invye = getMus(prof.species, xis)
        yes.append(1.0/invye)
        yis.append(1.0/invyi)
    
    yeax = plt.subplot2grid(layout, (1, 0), colspan=2, sharex=spax)
    pl1 = yeax.semilogy(xs, yes, label='$Y_e$', color='#0082c8', alpha=0.8)
    yeax.set_ylabel('$Y_e$', size=13, rotation=0, labelpad=10)
    yeax.yaxis.set_major_formatter(customFormatter(0, prec=4, width=4))
    offset = yeax.get_yaxis().get_offset_text()
    offset.set_visible(False)
    yeax.set_ylim(0.1, 1.2)
    # yeax.get_yticklabels()[1].set_visible(False)
    # print(yeax.get_yticklabels())
    # yeax.get_yticklabels()[0].set_visible(False)
    # temps
    tdax = yeax.twinx()
    # get fermi temperatures through the lineout
    # fermT = [ extRelFermi(d)/ut.kb for d in prof.density ]
    # tdax.plot(xs, prof.temperature/fermT, label='extreme fermT')#, color='k')
    fermT = [ nonRelFermi(d, ye=y)/ut.kb for d, y in zip(prof.density, yes) ]
    pl2 = tdax.semilogy(xs, prof.temperature/fermT, label=r'$\eta=\frac{T}{T_{f}}$', color='#3cb44b', alpha=0.8)
    tdax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    tdax.set_ylabel('$\eta$', size=13, rotation=0, labelpad=8)
    plts = pl1 + pl2
    labels = [ l.get_label() for l in plts ]
    lgd = tdax.legend(plts, labels, ncol=1, loc='upper left', labelspacing=0.0, markerfirst=False, 
                      numpoints=3, handletextpad=0.0, frameon=False)
    yeax.set_xlabel(xlab)
    if sum(xrange)!=0.0:
        spax.set_xlim(xrange)
    plt.subplots_adjust(hspace=0.001, wspace=0.0)
    plt.subplots_adjust(left=0.13, right=0.80)
    plt.subplots_adjust(top=0.99, bottom=0.10)
    fig.set_size_inches(8.5, 7.5, forward=True)
    return fig


def plotSpecies(ax, dmatr, byMass=True, thresh=1e-4, plotall=False, alpha=1.0):
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
            line = simplePlot(ax, dmatr, absc, s, log=True, alpha=alpha)
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
                                  color=props['color'], ls=props['linestyle'], alpha=alpha)
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

def drawGrid(ax, gridpoints, alpha=0.6, color='salmon', lw=2.0):
    """draws gridpoints as a line grid on an axes."""
    print('Gridpoints: {}'.format(len(gridpoints)))
    for r in gridpoints:
        ax.axvline(r, alpha=alpha, color=color, lw=lw)