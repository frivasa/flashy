"""nuclide aggregate plots, based on yield files."""
from ..nuclear import (split_species, convertYield2Abundance,
                       getMassFractions, getAbundances,
                       Avogadro, readSunComp, decayYield,
                       readYield, sortNuclides,
                       elemSplit, AGSS09, readIsotopicSolar)
from .globals import (np, os, mpl, plt, writeFig, io, log,
                      StrMethodFormatter, LogNorm)
from matplotlib.patches import Rectangle
import flashy.utils as ut
from ..datahaul.hdf5yt import getYields, getLineout
from ..datahaul.hdfdirect import directMeta, getBlockLineout


def flash_elementalYields(fname, tag='elem_solar', yrange=[-3, 15],
                          norm='Si', offset=6, batch=False):
    """plots decayed elemental yields from a checkpoint vs sun.

    Args:
        fname(str): filename path.
        tag(str): output tag.
        yrange(int list): abundance range.
        norm(str): abundance reference.
        offset(int): abundance scale offset.
        batch(bool): write to file toggle.

    Returns:
        (mpl.figure)

    """
    fig, ax = plt.subplots(figsize=(10, 5))
    name = os.path.dirname(os.path.dirname(fname)).split('/')[-2:]
    title = "/".join(name)
    time, species, masses = getYields(fname)
    # decay yields so that they're comparable to the Sun
    species, masses = decayYield(species, masses)
    fakef = ''
    for s, m in zip(species, masses):
        fakef += '{} {:.10f}\n'.format(s, m)
    ryield = readYield(io.StringIO(fakef))
    plabels = []
    plabels.append(plotAbun(ax, ryield, norm=norm, offset=offset,
                            label=title, tagspecies=True, yrange=yrange))
    # plot sun
    zs, names, mix = readSunComp(AGSS09)
    massF = getMassFractions(names, mix)
    p = ax.plot(zs, getAbundances(names, massF, scale=norm, offset=offset),
                marker='x', markersize=3, ls=':', lw=1, color='brown',
                label='Sun(AGSS09)')
    ax.set_ylim(yrange)
    ax.set_title("{:.5f} s".format(time), loc='left')
    lg = ax.legend(loc='upper right')
    if batch:
        writeFig(fig, fname, tag)
    else:
        return fig


def elementalYields(yieldfiles, tags=[], norm='Si', offset=6):
    """plots abundances for a list of files
    or a single file vs solar composition.

    Args:
        yieldfiles(str list): yield filenames.
        tags(str list): yield source names.
        norm(str): normalization element.
        offset(float): log of the abundance offset
            (generally 12 for H, 6 for Si).

    Returns:
        (mpl.figure)

    """
    fig, ax = plt.subplots(figsize=(10, 5))
    if isinstance(yieldfiles, list):
        if not tags:
            print('no tags specified.')
            tags = [None]*len(yieldfiles)
        plabels = []
        for i, f in enumerate(yieldfiles):
            ryield = readYield(f)
            props = next(ax._get_lines.prop_cycler)
            col = props['color'] if i else 'k'
            tagspecies = True if not i else False
            plabels.append(plotAbun(ax, ryield, norm=norm, offset=offset,
                                    color=col, label=tags[i],
                                    tagspecies=tagspecies))
            lg = ax.legend(*zip(*plabels))
        # speciesGrid(ax, modspecies, yoffset=5.0)
    else:
        plabels = []
        ryield = readYield(yieldfiles)
        tag = tags
        plabels.append(plotAbun(ax, ryield, norm=norm,
                                offset=offset, label=tag, tagspecies=True))
        # plot sun
        zs, names, mix = readSunComp(AGSS09)
        massF = getMassFractions(names, mix)
        p = ax.plot(zs, getAbundances(names, massF, scale=norm, offset=offset),
                    marker='x', markersize=3, ls=':', lw=1, color='brown',
                    label='Sun(AGSS09)')
        lg = ax.legend()
    return fig


def flash_productionFac(fnames, tag='netname', norm='Si',
                        offset=6, batch=False):
    fig, ax = plt.subplots(figsize=(10, 5))
    if isinstance(fnames, list):
        if len(fnames) != 2:
            print('Input files must be exactly two.')
            return None
        plabels = []
        name = os.path.dirname(os.path.dirname(fnames[0]))
        name = name.split('/')[-2:]
        time, species, masses = getYields(fnames[0])
        # decay yields so that they're comparable to the Sun
        species, masses = decayYield(species, masses)
        fakef = ''
        for s, m in zip(species, masses):
            fakef += '{} {:.10f}\n'.format(s, m)
        ryield = readYield(io.StringIO(fakef))
        plabels.append(plotPfac(ax, ryield, refname=fnames[1],
                                norm=norm, offset=6, label=tag,
                                reftype='checkpoint', tagspecies=True))
        lg = ax.legend()
    else:
        time, species, masses = getYields(fnames)
        # decay yields so that they're comparable to the Sun
        species, masses = decayYield(species, masses)
        ryield = readYield(zip(species, masses))
        plabels = []
        label = '{} vs AGSS09'.format(tag)
        plabels.append(plotPfac(ax, ryield, label=label,
                                norm=norm, offset=offset, tagspecies=True))
        lg = ax.legend()
    if batch:
        writeFig(fig, fname, 'pfac_' + tag)
    else:
        return fig


def productionFactor(yieldfiles, tag='Sim vs X', norm='Si', offset=6):
    """plots profuction factors for comparable
    species in a pair of yields, or vs solar compositon
    for a single file.
    yieldfiles = [input, reference]

    Args:
        yieldfiles(str list): yield filenames.
        tags(str list): yield source names.
        norm(str): normalization element.
        offset(float): log of the abundance offset
            (generally 12 for H, 6 for Si).

    Returns:
        (mpl.figure)

    """
    fig, ax = plt.subplots(figsize=(10, 5))
    if isinstance(yieldfiles, list):
        if len(yieldfiles) != 2:
            print('Input files must be exactly two.')
            return None
        plabels = []
        ryield = readYield(yieldfiles[0])
        plabels.append(plotPfac(ax, ryield, refname=yieldfiles[1],
                                norm=norm, offset=6, label=tag,
                                reftype='yield', tagspecies=True))
        lg = ax.legend()
        # speciesGrid(ax, modspecies, yoffset=5.0)
    else:
        plabels = []
        ryield = readYield(yieldfiles)
        plabels.append(plotPfac(ax, ryield, label='Sim vs AGSS09',
                                norm=norm, offset=offset, tagspecies=True))
        lg = ax.legend()
    return fig


def multi_nuclideYields(fnames, tags=['A'], xrange=[0, 70], yrange=[1e-5, 1.0],
                        addsun=True, filetag='nuclide_solar', batch=False, subsetR=0.0):
    """plots all nuclide yields in a checkpoint list. Adds sun as option.
    deprecates 'flash_nuclideYields' since this can be used for a single file.

    Args:
        fnames(str list): filename paths.
        tags(str list): custom label.
        xrange(int list): nuclear mass range.
        yrange(float list): solar mass range.
        addsun(str): plot sun composition.
        thresh(float): mass threshold for plot.
        tag(str): output tag.
        batch(bool): write to file toggle.
        subsetR(float): sum masses <= subsetR (0.0=whole domain).

    Returns:
        (mpl.figure)

    """
    fig, ax = plt.subplots(figsize=(10, 5))
    plabels = []
    syms = ['.', 'x', '+', '1', '2']
    cols = ['black', 'royalblue', 'forestgreen', 'darkorange', 'firebrick']
    if len(fnames) > len(tags):
        tags = [tags[0]]*len(fnames)
    for i, fname in enumerate(fnames):
        time, species, masses = getYields(fname, subsetR=subsetR)
        print(fname)
        if not i:
            lines = ["# {:.8f}".format(time)]
            ltag = "{} {:.8f}"
            lines += [ltag.format(x, y) for x,y in zip(species, masses)]
            print("\n".join(lines))
            meta = "\n".join(lines)
        title = "{} ({:.3f} s)".format(tags[i], time)
        fakef = ''
        for s, m in zip(species, masses):
            fakef += '{} {:.10f}\n'.format(s, m)
        ryield = readYield(io.StringIO(fakef))
        pp = next(ax._get_lines.prop_cycler)
        tag = True if i else False
        plabels.append(plotIsoMasses(ax, ryield, notag=tag, ls='None',
                                     ylims=yrange, xlims=xrange,
                                     marker=syms[i%len(syms)],
                                     color=cols[i%len(cols)],
                                     label=title))
    if addsun:
        ryield2 = readIsotopicSolar()
        plabels.append(plotIsoMasses(ax, ryield2, notag=False,
                                     ylims=yrange, xlims=xrange,
                                     label='AGSS09 (Solar)',
                                     color='#FF8C00'))
    ax.set_xlim(xrange)
    ax.set_ylim(yrange)
    legdict = {'ncol': 1, 'loc': 'best', 'columnspacing': 0.0,
               'labelspacing': 0.1, 'numpoints': 2, 'handletextpad': 0.2}
    lg = ax.legend(*zip(*plabels), **legdict)
    if batch:
        stem = fnames[0]
        writeFig(fig, stem, filetag, meta)
    else:
        return fig


# def flash_nuclideYields(fname, tag='nuclide_solar', xrange=[0, 70],
#                         compfile='sun', thresh=1e-5, ratio=False, batch=False):
#     """plots all nuclide yields in a checkpoint vs sun or any other
#     simulation specified by its runfolder.
#
#     Args:
#         fname(str): filename path.
#         tag(str): output tag.
#         xrange(int list): nuclear mass range.
#         compfile(str): comparison chekpoint path.
#         thresh(float): mass threshold for plot.
#         ratio(bool): print out compariosn ratios.
#         batch(bool): write to file toggle.
#
#     Returns:
#         (mpl.figure)
#
#     """
#     fig, ax = plt.subplots(figsize=(10, 5))
#     name = os.path.dirname(os.path.dirname(fname)).split('/')[-2:]
#     title = "/".join(name)
#     time, species, masses = getYields(fname)
#     fakef = ''
#     for s, m in zip(species, masses):
#         fakef += '{} {:.10f}\n'.format(s, m)
#     ryield = readYield(io.StringIO(fakef))
#     plabels = []
#     plabels.append(plotIsoMasses(ax, ryield, notag=False,
#                                  ylims=[thresh, 1.0],
#                                  label=title, color='k'))
#     if ratio:
#         vals = []
#         for k in ryield.keys():
#             zz = ryield[k]['z']
#             for n in ryield[k]['n'].keys():
#                 vals.append((k, n + zz, ryield[k]['n'][n]))
#         tags = ['{} {} {}'.format(*v) for v in vals]
#         print('\n'.join(tags))
#     if compfile == 'sun':
#         ryield2 = readIsotopicSolar()
#         plabels.append(plotIsoMasses(ax, ryield2, notag=False,
#                                      ylims=[thresh, 1.0], label='AGSS09',
#                                      color='orange'))
#         mtitle = "{:.5f} s".format(time)
#     else:
#         name = os.path.dirname(os.path.dirname(compfile))
#         name = name.split('/')[-2:]
#         title = "/".join(name)
#         time2, species, masses = getYields(compfile)
#         mtitle = "{:.3f}s vs\n {:.3f}s".format(time, time2)
#         fakef = ''
#         for s, m in zip(species, masses):
#             fakef += '{} {:.10f}\n'.format(s, m)
#         ryield2 = readYield(io.StringIO(fakef))
#         plabels.append(plotIsoMasses(ax, ryield2, notag=False,
#                                      label=title, color='#e68053'))
#     if ratio:
#         vals2 = []
#         for k in ryield2.keys():
#             zz = ryield2[k]['z']
#             for n in ryield2[k]['n'].keys():
#                 vals2.append((k, n + zz, ryield2[k]['n'][n]))
#         ratios = [v[-1]/vv[-1] for (v, vv) in zip(vals, vals2)]
#         tags = ['{} {} {}'.format(*a[:2], b) for (a, b) in zip(vals, ratios)]
#         print('\n'.join(tags))
#     ax.set_xlim(xrange)
#     ax.set_ylim([thresh, 2])
#     ax.set_title(mtitle, loc='left')
#     legdict = {'ncol': 1, 'loc': 'upper right', 'columnspacing': 0.0,
#                'labelspacing': 0.1, 'numpoints': 3, 'handletextpad': 0.2,
#                'bbox_to_anchor': (1.0, 1.15)}
#     lg = ax.legend(*zip(*plabels), **legdict)
#     if batch:
#         writeFig(fig, fname, tag)
#     else:
#         return fig


def nuclideYields(files, tags):
    """plot mass fractions for each nuclide in a yield file.

    Args:
        files(str list or str): yield filename(s).
        tags(str list or str): yield source name(s).

    Returns:
        (mpl.figure)

    """
    if isinstance(files, list):
        fig, ax = plt.subplots(figsize=(10, 5))
        plabels = []
        for i, (f, t) in enumerate(zip(files, tags)):
            ryield = readYield(f)
            props = next(ax._get_lines.prop_cycler)
            col = props['color'] if i else 'k'
            plabels.append(plotIsoMasses(
                ax, ryield, notag=False, label=t, color=col))
        lg = ax.legend(*zip(*plabels))
        return fig
    else:
        fig, ax = plt.subplots(figsize=(10, 5))
        plabels = []
        ryield = readYield(files)
        plabels.append(plotIsoMasses(
            ax, ryield, notag=False, label=tags, color='k'))
        ryield = readIsotopicSolar()
        plabels.append(plotIsoMasses(ax, ryield, notag=False,
                                     label='AGSS09', color='orange'))
        lg = ax.legend(*zip(*plabels))
        return fig


def plotFlowMatrix(fname, dpi=100, edge=1e34,
                   linthresh=1e32, batch=False):
    """RESTRUCT to .dat output
    builds a flow matrix with a given checkpoint and its succesor in number.

    Args:
        fname(str): checkpoint to be used.
        dpi(int): dpi of output figure.
        edge(float): limits of flow of particles.
        linthresh(float): SymLog scale "start"
            (-edge  -linthresh 0 +linthresh +edge).
        batch(bool): write to file(true) or return figure(false).

    Returns:
        (mpl.figure or None)

    """
    # get the mass yields (corrected for volume and density) for both files
    t1, species, masses1 = getYields(fname)
    shift = '{}{:04d}'.format(fname[:-4], int(fname[-4:])+1)
    # print(shift)
    t2, _, masses2 = getYields(shift)
    deltaT = t2 - t1
    deltaM = np.array(masses2) - np.array(masses1)
    # convert to particles via species masses
    names, prot, neut, nucmass = split_species(species, standardize=True)
    fluxes = ut.msol*deltaM/nucmass*Avogadro
    f, a = plt.subplots(dpi=dpi)
    sq = plotNuclideGrid(a, species, xmass=fluxes, sym=True, time=t1,
                         cmmin=-edge, cmmax=edge,
                         cmap='coolwarm', linthresh=linthresh)
    cax = f.add_axes([0.12, 0.9, 0.8, 0.037])  # left, bot, width, height
    norm = mpl.colors.SymLogNorm(vmin=-edge, vmax=edge, linthresh=linthresh)
    cbar1 = mpl.colorbar.ColorbarBase(
        cax, cmap='coolwarm', norm=norm, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    datablock, sps = getLineout(fname, species=False, geom='spherical')
    dens = datablock[1][0]
    temp = datablock[2][0]

    imgtag = u'$\\rho$ = {:.2E}\nT = {:.2E}'.format(dens, temp)
    # trick to make an empty handle
    sq = plt.Rectangle((0, 0), 1, 1, fill=False,
                       edgecolor='none', visible=False)
    leg = a.legend(handles=[sq], labels=[imgtag],
                   loc=(0.0, 0.83), frameon=False)

    if batch:
        writeFig(f, fname, 'flows')
    else:
        return f


def plotGridYield(yieldfile, dpi=150, cmap='Oranges',
                  filetag='gridspec', imgtag='', batch=False):
    """plots mass yields over whole species grid for a yield file.

    Args:
        yieldfile(str): path of file.
        dpi(float): dpi of figure.
        cmap(str): colormap name.
        filetag(str): prefix for batch mode.
        imgtag(str): label for plot (appears within the axes)
        batch(bool): skips returning figure,
        saving it to a structured directory instead.

    Returns:
        (mpl figure) or (None)

    """
    names = np.genfromtxt(yieldfile, dtype='|U4', usecols=(0,))
    vals = np.genfromtxt(yieldfile, usecols=(1,))
    mass = np.sum(vals)
    with open(yieldfile, 'r') as tf:
        allrows = tf.readlines()
    time = float(allrows[1].split()[-1])
    f, a = plt.subplots(dpi=dpi)
    sq = plotNuclideGrid(a, names, xmass=vals, cmap=cmap, time=time)
    cax = f.add_axes([0.12, 0.9, 0.8, 0.037])  # left, bot, width, height
    cmap = mpl.cm.Oranges
    # norm = mpl.colors.Normalize(vmin=0.0, vmax=0.5)
    norm = mpl.colors.LogNorm(vmin=1e-3, vmax=0.5)
    cbar1 = mpl.colorbar.ColorbarBase(
        cax, cmap=cmap, norm=norm, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    if imgtag:
        # trick to make an empty handle
        sq = plt.Rectangle((0, 0), 1, 1, fill=False,
                           edgecolor='none', visible=False)
        leg = a.legend(handles=[sq], labels=[imgtag],
                       loc=(0.0, 0.85), frameon=False)
    if batch:
        folder = os.path.dirname(yieldfile)
        filename = os.path.basename(yieldfile).split('.')[0]
        paths = (folder, filename)
        writeFig(f, paths, filetag)
    else:
        return f


def plotNuclideGrid(ax, species, xmass=[], time=0.0, z_range=[-2, 35],
                    n_range=[-2, 38], boxsize=1.0, cmmin=1.0e-5,
                    cmmax=0.9, cmap='Blues', addtags=True, tagsonright=False,
                    noedge=False, frameless=False, forceColor=False,
                    fColor='peru', sym=False, linthresh=1e10):
    """Plots a list of species on a grid.

    Args:
        ax(mpl.axes): axes to draw on.
        species(str list): list of nuclide codes (e.g.: Ca40).
        xmass(float list): list of mass fractions for each nuclide.
        time(float): timestamp for plot (for use with xmass).
        z_range(float list): atomic number range for plot.
        n_range(float list): neutron number range for plot.
        boxsize(float): size of nuclide marker box.
        cmmin(float): colormap min value.
        cmmax(float): colomap max value.
        cmap(str): colormap name.
        addtags(bool): print element names.
        tagsonright(bool): move tags to the rightside of the network.
        noedge(bool): remove nuclide marker edges.
        frameless(bool): remove frames, only show arrows.
        forceColor(bool): force a single color for all species.
        fColor(str): Hex or name of color to force.
        sym(bool): use symmetrical log (for neg values).
        linthresh(float): linthresh for symlognorm.

    Returns:
        (mpl.rectangle): handle for mpl legend.

    """
    cmap = mpl.cm.get_cmap(cmap)
    if sym:
        norm = mpl.colors.SymLogNorm(vmin=cmmin, vmax=cmmax,
                                     linthresh=linthresh)
    else:
        norm = mpl.colors.LogNorm(vmin=cmmin, vmax=cmmax)
    if len(xmass) > 0:
        xis = xmass
        a = ax.annotate("{:.5f} s".format(time),
                        xy=(0.0, 0.0), xytext=(0.72, 0.18), size=12,
                        textcoords='figure fraction',
                        xycoords='figure fraction')
        # , bbox=dict(boxstyle='round', fc='w', ec='k'))
    else:
        xis = [1e-2]*len(species)
    if noedge:
        clr = 'None'
    else:
        clr = 'black'
    nam, zs, ns, As = split_species(species, standardize=True)
    for (z, n, xi) in zip(zs, ns, xis):
        if forceColor:
            facecolor = fColor
        else:
            facecolor = cmap(norm(xi))
        square = plt.Rectangle((n-0.5*boxsize, z-0.5*boxsize),
                               boxsize, boxsize, facecolor=facecolor,
                               edgecolor=clr)
        ax.add_patch(square)
        mainsquare = square
    if addtags:
        tags, xs, ys = getNuclideGridNames(
            list(nam), zs, ns, rightside=tagsonright)
        for (t, x, y) in zip(tags, xs, ys):
            if tagsonright:
                ax.text(x, y, t, fontsize=8, verticalalignment='center',
                        horizontalalignment='left')
            else:
                ax.text(x, y, t, fontsize=8, verticalalignment='center',
                        horizontalalignment='right')

    ax.set_ylabel('Z', rotation=0, labelpad=12)
    ax.set_xlabel('N')
    ax.set_ylim(z_range)
    ax.set_xlim(n_range)

    # remove spines/ticks (arrow and named axes)
    if frameless:
        for side in ['bottom', 'right', 'top', 'left']:
            ax.spines[side].set_visible(False)
        ax.set_xticks([])  # labels
        ax.set_yticks([])
        ax.arrow(0, 12, 0, 10, head_width=0.5, head_length=1, fc='k', ec='k')
        ax.arrow(13, 0, 10, 0, head_width=0.5, head_length=1, fc='k', ec='k')

    # numbered axes
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.xaxis.set_minor_formatter(StrMethodFormatter(''))
    return mainsquare


def plotReacNet(ax, sunet, matr_shape, step=4, xoffset=0,
                color='#FF8200', tagsize=8, dotsize=15, aspect=1.0):
    """plots a reaction network based on
    default files at cdx/Network/netname

    Args:
        ax(mpl.axes): axes to draw on.
        sunet(str): filename with species names.
        matr_shape(str): filename with reaction tags.
        step(int): element name skip size.
        xoffset(int): element name tag start offset.
        color(str): color for dots.
        tagsize(int): font size for element tags.
        dotsize(int): plot circle size.
        aspect(float): overall ratio of axes.

    Returns:
        (int, int): number of species and number of rates found.

    """
    with open(matr_shape) as f:
        for i, line in enumerate(f):
            if len(line) != 22:
                continue
            else:
                break
    data = np.genfromtxt(matr_shape, skip_header=i)
    x, y, z = np.hsplit(data, 3)
    with open(sunet) as f:
        species = f.readlines()
    species = [s.strip() for s in species]
    names = [elemSplit(s, invert=True) for s in species]
    nsp = len(names)
    rates = len(x)
    im = ax.scatter(x, y, c=color, edgecolors='#000000', s=dotsize)
    ax.set_aspect(aspect)
    # set labels
    ax.set_yticks(np.arange(1, nsp+1, 1))  # all tags in y-axis
    ax.set_yticklabels(['$^{{{}}}{{{}}}$'.format(*t) for t in names],
                       fontsize=tagsize, ha='right')
    ax.set_xticks(np.arange(xoffset+1, nsp+1, step))
    labels = names[xoffset::step]
    ax.set_xticklabels(['$^{{{}}}{{{}}}$'.format(*t) for t in labels],
                       fontsize=tagsize, va='baseline')
    ax.tick_params('both', length=5, width=2, which='major')
    ax.tick_params('x', length=5, width=2, labeltop=True, labelbottom=False)
    t = '{} Isotopes\n{} Rates'
    note = ax.annotate(t.format(nsp, rates), xy=(100, 30), fontsize=8)
    plt.gca().invert_yaxis()
    return nsp, rates


def doYouBelieveInMagic(ax, color='brown'):
    """Draws a magical grid of stability.

    Args:
        ax(mpl.axes): axes to draw upon.
        color(str): color of grid.

    Returns:
        (mpl.line2D): last grid line for legend use.

    """
    cadabra = [2, 8, 20, 28, 50, 82, 126]
    for abra in cadabra:
        ax.axhline(abra, alpha=0.4, color=color, ls=':')
        p = ax.axvline(abra, alpha=0.4, color=color, ls=':')
    return p


def plotPfac(ax, querym, refname=AGSS09, label='Sun vs Ref',
             norm='Si', offset=6, reftype='solar', tagspecies=False):
    """draws abundquery/abundref from a massdict and a filename,
    types of reference are 'solar'(for solar composition)
    and 'yield' (for other sims)
    returns label and line element for legend

    Args:
        ax(mpl.axes): axes to draw upon.
        querym(str): yield file to use.
        refname(str): solar yield file to use.
        label(str): tag for plot.
        norm(str): reference element.
        offset(float): log of the abundance offset
            (generally 12 for H, 6 for Si).
        reftype(str): is the ref 'solar' or a pure 'yield'?
        tagspecies(bool): plot species names on plot.

    Returns:
        (mpl.line2D[0], str): line object for legend, tag of plot.

    """
    zs, ns, ab = convertYield2Abundance(querym, norm=norm, offset=offset)
    if reftype == 'solar':
        zss, nss, solab = readSunComp(refname)
        massF = getMassFractions(nss, solab)
        rescaledsolab = getAbundances(nss, massF, scale=norm, offset=offset)
        soldict = dict(zip(nss, rescaledsolab))
    elif reftype == 'yield':
        massdict = readYield(refname)
        zss, nss, solab = convertYield2Abundance(
            massdict, norm=norm, offset=offset)
        soldict = dict(zip(nss, solab))
    elif reftype == 'checkpoint':
        time, species, masses = getYields(refname)
        # decay yields so that they're comparable to the Sun
        species, masses = decayYield(species, masses)
        massdict = readYield(zip(species, masses))
        zss, nss, solab = convertYield2Abundance(
            massdict, norm=norm, offset=offset)
        soldict = dict(zip(nss, solab))
    pfacs = []
    for i, n in enumerate(ns):
        if n in soldict:
            pfacs.append((zs[i], n, ab[i]-soldict[n]))
            # pfacs.append((zs[i], n, ab[i]/soldict[n]))
        else:
            print('{} not found in ref.'.format(n))
    x, ns, y = zip(*sorted(pfacs))
    ax.axhline(0, ls=':', lw=1, color='green')
    line = ax.plot(x, y, label=label, ls='--',
                   lw=0.5, marker='.', color='black')
    # Prettify
    if tagspecies:
        for x, n, y in sorted(pfacs):
            ax.text(x, y, '${}$'.format(n), color='black',
                    size=8, horizontalalignment='right',
                    verticalalignment='bottom')
    ax.set_xlabel(u'Atomic Number (Z)')
    ax.set_ylabel('$[X/{}]- [X/{}]_{{ref}}$'.format(norm, norm))
    ax.autoscale()
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def plotIsoMasses(ax, mdict, label='W7', color='black',
                  xlims=[-2, 78], ylims=[1e-18, 1.0],
                  marker='.', notag=False, ls='--'):
    """draws isotopic masses vs atomic mass,
    returns label and line element for legend

    Args:
        ax(mpl.axes): axes to draw on.
        mdict(dict): dictionary of species and properties.
        label(str): plot label.
        color(str): line color.
        xlims(float list): atomic mass limits (avoids tagging points).
        ylims(float list): mass fraction limits.
        notag(bool): skip species tag on line.

    Returns:
        (mpl.line2D[0], str): line object for legend, tag of plot.

    """
    for k in mdict.keys():
        zz = mdict[k]['z']
        vals = []
        for n in mdict[k]['n'].keys():
            vals.append((n + zz, mdict[k]['n'][n]))
        purgevals = [v for v in vals if ylims[0] <= v[1] <= ylims[1]]
        if len(purgevals) == 0:
            continue
        xs, ys = zip(*sorted(purgevals))
        line = ax.semilogy(xs, ys, ls=ls, lw=0.5,
                           marker=marker, label=label, color=color)
        if notag:
            continue
        elif min(xs) < xlims[0] or max(xs) > xlims[1]:
            continue
        col = plt.gca().lines[-1].get_color()
        ax.text(xs[0] - 0.2, ys[0], '${}$'.format(k), color=color,  # , color=col
                size=8, horizontalalignment='right',
                verticalalignment='bottom')
    # Prettify
    ax.set_xlabel(u'Atomic Mass (A)')
    ax.set_ylabel('Mass ($M_{\odot}$)')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.0e}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def plotAbun(ax, mdict, norm='H', offset=12.0, label='W7', color='black',
             tagspecies=False, yrange=[-3, 15]):
    """draws abundances vs atomic number,
    returns label and line element for legend

    Args:
        ax(mpl.axes): axes to draw on.
        mdict(dict): dictionary of species and properties.
        norm(str): reference element.
        offset(float): abundance offset
        label(str): plot label.
        color(str): line color.
        tagpsecies(bool): plot element names.
        yrange(int list): range of plot for omitting labels.

    Returns:
        (mpl.line2D[0], str): line object for legend, tag of plot.

    """
    zs, names, mix = convertYield2Abundance(mdict, norm=norm, offset=offset)
    line = ax.plot(zs, mix, color=color, label=label, marker='.', ls=':', lw=1)
    # Prettify
    if tagspecies:
        for x, n, y in zip(zs, names, mix):
            if yrange[0] < y < yrange[1]:
                ax.text(x, y, '${}$'.format(n), color='black',
                        size=8, horizontalalignment='right',
                        verticalalignment='bottom')
            else:
                log.warning('skipped tag: {}'.format(n))
                continue
    ax.set_xlabel(u'Atomic Number (Z)')
    ax.set_ylabel('[X/{}] + {:2.1f}'.format(norm, offset))
    # ax.set_xlim([-2, 78])
    # ax.set_ylim(ylims)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def speciesGrid(ax, spcodes, yoffset=5.0):
    """draws a grid of the specified species on target ax.

    Args:
        ax(mpl.axes): axes to draw on.
        spcodes(str list): list of species.
        yoffset(float): line tag offset in plot.

    """
    Sp, Zs, Ns, As = split_species(spcodes)
    for i, sp in enumerate(Sp):
        props = next(ax._get_lines.prop_cycler)
        ax.axvline(Zs[i], alpha=1.0, lw=1, ls='--', color=props['color'])
        ax.text(Zs[i]+0.1, ax.get_ylim()[0]+yoffset,
                '${}$'.format(sp), size=10, horizontalalignment='left',
                verticalalignment='bottom', color=props['color'])


def getNuclideGridNames(names, zs, ns, rightside=False):
    """Retuns tags and their positions for
    a nuclide rectangle grid.

    Args:
        names(str list): species names.
        zs(int list): proton numbers.
        ns(int list): neutron numbers.
        rightside(bool): tags on the right of rectangles.

    Returns:
        (str list, float list, float list):
            unique tags and their plot positions.

    """
    if 'p' in names and 'H' in names:
        names[names.index('p')] = 'H'  # change proton name to void H tag
    if rightside:
        # make a zip to get the max N to each Z
        pairs = zip(zs, ns)
    tags = set()
    nams, x, y = [], [], []
    for i, n in enumerate(names):
        if n not in tags:
            tags.add(n)
            nams.append(n)
            if rightside:
                rightmostN = [n for (z, n) in zip(zs, ns) if z == zs[i]]
                y.append(zs[i]-0.2)
                x.append(rightmostN[-1]+0.7)
            else:
                y.append(zs[i])
                x.append(ns[i]-0.7)
    return nams, x, y
