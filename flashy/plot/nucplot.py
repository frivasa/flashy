from ..nuclear import \
splitSpecies, convertYield2Abundance, getMassFractions,\
getAbundances, readSunComp, readYield, sortNuclides,\
elemSplit, AGSS09, readIsotopicSolar
from .globals import *
from ..datahaul.hdfdirect import directMeta
from matplotlib.patches import Rectangle
import flashy.utils as ut


def speciesYields(yieldfiles, tags=[], norm='Si', offset=6):
    """plots abundances for a list of files or a single file vs solar composition."""
    fig, ax = plt.subplots(figsize=(10,5))
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
                                    color=col, label=tags[i], tagspecies=tagspecies))
            lg = ax.legend(*zip(*plabels))
        # speciesGrid(ax, modspecies, yoffset=5.0)
    else:
        plabels = []
        ryield = readYield(yieldfiles)
        tag = tags
        plabels.append(plotAbun(ax, ryield, norm=norm, offset=offset, label=tag, tagspecies=True))
        # plot sun
        zs, names, mix = readSunComp(AGSS09)
        massF = getMassFractions(names, mix)
        p = ax.plot(zs, getAbundances(names, massF, scale=norm, offset=offset), 
                    marker='x', markersize=3, ls=':', lw=1, color='brown', label='Sun(AGSS09)')
        lg = ax.legend()
    return fig


def productionFactor(yieldfiles, tag='Sim vs X', norm='Si', offset=6):
    """plots profuction factors for comparable species in a pair of yields, or vs solar compositon 
    for a single file.
    yieldfiles = [input, reference]
    """
    fig, ax = plt.subplots(figsize=(10,5))
    if isinstance(yieldfiles, list):
        if len(yieldfiles)!=2:
            print('Input files must be exactly two.')
            return None
        plabels = []
        ryield = readYield(yieldfiles[0])
        plabels.append(plotPfac(ax, ryield, refname=yieldfiles[1], norm=norm, offset=6,
                                reftype='yield', label=tag, tagspecies=True))
        lg = ax.legend()
        # speciesGrid(ax, modspecies, yoffset=5.0)
    else:
        plabels = []
        ryield = readYield(yieldfiles)
        plabels.append(plotPfac(ax, ryield, label='Sim vs AGSS09', 
                                norm=norm, offset=offset, tagspecies=True))
        lg = ax.legend()
    return fig


def nuclideYields(files, tags):
    """plot mass fractions for each nuclide in a yield file."""
    if isinstance(files, list):
        fig, ax = plt.subplots(figsize=(10,5))
        plabels = []
        for i, (f, t) in enumerate(zip(files, tags)):
            ryield = readYield(f)
            props = next(ax._get_lines.prop_cycler)
            col = props['color'] if i else 'k'
            plabels.append(plotIsoMasses(ax, ryield, notag=False, label=t, color=col))
        lg = ax.legend(*zip(*plabels))
        return fig
    else:
        fig, ax = plt.subplots(figsize=(10,5))
        plabels = []
        ryield = readYield(files)
        plabels.append(plotIsoMasses(ax, ryield, notag=False, label=tags, color='k'))
        ryield = readIsotopicSolar()
        plabels.append(plotIsoMasses(ax, ryield, notag=False, label='AGSS09', color='orange'))
        lg = ax.legend(*zip(*plabels))
        return fig


def plotGridYield(yieldfile, dpi=150, cmap='Oranges', 
                  filetag='gridspec', tag='', batch=False):
    """plots mass yields over whole species grid for a yield file.
    
    Args:
        yieldfile(str): path of file.
        dpi(float): dpi of figure.
        cmap(str): colormap name.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        
    Returns:
        (mpl figure) or (None)
    
    """
    names = np.genfromtxt(yieldfile, dtype='|U4', usecols=(0,))
    vals = np.genfromtxt(yieldfile, usecols=(1,))
    mass = np.sum(vals)
    f, a = plt.subplots(dpi=dpi)
    sq = plotNuclideGrid(a, names, xmass=vals, cmap=cmap)
    cax = f.add_axes([0.12, 0.9, 0.8, 0.037])  # left, bot, width, height
    cmap = mpl.cm.Oranges
    # norm = mpl.colors.Normalize(vmin=0.0, vmax=0.5)
    norm = mpl.colors.LogNorm(vmin=1e-3, vmax=0.5)
    cbar1 = mpl.colorbar.ColorbarBase(cax, cmap=cmap, norm=norm, orientation='horizontal')
    cax.xaxis.set_ticks_position('top')
    if tag:
        # trick to make an empty handle 
        sq = plt.Rectangle((0,0), 1, 1, fill=False, edgecolor='none', visible=False)
        leg = a.legend(handles=[sq], labels=[tag], loc=(0.0, 0.85), frameon=False)
    if batch:
        folder = os.path.dirname(yieldfile)
        filename = os.path.basename(yieldfile).split('.')[0]
        paths = (folder, filename)
        writeFig(f, paths, filetag)
    else:
        return f


def plotNuclideGrid(ax, species, xmass=[], time=0.0, z_range=[-2, 35], n_range=[-2, 38], 
                    boxsize=1.0, cmmin=1.0e-5, cmmax=0.9, cmap='Blues', addtags=True, tagsonright=False,
                    noedge=False, frameless=False):
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
        cmap(str): colormap name (https://matplotlib.org/examples/color/colormaps_reference.html)
        addtags(bool): print element names.
        tagsonright(bool): move tags to the rightside of the network.
        noedge(bool): remove nuclide marker edges.
        frameless(bool): remove frames, only show arrows.
    
    Returns:
        (mpl.rectangle): handle for mpl legend.
    
    """
    cmap = mpl.cm.get_cmap(cmap)
    norm = mpl.colors.LogNorm(vmin=cmmin, vmax=cmmax)
    if len(xmass)>0:
        xis = xmass
        a = ax.annotate("{:.5f} s".format(time),
                xy=(0.0, 0.0), xytext=(0.72, 0.18), size=12,
                textcoords='figure fraction', xycoords='figure fraction')
                        #, bbox=dict(boxstyle='round', fc='w', ec='k'))
    else:
        xis = [1e-2]*len(species)
    if noedge:
        clr = 'None'
    else:
        clr = 'black'
    nam, zs, ns, As = splitSpecies(species, standardize=True)
    for (z, n, xi) in zip(zs, ns, xis):
        square = plt.Rectangle((n-0.5*boxsize, z-0.5*boxsize),
                boxsize, boxsize, facecolor=cmap(norm(xi)),
                edgecolor=clr)
        ax.add_patch(square)
        mainsquare = square
    if addtags:
        tags, xs, ys = getNuclideGridNames(list(nam), zs, ns, rightside=tagsonright)
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
        for side in ['bottom','right','top','left']:
            ax.spines[side].set_visible(False)
        ax.set_xticks([]) # labels 
        ax.set_yticks([])
        ax.arrow(0, 12, 0, 10, head_width=0.5, head_length=1, fc='k', ec='k')
        ax.arrow(13, 0, 10, 0, head_width=0.5, head_length=1, fc='k', ec='k')
    
    # numbered axes
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.xaxis.set_minor_formatter(StrMethodFormatter(''))
    return mainsquare


def plotReacNet(ax, sunet, matr_shape, forcedZ=0, step=6, xoffset=1, cmap='Blues'):
    """plots a reaction network based on default files at cdx/Network/netname"""
    with open(matr_shape) as f:
        for i, line in enumerate(f):
            if len(line)!=22:
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
    mat = np.zeros((nsp, nsp))
    mat[x.astype(int)-1, y.astype(int)-1] = z+forcedZ
    ax.imshow(mat, cmap=cmap)
    # set labels
    ax.set_yticks(np.arange(0, nsp, step))
    labels = names[::step]
    ax.set_yticklabels(['$^{{{}}}{{{}}}$'.format(*t) for t in labels], 
                       fontsize=5, ha='right')
    ax.set_xticks(np.arange(xoffset, nsp, step))
    labels = names[xoffset::step]
    ax.set_xticklabels(['$^{{{}}}{{{}}}$'.format(*t) for t in labels], 
                       fontsize=5, va='baseline')

    t = '{} Isotopes\n{} Rates'
    note = ax.annotate(t.format(nsp, rates), xy=(100,30), fontsize=10)
    return nsp, rates


def doYouBelieveInMagic(ax, color='brown'):
    """Draws a magical grid of stability."""
    cadabra = [2, 8, 20, 28, 50, 82, 126]
    for abra in cadabra:
        ax.axhline(abra, alpha=0.4, color=color, ls=':')
        p = ax.axvline(abra, alpha=0.4, color=color, ls=':')
    return p


def plotPfac(ax, querym, refname=AGSS09, label='Sun vs Ref',#  ylims=[1e-9, 1],
             norm='Si', offset=6, reftype='solar', tagspecies=False):
    """draws abundquery/abundref from a massdict and a filename,
    types of reference are 'solar'(for solar composition) and 'yield' (for other sims)
    returns label and line element for legend"""
    zs, ns, ab = convertYield2Abundance(querym, norm=norm, offset=offset)
    if reftype=='solar':
        zss, nss, solab = readSunComp(refname)
        massF = getMassFractions(nss, solab)
        rescaledsolab = getAbundances(nss, massF, scale=norm, offset=offset)
        soldict = dict(zip(nss, rescaledsolab))
    elif reftype=='yield':
        massdict = readYield(refname)
        zss, nss, solab = convertYield2Abundance(massdict, norm=norm, offset=offset)
        soldict = dict(zip(nss, solab))
    pfacs = []
    for i, n in enumerate(ns):
        if n in soldict:
            pfacs.append((zs[i], n, ab[i]-soldict[n]))
            #pfacs.append((zs[i], n, ab[i]/soldict[n]))
        else:
            print('{} not found in ref.'.format(n))
    x, ns, y = zip(*sorted(pfacs))
    ax.axhline(0, ls=':', lw=1, color='green')
    line = ax.plot(x, y, label=label, ls='--', lw=0.5, marker='.', color='black')
    # Prettify
    if tagspecies:
        for x, n, y in sorted(pfacs):
            ax.text(x, y, '${}$'.format(n), color='black',
                    size=8, horizontalalignment='right', verticalalignment='bottom')
    ax.set_xlabel(u'Atomic Number (Z)')
    ax.set_ylabel('$[X/{}]- [X/{}]_{{ref}}$'.format(norm, norm))
    ax.autoscale()
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def plotIsoMasses(ax, mdict, label='W7', color='black', ylims=[1e-18, 1.0], notag=False):
    """draws isotopic masses vs atomic mass, returns label and line element for legend"""
    for k in mdict.keys():
        zz = mdict[k]['z']
        vals = []
        for n in mdict[k]['n'].keys():
            vals.append((n+zz, mdict[k]['n'][n]))
        purgevals = [v for v in vals if ylims[0]<= v[1]<=ylims[1]]
        if len(purgevals)==0:
            continue
        xs, ys = zip(*sorted(purgevals))
        line = ax.semilogy(xs, ys, ls='--', lw=0.5, marker='.', label=label, color=color)
        if notag:
            continue
        ax.text(xs[0], ys[0], '${}$'.format(k), color=color,
                size=8, horizontalalignment='right', verticalalignment='bottom')
    # Prettify
    ax.set_xlabel(u'Atomic Mass (A)')
    ax.set_ylabel('Mass ($M_{\odot}$)')
    ax.set_xlim([-2, 78])
    ax.set_ylim(ylims)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.0e}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def plotAbun(ax, mdict, norm='H', offset=12.0, label='W7', color='black', tagspecies=False):
    """draws abundances vs atomic number, returns label and line element for legend"""
    zs, names, mix = convertYield2Abundance(mdict, norm=norm, offset=offset)
    line = ax.plot(zs, mix, color=color, label=label, marker='.', ls=':', lw=1)
    # Prettify
    if tagspecies:
        for x, n, y in zip(zs, names, mix):
            ax.text(x, y, '${}$'.format(n), color='black',
                    size=8, horizontalalignment='right', verticalalignment='bottom')
    ax.set_xlabel(u'Atomic Number (Z)')
    ax.set_ylabel('[X/{}] + {:2.1f}'.format(norm, offset))
    #ax.set_xlim([-2, 78])
    #ax.set_ylim(ylims)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def speciesGrid(ax, spcodes, yoffset=5.0):
    """draws a grid of the specified species on target ax."""
    Sp, Zs, Ns, As = splitSpecies(spcodes)
    for i, sp in enumerate(Sp):
        props = next(ax._get_lines.prop_cycler)
        ax.axvline(Zs[i], alpha=1.0, lw=1, ls='--', color=props['color'])
        ax.text(Zs[i]+0.1, ax.get_ylim()[0]+yoffset, '${}$'.format(sp), size=10,
                horizontalalignment='left', verticalalignment='bottom', color=props['color'])
    return 0


def getNuclideGridNames(names, zs, ns, rightside=False):
    """Retuns tags and positions for a nuclide grid."""
    if 'p' in names and 'H' in names:
        names[names.index('p')] = 'H'  # change proton name to void H tag
    if rightside:
        # make a zip to get the max N to each Z
        pairs = zip(zs,ns)
    tags = set()
    nams, x, y = [], [], []
    for i, n in enumerate(names):
        if n not in tags:
            tags.add(n)
            nams.append(n)
            if rightside:
                rightmostN = [n for (z,n) in zip(zs, ns) if z==zs[i]]
                y.append(zs[i]-0.2)
                x.append(rightmostN[-1]+0.7)
            else:
                y.append(zs[i])
                x.append(ns[i]-0.7)
    return nams, x, y