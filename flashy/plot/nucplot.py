from flashy.nuclear import \
splitSpecies, convertYield2Abundance, getMassFractions,\
getAbundances, readSunComp, readYield, sortNuclides, elemSplit, AGSS09
from .globals import *
import flashy.utils as ut
from flashy.datahaul.hdf5yt import getLineout, getYields
from flashy.datahaul.hdfdirect import directMeta
from matplotlib.patches import Rectangle


# class patchN(object):
#     def legend_artist(self, legend, orig_handle, fontsize, handlebox):
#         x0, y0 = handlebox.xdescent, handlebox.ydescent
#         width, height = handlebox.width, handlebox.height
#         patch = Rectangle((x0,y0), height, height)
#         patch.update_from(orig_handle)
#         patch.set_transform(handlebox.get_transform())
#         handlebox.add_artist(patch)
#         return patch

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
        plabels.append(plotPfac(ax, ryield, refname=yieldfiles[1], 
                                reftype='yield', label=tag, tagspecies=True))
        lg = ax.legend()
        # speciesGrid(ax, modspecies, yoffset=5.0)
    else:
        plabels = []
        ryield = readYield(yieldfiles)
        plabels.append(plotPfac(ax, ryield, label='Sim vs AGSS09', tagspecies=True))
        lg = ax.legend()
    return fig


def nuclideYields(files, tags):
    """plot mass fractions for each nuclide in a yield file."""
    fig, ax = plt.subplots(figsize=(10,5))
    plabels = []
    for i, (f, t) in enumerate(zip(files, tags)):
        ryield = readYield(f)
        props = next(ax._get_lines.prop_cycler)
        col = props['color'] if i else 'k'
        plabels.append(plotIsoMasses(ax, ryield, notag=False, label=t, color=col))
    lg = ax.legend(*zip(*plabels))
    return fig


def plotNuclideGrid(ax, species, xmass=[], z_range=[-0.5,35], n_range=[-0.5,38], 
                    boxsize=1.0, cmmin=1.0e-5, cmmax=1.0, cmap='Blues', 
                    log=False, addtags=True, invert=False):
    """plots a list of species on a grid."""
    cmap = mpl.cm.get_cmap(cmap)
    if len(xmass)>0:
        xis = xmass
    else:
        fakeXi = 1.0-np.random.rand(1)
        xis = [fakeXi[0]]*len(species)
    if invert:
        clr = 'white'
    else:
        clr = 'black'
    nam, zs, ns, As = splitSpecies(species, standardize=True)
    for (z, n, xi) in zip(zs, ns, xis):
        if log:
            square = plt.Rectangle((n-0.5*boxsize, z-0.5*boxsize),
                    boxsize, boxsize, facecolor=cmap(ut.x2clog(xi)),
                    edgecolor=clr)
        else:
            square = plt.Rectangle((n-0.5*boxsize, z-0.5*boxsize),
                    boxsize, boxsize, facecolor=cmap(xi),
                    edgecolor=clr)
        ax.add_patch(square)
        mainsquare = square
    if addtags:
        tags, xs, ys = getNuclideGridNames(list(nam), zs, ns)
        for (t, x, y) in zip(tags, xs, ys):
            ax.text(x, y, t, fontsize=8, verticalalignment='center', 
                    horizontalalignment='right', color=clr)
    #ax.set_ylabel('Z ($p^+$)', color=clr)
    #ax.set_xlabel('N ($n^0$)', color=clr)
    ax.set_ylabel('Proton Number', color=clr)
    ax.set_xlabel('Neutron Number', color=clr)
    ax.set_ylim(z_range)
    ax.set_xlim(n_range)
    #ax.xaxis.tick_top()
    #ax.xaxis.set_label_position('top')
    # remove spines/ticks
    for side in ['bottom','right','top','left']:
        ax.spines[side].set_visible(False)
    ax.set_xticks([]) # labels 
    ax.set_yticks([])
    ax.arrow(0, 12, 0, 10, head_width=0.5, head_length=1, fc=clr, ec=clr)
    ax.arrow(13, 0, 10, 0, head_width=0.5, head_length=1, fc=clr, ec=clr)
    #ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    #ax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    #ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    #ax.xaxis.set_minor_formatter(StrMethodFormatter(''))
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
    """Draws a magical grid of stability (?)"""
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


def getNuclideGridNames(names, zs, ns):
    """Retuns tags and positions for a nuclide grid."""
    if 'p' in names and 'H' in names:
        names[names.index('p')] = 'H'  # change proton name to void H tag
    tags = set()
    nams, x, y = [], [], []
    for i, n in enumerate(names):
        if n not in tags:
            tags.add(n)
            nams.append(n)
            x.append(ns[i]-0.7)
            y.append(zs[i])
    return nams, x, y