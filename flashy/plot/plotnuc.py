from flashy.nuclear import *

def plotPfac(ax, querym, refname=AGSS09, label='Sun vs Ref',#  ylims=[1e-9, 1],
             norm='Si', offset=6, reftype='solar'):
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
            print '{} not found in ref.'.format(n)
    x, _, y = zip(*sorted(pfacs))
    ax.axhline(0, ls=':', lw=1, color='green')
    line = ax.plot(x, y, label=label, ls='--', lw=0.5, marker='.', color='black')
    # Prettify
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
                size=8, horizontalalignment='right', verticalalignment='bottom',)
    # Prettify
    ax.set_xlabel(u'Atomic Mass (A)')
    ax.set_ylabel('Mass ($M_{\odot}$)')
    ax.set_xlim([-2, 78])
    ax.set_ylim(ylims)
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:0.0e}'))
    ax.xaxis.set_major_formatter(StrMethodFormatter('{x:2.0f}'))
    return line[0], label


def plotAbun(ax, mdict, norm='H', offset=12.0, label='W7', color='black'):
    """draws abundances vs atomic number, returns label and line element for legend"""
    zs, names, mix = convertYield2Abundance(mdict, norm=norm, offset=offset)
    line = ax.plot(zs, mix, color=color, label=label, marker='.', ls=':', lw=1)
    # Prettify
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
    clrs = colIter()
    for i, sp in enumerate(Sp):
        col = clrs.next()[0]
        ax.axvline(Zs[i], alpha=1.0, lw=1, ls='--', c=col)
        ax.text(Zs[i]+0.1, ax.get_ylim()[0]+yoffset, '${}$'.format(sp), color=col, size=10,
                horizontalalignment='left', verticalalignment='bottom')
    return 0