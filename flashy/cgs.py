"""CGS and Nuclide module. This module can handle abundances,
mass fractions and Isotope weighing.

Todo:
    * read both compositions and mass fractions and relate them through atomic masses.
    * plotting functions.

"""
import numpy as np
import os
from scipy.constants import *
from astropy.constants import M_sun, L_sun, R_sun
import pkg_resources
from matplotlib.ticker import StrMethodFormatter
from flashy.plot._globals import colIter

AMDC = pkg_resources.resource_filename('flashy', 'data/nuclides/ame2017.masses')
CIAAW = pkg_resources.resource_filename('flashy', 'data/nuclides/ciaaw2013.masses')
AGSS09 = pkg_resources.resource_filename('flashy', 'data/suncomp/AGSS09.comp')

# CODATA 2014
# tex source: https://arxiv.org/abs/1507.07956
# web source: https://physics.nist.gov/cuu/Constants/index.html
# enabled through scipy.constants
# WARNING: EVERYTHING IS ADJUSTED TO CGS
G = G*hecto*hecto*hecto/kilo  # Newtonian constant of gravitation cm3/g/s2
c = c*hecto  # speed of light in vacuum cm/s
h = Planck/erg  # Planck constant erg/s 
sigma = value('Stefan-Boltzmann constant')/hecto/hecto/erg
msol = M_sun*kilo
lsol = L_sun/erg
rsol = R_sun*hecto


def readStandardWeights():
    """Returns a species dict with standard atomic weights from:
    http://www.ciaaw.org/atomic-weights.htm:
    http://www.ciaaw.org/pubs/TSAW2013_xls.xls
    # CIAAW only tabulates 2013 values  so no Tennessine :( and only in xls.
    # Abridged values. Ranges are picked on their lower end.
    """
    raw = np.genfromtxt(CIAAW, dtype='i4,|S4,f8')
    dblock = [x for x in raw if x[2]>0.0]
    keys = [a[1] for a in dblock] + [a[0] for a in dblock]  # set both name and Z as keys
    vals = [a[2]for a in dblock]
    masses = dict(zip(keys, 2*vals))
    return masses
    

def readNuclideMasses():
    """Reads an AMDC format table (Wang 2017 41 030003).
    AMDC atomic masses (2017): http://amdc.in2p3.fr/web/masseval.html
    Reads only masses, there's more data within it so 
    this function is extensible.
    output is a dict with capitalized element names, then a 
    z and n key to access z and a keyset of N numbers.
    """
    massdict = {}
    with open(AMDC, 'r') as f:
        # skip header
        for _ in xrange(39):
            next(f)
        for line in f:
            l = line.replace('#', '.')
            n, z, name = int(l[6:10]), int(l[12:15]), l[20:23]
            excess, exunc = float(l[30:42]), float(l[43:54])  # [M(in u) - A] keV
            bind, biunc = float(l[55:64]), float(l[66:73])  # B/A in keV
            try:
                beta, beunc = float(l[77:87]), float(l[88:97])
            except:
                beta, beunc = np.nan, np.nan
            mass, maunc = float(l[96:114].replace(' ','')), float(l[115:])
            # neutrons have the same name as nitrogen x_x 
            # dictionaries are case sensitive though.
            if ' n' in name:
                name = 'n'
            else:
                name = name.strip().capitalize()
            if name in massdict:
                massdict[name]['n'][n] = mass*1e-6  # source mass is in micro u
            else:
                massdict[name] = {}
                massdict[name]['n'] = {}
                massdict[name]['n'][n] = mass*1e-6  # source mass is in micro u
                massdict[name]['z'] = z
    return massdict


def readYield(filename):
    """Reads a mass yield file or a zipped list (for simulation yields) 
    Returns a mass dictionary."""
    if type(filename)==list:
        codes, values = zip(*filename)
    else:
        comp = np.genfromtxt(filename, comments='#', dtype='|S4,f8')
        codes, values = zip(*comp)
    sp, z, n, a = splitSpecies(codes)
    mdict = {}
    for i, s in enumerate(sp):
        if s in mdict:
            mdict[s]['z'] = z[i]
            mdict[s]['n'][n[i]] = values[i]
        else:
            mdict[s] = {}
            mdict[s]['n'] = {}
            mdict[s]['n'][n[i]] = values[i]
            mdict[s]['z'] = z[i]
    return mdict


def convertYield2Abundance(mdict, norm='H', offset=12.0):
    """Turns mass yields into abundances."""
    nucmass = readNuclideMasses()
    partdict = {}
    for k in mdict.keys():
        zz = mdict[k]['z']
        particles = 0
        for n in mdict[k]['n']:
            particles += mdict[k]['n'][n]*msol.value/nucmass[k]['n'][n]*Avogadro
        partdict[k] = particles
    print partdict
    nrm = partdict[norm]
    otp = []
    for k, v in partdict.items():
        if v==0:
            continue
        else:
            otp.append((mdict[k]['z'], k, np.log(partdict[k]/nrm)+offset))
    zs, names, mix = zip(*sorted(otp))
    return zs, names, mix


def readSunComp(filename):
    """Read a 6 column file with solar composition
    (Z name abundance abunerror meteoriticAb mAberror)
    mixes meteor and spectral/inferred abundances for completion.
    """
    comp = np.genfromtxt(filename, comments='#', dtype=None)
    zs, names, abun, aerr, met, emet = zip(*comp)
    # pick meteorite values when photospheric is unavailable (AGSS09)
    mix = []
    for x,y in zip(abun,met):
        if x==-9.99:
            mix.append(y)
        else:
            mix.append(x)
    return zs, [n.capitalize() for n in names], mix


def getAbundances(names, massfrac, scale='H', offset=12):
    """Returns abundances for a list of species and mass Fractions, 
    scaled by a set species.
    """
    if scale not in names:
        print 'Scaling species not in the specified species group.'
        return [np.nan]
    else:
        totmass = np.sum(massfrac)
        wgtdict = readStandardWeights()
        wgts = [wgtdict[x] for x in names]
        normnum = names.index(scale)
        allabun = [totmass*massfrac[i]/wgts[i] for i in range(len(names))]
        norm = allabun[normnum]
        allabun = [x/norm for x in allabun]
        return np.log10(allabun)+offset


def getMassFractions(names, abun):
    """Takes log abundances and species names to weigh them and return 
    mass fractions.
    returns mass fractions."""
    wgtdict = readStandardWeights()
    wgts = [wgtdict[x] for x in names]
    powabun = np.power(10, abun)
    masses = powabun*wgts
    totmass = sum(masses)
    return masses/totmass


def getXYZ(masses):
    """Returns Hydrogen, Helium and Metal Fractions from 
    given mass fractions (assumes H and He at start of list).
    """
    atot = reduce(lambda x, y: x + y, masses)
    z = reduce(lambda x, y: x + y, masses[2:])/atot
    x = masses[0]/atot
    y = masses[1]/atot
    return x, y, z


def getTotalMass(massdict):
    """Sums all isotope masses in a massdict"""
    mass = 0.0
    for k in massdict.keys():
        for n in massdict[k]['n']:
            mass += massdict[k]['n'][n]
    print sum(mass)

    
def fetchNuclides(flist, lowercase=True):
    """filters flash otp field list, extracting species found 
    in the checkpoint.
    """
    species = []
    for (t, field) in flist:
        if any(char.isdigit() for char in field):
            species.append(field)
    ss, zs, ns, ms = splitSpecies(species, trueA=False)
    otp = sorted(zip(ms, ss))
    return ['{}{}'.format(s.lower(), m) for (m, s) in otp]


def convXmass2Abun(species, xmasses):
    """Returns abundances, abar and zbar from a list of nuclide 
    codes and mass fractions."""
    _, Zs, _, As = splitSpecies(species, trueA=True)
    abar = 1.0e0/sum([x/a for (x,a) in zip(xmasses, As)])
    zbar = abar * sum([x*z/a for (z,x,a) in zip(Zs, xmasses, As)])
    return [x/a for (x,a) in zip(xmasses, As)], abar, zbar


def splitSpecies(Spcodes, trueA=True):
    """returns list of symbols, Z, N, A from a list of 
    nuclide codes.(Ni56, He4, U238, ...)
    """
    Sp, As = zip(*[elemSplit(s) for s in Spcodes])
    As = np.array(As)
    mdict = readNuclideMasses()
    Zs = np.array([mdict[n]['z'] for n in Sp])
    Ns = As - Zs
    if trueA:
        As = [mdict[s]['n'][n] for (s,n) in zip(Sp,Ns)] 
    return Sp, Zs, Ns, As


def elemSplit(s):
    """Standalone element name spliter.
    (A, name)
    he4 -> (4, He)
    """
    if len(s)==1:
        if s.lower()=='n':
            return 'n', 1
        elif s.lower()=='p':
            return 'H', 1
        elif s.lower()=='d':
            return 'H', 2
        elif s.lower()=='t':
            return 'H', 3
    else:
        sym = s.strip('0123456789 ')
        A = s[len(sym):]
    return sym.capitalize(), int(A)

#Plotting. this should be elsewhere... maybe this whole module can be isolated
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