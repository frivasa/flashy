"""CGS and Nuclide module. This module can handle abundances,
mass fractions and Isotope weighing.

"""
from .utils import np, msol, Avogadro, m_e, m_p, amu
import pkg_resources

AMDC = pkg_resources.resource_filename('flashy', '../data/nuclides/ame2017.masses')
CIAAW = pkg_resources.resource_filename('flashy', '../data/nuclides/ciaaw2013.masses')
AGSS09 = pkg_resources.resource_filename('flashy', '../data/suncomp/AGSS09.comp')
solDecays = pkg_resources.resource_filename('flashy', '../data/yields/decay_paths.dat')

def readStandardWeights():
    """Returns a species dict with standard atomic weights from:
    http://www.ciaaw.org/atomic-weights.htm:
    http://www.ciaaw.org/pubs/TSAW2013_xls.xls
    # CIAAW only tabulates 2013 values  so no Tennessine :( and only in xls.
    # Abridged values. Ranges are picked on their lower end.
    """
    raw = np.genfromtxt(CIAAW, dtype='i4,|U4,f8')
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
        for _ in range(39):
            next(f)
        for line in f:
            l = line.replace('#', '.')
            n, z, name = int(l[6:10]), int(l[11:15]), l[20:23]
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
    # AMDC doesn't list m_e or m_p so add them by hand from scipy constants
    massdict['p'] = {}
    massdict['p']['n'] = {}
    massdict['p']['n'][0] = m_p/amu
    massdict['p']['z'] = 1
    massdict['e'] = {}
    massdict['e']['n'] = {}
    massdict['e']['n'][0] = m_e/amu
    massdict['e']['z'] = -1
    return massdict


def readDecays():
    sunData = np.genfromtxt(solDecays, dtype="U5,i4,i4,f8,i4")
    snames, _, _, _, scodes = zip(*sunData)
    _, solZs, _, solAs = splitSpecies([s.capitalize() for s in snames], trueA=False)
    return snames, solZs, solAs, scodes


def decayYield(names, masses):
    snames, solZ, solA, scodes = readDecays()
    sunsp = len(snames)
    _, Zs, Ns, As = splitSpecies([s.capitalize() for s in names], standardize=True, trueA=False)
    tmass = sum(masses)
    normedm = masses/tmass
    m_out = np.zeros(sunsp)
    for (name, zin, ain, nmass) in zip(names, Zs, As, normedm):
        for i in range(sunsp):
            if ain==solA[i] and \
            ((scodes[i] == 0) or \
             (scodes[i] == 1 and zin>=solZ[i]) or \
             (scodes[i] == 2 and zin<=solZ[i]) or \
             (scodes[i] == 3 and zin==solZ[i])):
                m_out[i] = m_out[i] + nmass
    return snames, m_out*tmass
    

def readYield(filename):
    """Reads a mass yield file or a zipped list (for simulation yields)
    Returns a mass dictionary."""
    if type(filename)==list:
        codes, values = zip(*filename)
    else:
        comp = np.genfromtxt(filename, comments='#', dtype='|U5,f8')
        codes, values = zip(*comp)
    sp, z, n, a = splitSpecies(codes, standardize=True)
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
            particles += mdict[k]['n'][n]*msol/nucmass[k]['n'][n]*Avogadro
        partdict[k] = particles
    #print partdict
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
    comp = np.genfromtxt(filename, comments='#', dtype=None, encoding=None)
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
        print('Scaling species not in the specified species group.')
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
    # atot = reduce(lambda x, y: x + y, masses)
    atot = sum(masses)
    # z = reduce(lambda x, y: x + y, masses[2:])/atot
    z = sum(masses[2:])/atot
    x = masses[0]/atot
    y = masses[1]/atot
    return x, y, z


def getTotalMass(massdict):
    """Sums all isotope masses in a massdict"""
    mass = 0.0
    for k in massdict.keys():
        for n in massdict[k]['n']:
            mass += massdict[k]['n'][n]
    print(sum(mass))


def convXmass2Abun(species, xmasses):
    """Returns abundances, abar and zbar from a list of nuclide
    codes and mass fractions.
    
    Args:
        species(list of str): nuclide code list.
        xmasses(list of float): mass fractions of each species.
    
    Returns:
        (list of float): MOLAR abundances for species.
        (float): inverse of molar sum (abar).
        (float): average charge (zbar).
    
    """
    _, Zs, _, As = splitSpecies(species, trueA=True)
    molar = [x/a for (x,a) in zip(xmasses, As)]
    abar = 1.0e0/sum(molar)
    zbar = abar * sum([x*z/a for (z,x,a) in zip(Zs, xmasses, As)])
    return molar, abar, zbar


def getMus(species, xmasses):
    """Returns abundances, abar and zbar from a list of nuclide
    codes and mass fractions.
    
    Args:
        species(list of str): nuclide code list.
        xmasses(list of float): mass fractions of each species.

    Returns:
        (float): molec weight per free particle (ions + e).
        (float): molec weight per free electron.
    
    """
    _, Zs, _, As = splitSpecies(species, trueA=True)
    mue = 1.0e0/sum([x*z/a for (z,x,a) in zip(Zs, xmasses, As)])
    muion = 1.0/sum([x*(z+1)/a for (z,x,a) in zip(Zs, xmasses, As)])
    return muion, mue


def splitSpecies(Spcodes, trueA=True, standardize=False):
    """returns list of symbols, Z, N, A from a list of
    nuclide codes.(Ni56, He4, U238, ...)
    
    Args:
        Spcodes(str list): nuclide code list (Ni56, He4, U238, ...).
        trueA(bool): return real weight instead of N+Z.
        standardize(bool): return element for special names (deuteron, proton, tritium).
    
    """
    Sp, As = zip(*[elemSplit(s) for s in Spcodes])
    As = np.array(As)
    mdict = readNuclideMasses()
    if standardize:
        # to allow the Nitrogen/neutron disctinction, both upper and lower cases must be checked
        # proton/Phosphorous
        stdnames = [('p', 'H'), ('d', 'H'), ('D', 'H'), ('t', 'H'), ('T', 'H'), ('h', 'H')]
        Sp = list(Sp)
        for (k, v) in stdnames:
            if k in Sp:
                indices = [i for i, s in enumerate(Sp) if s==k]  # find all indices for repeated species
                for j in indices:
                    Sp[j] = v
    Zs = np.array([mdict[n]['z'] for n in Sp])
    Ns = As - Zs
    if trueA:
        As = [mdict[s]['n'][n] for (s,n) in zip(Sp,Ns)]
    return Sp, Zs, Ns, As


def sortNuclides(spcodes, capitalize=False):
    """sorts a list of nuclides by atomic number."""
    tuples = [elemSplit(s, True) for s in spcodes]
    nucs = ['{}{}'.format(s, a) for (a, s) in sorted(tuples)]
    if capitalize:
        return nucs
    else:
        return [n.lower() for n in nucs]


def elemSplit(s, invert=False):
    """Standalone element name spliter.
    (A, name)
    he4 -> (He, 4)
    """
    if len(s.strip())==1:
        A = {
        'n':1,
        'p':1,
        'd':2,
        't':3
        }[s.strip().lower()]
        sym = s.strip().lower()
    else:
        sym = s.strip('0123456789 ')
        A = s[len(sym):]
        if int(A)>2:
            sym = sym.capitalize()
    if invert:
        return int(A), sym
    else:
        return sym, int(A)
