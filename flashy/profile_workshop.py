from .utils import np
from copy import deepcopy
from .datahaul.helmholtz import getPres
from .datahaul.plainText import spliceProfs, snipProf
from .wdprofiles.coldFermi import buildFermiHelmhotz
from .wdprofiles.polytropes import buildPolytropic


def goldenIdolShell(dmatr, limit, reference='rho'):
    """replaces composition in a profile to pure helium from outer edge to
    limit. This can be a density or an amount of mass to turn to helium based
    of the reference key.

    Args:
        dmatr(dataMatrix): starting profile.
        limit(float): cutoff for replacing composition.
        reference(str): cutoff reference ('rho' or 'mass').

    Returns:
        (dataMatrix): shell profile obj.

    """
    if reference == 'mass':
        shellc = np.where(dmatr.masses < (dmatr.masses[-1]-limit))[0][-1]
    elif reference == 'rho':
        shellc = np.where(dmatr.density < limit)[0][0]
    else:
        print('profile_workshop.goldenIdolShell: '
              'invalid reference, use "rho" or "mass"')
        return -1
    shell = snipProf(dmatr, dmatr.radius[shellc], left=False)
    npnts = len(shell.radius)
    # wipe other species
    for s in shell.species:
        setattr(shell, s, np.zeros(npnts))
    # pour the helium
    try:
        setattr(shell, 'he4', np.ones(npnts))
    except AtrributeError:
        shell.species.append('he4')
        setattr(shell, 'he4', np.ones(npnts))
    return shell


def smudgeShell(dmatr, csp='he4'):
    """extends xmasses at the maximum mass fraction of a determined
    species throughout the profile, conserving bulk properties.

    Args:
        dmatr(dataMatrix): profile to modify.
        csp(str): reference species.

    Returns:
        (dataMatrix): modified profile.

    """
    for k, pos in getMaximaPositions(dmatr).items():
        if k == csp:
            cutPos = pos
            cutRadius = dmatr.radius[pos]

    cr, cd, cp, ct, xmass, specs = getCond(dmatr, cutPos, verbose=False)
    shellc = np.where(dmatr.radius > cutRadius)[0]

    for sp, val in zip(specs, xmass):
        array = getattr(dmatr, sp)
        array[shellc] = val
        setattr(dmatr, sp, array)

    return dmatr


def fifty50Shell(dmatr, eps=1e-3, limdens=1e2):
    """build a helium shell starting from the first point where C/O = 0.5,
    keeping bulk properties fixed out to the "real" helium/hydrogen shell.
    This is done to ensure a high enough density for detonation and a
    smooth Abar transition at the core/shell interface
    (other species may skew the 5050 but they're 1e-1 smaller than C or O).
    Physically, this should happen through accretion compressing helium
    over the rho = 1e5 needed and ensuring a continuous transition in Abar
    from any surface C/O ratio.

    Args:
        dmatr(dataMatrix): starting profile.
        eps(float): tolerance for C/0 = 1 ratio.
        limdens(float): density for shell cutoff.

    Returns:
        (dataMatrix): shell profile obj.

    """
    # find the 5050 pos starting from the middle outwards
    for pos in range(int(0.5*len(dmatr.radius))):
        ratio = dmatr.o16[pos]/dmatr.c12[pos]
        if abs(ratio-1.0) < eps:
            break
    coreend = pos
    cutRadius = dmatr.radius[pos]
    envelope = snipProf(dmatr, cutRadius, left=False)
    # wipe other species and pour the helium
    npnts = len(envelope.radius)
    for s in envelope.species:
        setattr(envelope, s, np.zeros(npnts))
    setattr(envelope, 'he4', np.ones(npnts))
    # now trim the shell at a defined density
    for pos in range(len(envelope.radius)):
        if envelope.density[pos] < limdens:
            break
    shellend = pos
    shell = snipProf(envelope, envelope.radius[pos], left=True)
    return shell, coreend, shellend


def polyShell(dmatr, cut, index=1.5, species=['he4'], xmass=[1.0]):
    """return a polytropic shell for a given wd, starting at 'cut' radius.

    Args:
        dmatr(dataMatrix): reference profile.
        cut(float): shell start radius .

    Returns:
        (dataMatrix): shell profile obj.

    """
    trimmed = snipProf(dmatr, cut)
    print("Conditions at cut:")
    cr, cd, cp, ct, _, _ = getCond(trimmed, -1)
    print("Ignoring cut composition...")
    # index = 3.0  # ultra-rel degenerate
    # index = 1.5  # non-rel degenerate
    shell = buildPolytropic(cr, cp, cd, gamma=(index+1)/index,
                            species=species, xmass=xmass)
    print('Generated shell mass: {:E}'.format(shell.meta['mass']))
    return shell


def getCond(dmatr, where, verbose=True):
    """returns all available conditions in a profile at a given cell.
    for dmatr, 0 is center and -1 is edge (MESA is reversed).

    Args:
        dmatr(dataMatrix): profile to modify.
        where(int): cell to probe.

    Returns:
        (float): radius.
        (float): density.
        (float): pressure.
        (float list): mass fractions.
        (str list): corresponding species names.

    """
    xmasses = [getattr(dmatr, x)[where] for x in dmatr.species]
    cdens = dmatr.density[where]
    cradi = dmatr.radius[where]
    splist = dmatr.species
    cpres = getPres(dmatr.density[where], dmatr.temp[where],
                    xmasses, dmatr.species)[0]
    ctemp = dmatr.temperature[where]
    if verbose:
        print('\n'.join(['{:5}: {:>16.5E}'.format(s, v)
                         for (s, v) in zip(splist, xmasses)]))
        print('Radius:{:>16.4E}\nPressure:{:>14.4E}\n'
              'Density:{:15.4E}\nTemperature:{:11.4E}'.format(cradi, cpres,
                                                              cdens, ctemp))
        CO = [v for (s, v) in zip(splist, xmasses) if s in ['c12', 'o16']]
        print('Central C/O: {:>10.4e}'.format(CO[0]/CO[1]))
    return cradi, cdens, cpres, ctemp, xmasses, splist


def getSummedMasses(dmatr):
    """Returns summed masses from all species in the profile.

    Args:
        dmatr(dataMatrix): reference profile.

    Returns:
        (dict): {species: summed masses}.

    """
    yields = []
    keys = dmatr.species
    cell_masses = np.diff(dmatr.masses)
    cell_masses = np.insert(cell_masses, 0, dmatr.masses[0])
    for k in keys:
        # check for nans
        line = np.nan_to_num(getattr(dmatr, k), copy=True)
        yields.append(sum(line*cell_masses))
    # print('Total Summed Mass: {:e}'.format(sum(yields)))
    # print('Last Mass Cell: {:e}'.format(dmatr.masses[-1]))
    return dict(zip(keys, yields))


def getMaximaPositions(dmatr):
    """returns cell-at-maximum for all properties in a profile.

    Args:
        dmatr(dataMatrix): reference profile.

    Returns:
        (dict): {property/species: pos of maximum}.

    """
    pos = []
    keys = dmatr.bulkprops+dmatr.species
    for k in keys:
        # check for nans
        line = np.nan_to_num(getattr(dmatr, k), copy=True)
        pos.append(np.argmax(line))
        # pos.append(np.argmax(getattr(dmatr, k)))
    return dict(zip(keys, pos))


def getInterfacePosition(dmatr, species='he4'):
    """returns interface array position for a species based on
    gradient > mean value.

    Args:
        dmatr(dataMatrix): reference profile.
        species(str): interface nuclide

    Returns:
        (int): position in the DataMatrix for the found interface.
              (0 if species is not found)

    """
    try:
        spec_data = getattr(dmatr, species)
    except AttributeError:
        return 0
    grad = np.gradient(spec_data)
    mean = np.mean(grad)
    deltapos = np.where(grad > mean)
    return deltapos[0][-1]


def getMaxima(dmatr):
    """returns dictionary with maximal values.

    Args:
        dmatr(dataMatrix): reference profile.

    Returns:
        (dict): {property/species: maximum value}.

    """
    maxd = {}
    for (k, v) in getMaximaPositions(dmatr).items():
        maxd[k] = getattr(dmatr, k)[v]
    return maxd
