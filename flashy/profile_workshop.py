from .utils import np
from copy import deepcopy
from .datahaul.helmholtz import getPres
from .datahaul.plainText import spliceProfs, snipProf
from .wdprofiles.coldFermi import buildFermiHelmhotz
from .wdprofiles.polytropes import buildPolytropicHelmholtz


def goldenIdolShell(dmatr, shellmass):
    """replaces composition to pure helium in a profile down to 'shellmass'"""
    shellc = np.where(dmatr.masses<(dmatr.masses[-1]-shellmass))[0][-1]
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


def polyShell(dmatr, cut, index=1.5):
    """return a polytropic shell for a given wd, starting at 'cut' radius."""
    trimmed = snipProf(dmatr, cut)
    cr, cd, cp, ct, xmass, specs = getCond(trimmed, -1)
    # index = 3.0  # ultra-rel degenerate
    # index = 1.5  # non-rel degenerate
    shell = buildPolytropicHelmholtz(cr, cp, cd, gamma=(index+1)/index, 
                                     species=specs, xmass=xmass)
    print ('Generated shell mass: {:E}'.format(shell.meta['mass']))
    return shell


def getCond(dmatr, where, verbose=True):
    """returns all available conditions in a profile at a given cell.
    for dmatr, 0 is center and -1 is edge (MESA is reversed).
    
    Args:
        dmatr(dataMatrix): profile to modify.
        where(int): cell to probe.
    
    Returns:
        (float): central radius.
        (float): central density.
        (float): central pressure.
        (float list): mass fractions.
        (str list): corresponding species names.
    
    """
    xmasses = [getattr(dmatr, x)[where] for x in dmatr.species]
    cdens = dmatr.density[where]
    cradi = dmatr.radius[where]
    splist = dmatr.species
    cpres = getPres(dmatr.density[where], dmatr.temp[where], xmasses, dmatr.species)[0]
    ctemp = dmatr.temperature[where]
    if verbose:
        print ('\n'.join(['{:5}: {:>16.5E}'.format(s, v) for (s, v) in zip(splist, xmasses)]))
        print ('Radius:{:>16.4E}\nPressure:{:>14.4E}\n'\
               'Density:{:15.4E}\nTemperature:{:11.4E}'.format(cradi, cpres, cdens, ctemp))
        CO = [v for (s, v) in zip(splist, xmasses) if s in ['c12', 'o16']]
        print ('Central C/O: {:>10.4e}'.format(CO[0]/CO[1]))        
    return cradi, cdens, cpres, ctemp, xmasses, splist
    

def getMaximaPositions(dmatr):
    """returns cell-at-maximum for all properties in a profile."""
    pos = []
    keys = dmatr.bulkprops+dmatr.species
    for k in keys:
        pos.append(np.argmax(getattr(dmatr, k)))
    return dict(zip(keys, pos))


def getMaxima(dmatr):
    """returns dictionary with maximal values."""
    maxd = {}
    for (k, v) in getMaximaPositions(dmatr).items():
        maxd[k] = getattr(dmatr, k)[v]
    return maxd




