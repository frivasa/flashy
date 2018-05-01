"""get mass yields from a checkpoint. (similar to cjvelocities module)
"""
import flashy.datahaul.hdf5yt as reader
from flashy.nuclear import convXmass2Abun
from scipy.integrate import trapz
import flashy.utils as ut

def getYields(fname, geom='spherical'):
    """returns time, names and masses for each species in a checkpoint file.
    Note: there's a small difference between sum(masses) and masses[-1] 
    (less than 1e-3 percent)
    
    Args:
        fname(str): filename to inspect.
        geom(str): geometry spec for lineout.
    
    Returns:
        time (float), 
        species names: (list of str),
        masses: (list of float)

    """
    data, species = reader.getLineout(fname, fields=['density'], species=True, geom='spherical')
    time, _, _, _, _ = reader.getMeta(fname)
    masses = ut.byMass(data[0], data[1])
    spms = []
    for i in range(len(species)):
        spdat = data[i+len(data)-len(species)]
        spmass = trapz(data[i+len(data)-len(species)], x=masses)
        spms.append(spmass)
    return time, species, spms