"""get mass yields  and cj velocities from a checkpoint."""
import flashy.datahaul.hdf5yt as reader
from scipy.integrate import trapz
import flashy.utils as ut
from flashy.utils import h, m_e, np, c, kb, Avogadro


def getYields(fname, geom='spherical', direction=[]):
    """returns time, names and masses for each species in a checkpoint file.
    Note: there's a small difference between sum(masses) and masses[-1] 
    (less than 1e-3 percent)
    direction is unphysical, but leaving it to test on multiD models
    
    Args:
        fname(str): filename to inspect.
        geom(str): geometry spec for lineout.
    
    Returns:
        time (float), 
        species names: (list of str),
        masses: (list of float)

    """
    data, species = reader.getLineout(fname, fields=['density'], species=True, geom=geom, direction=direction)
    time, _, _, _, _ = reader.getMeta(fname)
    masses = ut.byMass(data[0], data[1])
    spms = []
    for i in range(len(species)):
        spdat = data[i+len(data)-len(species)]
        spmass = trapz(data[i+len(data)-len(species)], x=masses)
        spms.append(spmass)
    return time, species, spms


def getVelocities(fname, geom='spherical'):
    """Returns positions of the shock, and both in an outer cj 
    velocities for a file, plus the starting x_match.
    (xin, xout, cjin, cjout, float(ray.ds.current_time), ray.ds.parameters['x_match'])
    """
    fields = ['sound_speed', 'density', 'pressure']
    data, _ = reader.getLineout(fname, fields=fields, species=False, geom=geom)
    time, params, _, _, _ = reader.getMeta(fname)
    rad, cs, dens, pres = data[0], data[1], data[2], data[3]
    shockin, shockout = ut.locateShock(rad, cs, params['x_match'], vvv=False)
    xin, xout = rad[shockin], rad[shockout]
    cjin = ut.roughCJ(dens, pres, shockin)
    cjout = ut.roughCJ(dens, pres, shockout)
    return xin, xout, cjin, cjout, time, params['x_match']


def nonRelFermi(dens, ye=0.5):
    """ Completely degenerate, non-relativistic Fermi energy.
    E_f = (hbar^2/(2m_e))(3pi^(2/3))(N_a \rho Y_e)^(2/3)
    """
    par1 = 0.5*np.power(0.5*h/np.pi, 2.0)/m_e
    par2 = np.power(3.0*np.pi*np.pi,2.0/3)
    par3 = np.power(Avogadro*dens*ye, 2.0/3)
    return par1*par2*par3