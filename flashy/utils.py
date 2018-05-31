import numpy as np
from scipy.constants import *
from astropy.constants import M_sun, L_sun, R_sun
# CODATA 2014
# tex source: https://arxiv.org/abs/1507.07956
# web source: https://physics.nist.gov/cuu/Constants/index.html
# enabled through scipy.constants
# WARNING: EVERYTHING IS ADJUSTED TO CGS
G = G*hecto*hecto*hecto/kilo  # Newtonian constant of gravitation cm3/g/s2
c = c*hecto  # speed of light in vacuum cm/s
h = Planck/erg  # Planck constant erg/s
sigma = value('Stefan-Boltzmann constant')/hecto/hecto/erg
kb = value('Boltzmann constant')/erg
amu = value('unified atomic mass unit')*kilo
msol = M_sun.value*kilo
lsol = L_sun.value/erg
rsol = R_sun.value*hecto
m_e = m_e*kilo
m_p = m_p*kilo
Rg = gas_constant

def byMass(radii, dens):
    """Returns a mass abscissa for plots.
    
    Args:
        radii (list of float): cell edge radii.
        dens (list of float): corresponding cell densities.
    
    Returns:
        (list of float): cell mass array.
    
    """
    xs = len(radii)
    dr = radii[0]
    vol = dr**3 *4.0*np.pi/3.0
    mass = vol*dens[0]/msol
    absc = []
    absc.append(mass)
    for i in range(1, xs):
        dr = radii[i] - radii[i-1]
        dvol = dr * ( 3.0*radii[i-1]*radii[i] + dr*dr ) * 4.0*np.pi/3.0
        mass = mass + dvol*dens[i]/msol
        absc.append(mass)
    return absc


def getBearing(angles, geom):
    if len(angles)==0:
        if geom=='cartesian':
            return 'x', np.array([1.0, 0.0, 0.0])
        else:
            return 'r', np.array([1.0, 0.0, 0.0])
    elif len(angles)==1:
        angles = np.radians(angles)
        y = np.sin(angles[0])
        x = np.sqrt(1.0-y**2)
        return 'radius', np.array([x, y, 0.0])
    elif len(angles)==2:
        angles = np.radians(angles)
        x = np.sin(angles[0])*np.cos(angles[1])
        y = np.sin(angles[0])*np.sin(angles[1])
        z = np.sqrt(1.0-y**2-x**2)
        return 'radius', np.array([x, y, z])
    else:
        print('utils.getBearing: too many angles.')
        return -1   
    # spherical (r, theta, phi)
    # cartesian (x, y, z)


def roughCJ(dens, pres, index):
    """returns rayleigh line velocity around index.
    
    Args:
        dens (list of float): density lineout.
        pres (list of float): pressure lineout.
        index (int): position of shock ("shocked cell").
    
    Returns:
        (float): calculated "rayleigh velocity".
    
    """
    dp = float(pres[index-1]-pres[index+1])
    nu1 = 1.0/float(dens[index+1])
    nu2 = 1.0/float(dens[index-1])
    dnu = (nu1-nu2)/(nu1**2)
    # print dp, nu1, nu2, dnu
    if dnu==0.0:
        dnu=1.0
    svel = np.sqrt(dp/dnu)
    return svel


def locateShock(radii, soundcs, xguess, vvv=True):
    """returns index within ray of detected shock.
    
    Args:
        radii (list of float): ordinates.
        soundcs (list of float): sound speed lineout.
        xguess (float): reference position for shock finding.
    
    Returns:
        (int, int): shocked cell indices (inwards and outwards from xguess).
    
    """
    filt, offs1 = split(radii, xguess, True)
    shockin = shock1D(radii[filt], soundcs[filt][:-1], True)
    filt, offs2 = split(radii, xguess, False)
    shockout = shock1D(radii[filt], soundcs[filt][:-1], False)
    if vvv:
        print('Ignition Center: ', xguess)
        print('Inward Shock at: {:E}'.format(float(radii[shockin+offs1])))
        print('Outward Shock at: {:E}'.format(float(radii[shockout+offs2])))
    return shockin+offs1, shockout+offs2


def shock1D(radii, soundspeeds, inward=True):
    """finds a shock in an array by detecting the last 
    large variation within it that is larger than the mean of deltas.
    
    Args:
        radii (list of float): ordinates.
        soundcs (list of float): sound speed lineout.
        inward (bool): return inward or outward shock position.
    
    Returns:
        (int): shocked cell index.
    
    """
    dr = np.diff(radii)[:-1]
    ds = np.diff(soundspeeds)
    div = np.nan_to_num(ds/dr)
    accel = np.abs(np.nan_to_num(ds/dr))
    mean = np.mean(accel)
    pertp = np.where(accel>mean)
    if inward:
        return pertp[0][0]
    else:
        return pertp[0][-1]
    

def split(x, xsplit, inward=True):
    """finds indices below or above xsplit and offset
    [0,1,2,3,4,5,6]
    inward True, xsplit 3: [0,1,2], 0
    inward False, xsplit 3: [4,5,6], 4
    
    Args:
        x (list of float): sorted list of values.
        xsplit (float): pivoting value.
        inward (bool): cut below/above(True/False) pivot.

    Returns:
        (list, float): cut indices, offset.
    
    """
    if inward:
        return np.where(x<xsplit), 0
    else:
        filt = np.where(x>xsplit)
        return filt, filt[0][0]


def estimateMatch(direction, paramd, vvv=True):
    """returns index within ray of detected shock. """
    igr = np.array((paramd['x_match'], paramd['y_match'], paramd['z_match']))
    if np.linalg.norm(igr)==0:  # match is a ring so return the middle of the ring
        spx = paramd['r_match_inner'] + 0.5*(paramd['r_match_outer']-paramd['r_match_inner'])
    else:
        name, r = getBearing(direction, paramd['geometry'])
        cross = np.linalg.norm(np.cross(igr,r))
        dot = np.linalg.norm(np.dot(igr, r))
        angle = np.arctan(cross/dot)
        spx = np.cos(angle)*np.linalg.norm(igr)
    if vvv:
        print('Ignition Center: ', igr)
        print('Estimated position', spx)
    return spx


def cart2sph(x, y, z):
    """returns radius, polar angle and azimuth for a vector."""
    r = np.sqrt(x**2+y**2+z**2)
    if r==0.0:
        return 0, 0, 0
    if x==0.0:
        phi = 90.0
    else:
        phi = np.arctan(y/x)
    tht = np.arccos(z/r)
    return r, thr, phi


def percentDiff(x1, y1, x2, y2):
    """returns the percentage difference between two abscissas 
    subject to the x range of the first via interpolation.
    
    Args:
        x1(float list): reference abscissa.
        y1(float list): comparison ordinate.
        x2(float list): interpolant ordinate.
        y2(float list): interpolant abscissa.
    
    Returns:
        (float list): percentage difference (vs y1, i.e., <0 implies y1<y2)
    
    """
    jake = np.interp(x1, x2, y2)
    diffs = np.array([j/y1-1.0 for (y1, j) in zip(y1, jake)])
    return 100*diffs


def x2clog(x, cmin=1e-5, cmax=1.0):
    """normalizes log10(input) to (cmin, cmax) range."""
    y = (np.log10(x)-np.log10(cmin))/(np.log10(cmax)-np.log10(cmin)) 
    if y > 1.0:
        return 1.0
    elif y < 0:
        return 0.0
    else:
        return y
