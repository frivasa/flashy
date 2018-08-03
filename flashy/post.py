"""calculates cj velocities for a checkpoint."""
from flashy.datahaul.hdf5yt import getLineout
from flashy.datahaul.hdfdirect import directMeta
import flashy.utils as ut
from .utils import h, m_e, np, c, kb, Avogadro
from scipy.optimize import newton, minimize


def getFittedVelocities(fname, **kwargs):
    """Analyze a filename, extracting shock position, time, and calculating 
    cj speed in both ends of the shock.
    
    Args:
        fname(str): filepath.
        **kwargs: arguments for lineout(geometry, direction, etc.).
    
    Returns:
        (list of float): xin, cjin, xout, cjout, time, matchhead position.
    
    """
    _, _, cjin, pmi, pm, time = getNewtonCJ(fname, inward=True, **kwargs)
    _, _, cjout, pmo, _, _ = getNewtonCJ(fname, inward=False, **kwargs)
    return pmi, cjin, pmo, cjout, pm, time


def getRayleighVelocities(fname, direction=[]):
    """Returns positions of the shock, and both inner an outer rayleigh line velocities joining 
    both sides of the shocked cell.
    (xin, xout, cjin, cjout, float(ray.ds.current_time), ray.ds.parameters['x_match'])
    """
#     fields = ['sound_speed', 'density', 'pressure']
#     data, _ = reader.getLineout(fname, fields=fields, species=False, geom=geom)
#     time, params, _, _, _ = reader.getMeta(fname)
#     rad, cs, dens, pres = data[0], data[1], data[2], data[3]
    
    fields = ['sound_speed', 'density', 'pressure']
    time, pars, _, _, paths = directMeta(fname)
    if len(direction)>(pars['dimensionality']-1):
        print("Direction doesn't match dimensionality: {}".format(pars['dimensionality']))
        return None
    data, _ = getLineout(fname, fields=fields, species=False, direction=direction, geom=pars['geometry'])
    #time, params, _, _, _ = directMeta(fname)
    rad, cs, dens, pres = data[0], data[1], data[2], data[3]
    
    # this fails for x_match = y_match = z_match = 0.0
    linepos = ut.estimateMatch(direction, pars, vvv=False)
    # print(linepos)
    shockin, shockout = ut.locateShock(rad, cs, linepos, vvv=False)
#     shockin, shockout = ut.locateShock(rad, cs, pars['x_match'], vvv=False)
    
    xin, xout = rad[shockin], rad[shockout]
    cjin = ut.roughCJ(dens, pres, shockin)
    cjout = ut.roughCJ(dens, pres, shockout)
    return xin, xout, cjin, cjout, time, pars['x_match']


def getNewtonCJ(fname, inward=False, width=0.8, **kwargs):
    """Calculate CJ velocity for a file.
    
    Args:
        fname(str): filepath.
        inward(bool): toggle for inward/outward bound shock.
    
    Returns:
        (float tuple): specific volume, pressure, CJVelocity, matchhead position, time.
    
    """
    #print('check newtonCJ')
    pos, dens, pres, gamc, cjest, pm, time = getShockConditions(fname, addvar='gamc', 
                                                                inward=inward, **kwargs)
    # set bulk properties
    fv, fp, fg = 1.0/dens[-1], pres[-1], gamc[-1]
    av, ap, ag = 1.0/dens[0], pres[0], gamc[0]
    try:
        v, p, cj = newtonCJ(cjest, fv, fp, fg, av, ap, ag, width=width)
    except:
        #print('getNewt error')
        v, p, cj = 0.0, 0.0, cjest
    return v, p, cj, pos[1], pm, time


def newtonCJ(cjest, fuelv, fuelp, fgam, ashv, ashp, agam, width=0.8):
    """fits CJ velocity to a pair of states by varying the speed of the rayleigh line.
    
    Args:
        cjest(float): starting estimate for velocity.
        fuelv(float): fuel state specific volume.
        fuelp(float): fuel state pressure.
        fgam(float): fuel state sp. heat ratio.
        ashv(float): ash state specific volume.
        ashp(float): ash state pressure.
        agam(float): ash state sp. heat ratio.
        width(float): fitting sp. volume range for minimizing function.
        
    Returns:
        (float tuple): specific volume, pressure, CJVelocity.
    
    """
    miniu = minimize(fun=lambda x: diffHRupper(x, ph=fuelp, vh=fuelv, gh1=fgam, gh2=agam, 
                     pr=ashp, vr=ashv, env=[ashv*(1.0-width), ashv*(1.0+width)]), 
                     x0=cjest/2, tol=1e-14)
    minil = minimize(fun=lambda x: diffHRlower(x, ph=fuelp, vh=fuelv, gh1=fgam, gh2=agam, 
                     pr=ashp, vr=ashv, env=[ashv*(1.0-width), ashv*(1.0+width)]), 
                     x0=cjest/2, tol=1e-14)
    cjposu, cjspdu = miniu.fun[0], miniu.x[0]
    cjposl, cjspdl = minil.fun[0], minil.x[0]
    cjpos = 0.5*(cjposu+cjposl)
    cjpres = shockhugoniot(cjpos, p1=fuelp, v1=fuelv, g1=fgam, g2=agam)
    cjspd = rayleighSpeed(ashp, ashv, cjpres, cjpos)
    return cjpos, cjpres, cjspd


def getShockConditions(fname, inward=False, addvar='temp', direction=[]):
    """Returns bulk conditions at both sides of shock.
    Conditions are sorted so that output has the form: [ash, shock, fuel]
    
    Args:
        fname(str): filepath
        inward(bool): toggle for inward/outward bound shock.
        addvar(float): extra variable to get from lineout.
        **kwargs: arguments for lineout (direction, geometry, etc.)
    
    Returns:
        (list): radii of states.
        (list): densities of states.
        (list): pressures of states.
        (list): addvar at each state.
            (float): direct Rayleigh speed for the ash state.
        (float): match head position
        (float): timestamp of file.
    
    """
    fields = ['sound_speed', 'density', 'pressure', addvar]
    time, pars, _, _, paths = directMeta(fname)
    if len(direction)>(pars['dimensionality']-1):
        print("Direction doesn't match dimensionality: {}".format(pars['dimensionality']))
        return None
    data, _ = getLineout(fname, fields=fields, species=False, direction=direction, geom=pars['geometry'])
    #time, params, _, _, _ = directMeta(fname)
    rad, cs, dens, pres, var = data[0], data[1], data[2], data[3], data[4]
    
    # this fails for x_match = y_match = z_match = 0.0
    linepos = ut.estimateMatch(direction, pars, vvv=False)
    #print(linepos)
    shockin, shockout = ut.locateShock(rad, cs, linepos, vvv=False)
    #print('check shockConditions')
    if inward:
        ind = shockin
        offset = -1
    else:
        ind = shockout
        offset = 1
    # directly calculate the rayleigh speed joining both points
    cjest = rayleighSpeed(pres[ind+offset], 1.0/dens[ind+offset], 
                          pres[ind-offset], 1.0/dens[ind-offset])
    pos = [rad[ind-offset], rad[ind], rad[ind+offset]]
    condd = [dens[ind-offset], dens[ind], dens[ind+offset]]
    condp = [pres[ind-offset], pres[ind], pres[ind+offset]]
    xvar =  [var[ind-offset], var[ind], var[ind+offset]]
    return pos, condd, condp, xvar, cjest, linepos, time


def nonRelFermi(dens, ye=0.5):
    """ Completely degenerate, non-relativistic Fermi energy.
    E_f = (hbar^2/(2m_e))(3pi^(2/3))(N_a \rho Y_e)^(2/3)
    
    Args:
        dens(float): input density.
        ye(float): electron fraction.
    
    Returns:
        (float)
    
    """
    par1 = 0.5*np.power(0.5*h/np.pi, 2.0)/m_e
    par2 = np.power(3.0*np.pi*np.pi,2.0/3)
    par3 = np.power(Avogadro*dens*ye, 2.0/3)
    return par1*par2*par3


def extRelFermi(dens, ye=0.5):
    """ Completely degenerate, extreme-relativistic Fermi energy.
    E_f = hbar(3/8pi)^(1/3)(N_a \rho Y_e)^(1/3)
    
    Args:
        dens(float): input density.
        ye(float): electron fraction.
    
    Returns:
        (float)
    
    """
    par1 = h
    par2 = np.power(3.0/8.0/np.pi,1.0/3)
    par3 = np.power(Avogadro*dens*ye, 1.0/3)
    #print(par1, par2, par3)
    return par1*par2*par3


def isotherm(v, p0=1e23, v0=0.02):
    """calculate isoterm pressure passing through (p0, v0)
    at v.
    
    Args:
        v(float): input specific volume.
        p0(float): fixed pressure.
        v0(float): fixed sp. volume.
    
    Returns:
        (float)
    
    """
    return v0*p0/v


def adiabat(v, p0=1e23, v0=0.02, gamma=1.666):
    """calculate adiabat pressure passing through (p0, v0)
    at v, with gamma.
    
    Args:
        v(float): input specific volume.
        p0(float): fixed pressure.
        v0(float): fixed specific volume.
        gamma(float): fixed sp. heat ratio.
    
    Returns:
        (float)
    
    """
    num = p0*np.power(v0, gamma)
    denom = np.power(v, gamma)
    return num/denom


def rayleighSpeed(p1, v1, p2, v2):
    """returns Rayleigh line speed for a pair of points"""
    nom = p2-p1
    denom = v1-v2
    if denom <0:
        print ('Negative specific volume difference: dP {:.2e} dnu {:.2e}'.format(nom, denom))
        denom = 1.0
    return v2*np.sqrt(nom/denom)


def shockhugoniot(v2, p1=1e23, v1=0.02, g1=1.6666, g2=1.6666):
    """returns the huigoniot adiabat pressure corresponding to a 
    given specific volume while passing through a set point (v1,p1)."""
    g1fac = g1*v1/(g1-1)
    g2fac = g2*v2/(g2-1)
    var = 0.5*(v1+v2)
    return p1*(g1fac-var)/(g2fac-var)


def rayleigh(v2, p1=1e23, v1=0.02, speed=1e5):
    """returns the Rayleigh line pressure for a line crossing (v1,p1)."""
    sq = np.power(speed/v1, 2)
    fac = sq*(v1-v2)
    return p1+ fac


def diff(v, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666, pr=1e22, vr=0.02,speed=2.4e10):
    """yields the difference between the hugoniot adiabat and a rayleigh line,
    both passing through (v1, p1)"""
    hug = shockhugoniot(v, p1=ph, v1=vh, g1=gh1, g2=gh2)
    ray = rayleigh(v, p1=pr, v1=vr, speed=speed)
    return abs(ray)- abs(hug)


def customFormatter(factor, prec=1, width=2):
    """create a mpl formatter which factors labels by 10^factor 
    for clearer axes labels.
    """
    fstr = '{:{width}.{prec}f}'
    exp = 10.0**factor
    return FuncFormatter(lambda x, pos:fstr.format(x/exp, width=width, prec=prec))


def diffHRupper(sp, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666, pr=1e22, vr=0.02, env=[0.04, 0.08]):
    """minimizes hugoniot adiabat - rayleigh line difference, starting near a strong detonation.
    (starting from the "top")."""
    w = newton(func=lambda x: diff(x, ph=ph, vh=vh, gh1=gh1, gh2=gh2, pr=pr, vr=vr, speed=sp), x0=env[0], tol=1e-14)
    return w


def diffHRlower(sp, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666, pr=1e22, vr=0.02, env=[0.04, 0.08]):
    """minimizes hugoniot adiabat - rayleigh line difference, starting near a weak detonation.
    (starting from the "bottom")."""
    s = newton(func=lambda x: diff(x, ph=ph, vh=vh, gh1=gh1, gh2=gh2, pr=pr, vr=vr, speed=sp), x0=env[1], tol=1e-10)
    return s
