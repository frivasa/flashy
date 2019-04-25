"""Direct hydrostatic eq subject to helmholtz eos.
Adding temp... this should be redundant.
"""
from ..utils import msol, rsol, G, c, h, m_e, Avogadro, np, kb
from ..nuclear import convXmass2Abun
from ..datahaul.plainText import dataMatrix
from ..datahaul.helmholtz import getTemps
from ..post import nonRelFermi
import helmholtz as hh
# Giant Hammer
from scipy.integrate import solve_ivp

def buildHelmholtz(denc, temp, xmass, species, pdens=0, hack=False):
    """Solves an IVP for a ???
    WARN: for high densities (over 1e8), helmholtz temp breaks 
    (reaches minimum 1e4 K and 'bounces').
    WARN: pdens maxes around 600 points, beyond that the helmholtz package breaks.

    Args:
        denc(float): initial central density.
        xmass(float list): mass fractions for the used species.
        species(str list): nuclide code for each species.
        pdens(int): force a set number of points in the output (WARNING: runs ivp twice).
        hack(bool): avoid helmholtz temperature and return a 'weighted' fermi temperature.

    Returns:
        (dataMatrix): tabbable profile object.

    """
    ymass, abar, zbar = convXmass2Abun(species, xmass)
    print('Pre BP')
    r, m, d, p = buildHelmholtzProfile(denc, abar=abar, zbar=zbar, pdens=pdens, temp=temp)
    print('post buildProf')
    t = getTemps(d, p, len(d)*[xmass], species)
    # otp dmatr
    keys = ['radius', 'dens', 'pres', 'temp']
    datablock = np.column_stack([r, d, p, t])
    mult = len(r)
    for i, x in enumerate(xmass):
        keys.append(species[i])
        datablock = np.column_stack((datablock, [x]*mult))
    return dataMatrix([keys, datablock])


def buildHelmholtzProfile(denc, abar=1.0, zbar=0.5, start=1e4, stop=1e10, pdens=400, debug=False, temp=1e7):
    """Solves an IVP for 
    a completely degenerate Fermi gas under hydrostatic equilibrium.

    Args:
        denc(float): initial central density.
        ye(float): electron fraction.
        start(float): initial radius (minimum radius).
        stop(float): ivp solver stopping radius.
        pdens(int): points/2 to evaulate. (uses other half for edge)

    Returns:
        (float list): radii
        (float list): masses
        (float list): densities
        (float list): pressures

    """
    y0 = setBC(denc, abar=abar, zbar=zbar, start=start, temp=temp)
    # run to find the edge
    pheidippides = solve_ivp(fun=lambda t, y: derv(t, y, abar=abar, zbar=zbar, temp=temp), method='BDF', jac=jac,
                             t_span=(start, stop), y0=y0)# events=athens)
    # run again to get the desired number of points and focus on the edge
    print('first run')
#     core = np.logspace(np.log10(start), np.log10(0.9*pheidippides.t[-1]), pdens)
#     # near the outermost 10%, use pdens points
#     edge = np.linspace(0.91*pheidippides.t[-1], 1.01*pheidippides.t[-1], pdens)
#     rads = np.append(core, edge)
#     stsize = edge[-1]-edge[-2]
#     pheidippides = solve_ivp(fun=lambda t, y: derv(t, y, ye=ye), jac=jac,
#                              method='BDF',
#                              max_step=stsize*0.5,
#                              t_span=(start, stop), y0=y0,
#                              t_eval=rads)
    rs, ms, ps = pheidippides.t, pheidippides.y[0], pheidippides.y[1]
    print('second run finished')
    ds = [invert_helm(denc, p, abar=abar, zbar=zbar, temp=temp)[0] for p in ps]
    print('final run, temps')
    if debug:
        print(pheidippides['message'])
        print('{:e}'.format(pheidippides['t'][-1]))
        print('{:e}'.format(pheidippides['y'][0][-1]))
        return pheidippides, ds
    else:
        return rs, ms, ds, ps


def setBC(dens, abar=4.0, zbar=2.0, start=1e4, temp=1e7):
    """starting point for solver."""
    helmobj = hh.helmeos(dens, temp, abar, zbar)
    # intial conditions: first cell's mass, and central pressure
    con1   = 4.0e0 * np.pi
    ms0 = con1 * start**3 * dens
    ps0 = helmobj.ptot - 0.5e0 * con1 * G * start**2 * dens**2
    return [ms0, ps0]


# def athens(t, y):
#     athens.terminal = True
#     if y[1]<=-1e6:
#         return 0.0
#     else:
#         return 1.0


def jac(x, y, denc=1e9, abar=4.0, zbar=2.0, temp=1e7):
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_helm(y[0]/vol, y[1], abar=abar, zbar=zbar, temp=temp)
    mdm = 0
    mdp = 0
    pdm = -G*den/x/x
    pdp = 0
    return np.array([[mdm, mdp],[pdm, pdp]])


def derv(x, y, denc=1e9, abar=4.0, zbar=2.0, temp=1e7, genrel=True):
    # this routine sets up the continuity and hydrostatic equilibrium ode's.
    # x is the radial coordinate, y(1) is the gravitational mass,
    # y(2) is the pressure
    con1   = 4.0e0 * np.pi
    c2     = c*c
    # map the input vector
    massg = y[0]
    pres  = y[1]

    # cold ideal fermi gas
    # guess through mean density
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_helm(y[0]/vol, pres, abar=abar, zbar=zbar)
    
    # here is d(massg)/dr
    dydx = [con1 * x*x * den]
    # here is d(press)/dr
    if genrel:
        cor = (1.0 + pres/(den*c2)) *\
              (1.0 + (con1*pres*x**3)/(massg*c2)) /\
              (1.0 - (2.0*G*massg)/(x*c2))
    else:
        cor = 1.0
    dydx.append(-G * massg/x/x* den * cor)
    return dydx


def helm(dens, abar=4.0, zbar=2.0, temp=1e7):
    """the eos for a partial relativistic completely degenerate (cold) fermi gas
    (\eta\beta ~ 1 and \eta=\infty)
    (Cox & Giuli section 24.6c)
    
    input is the density dens and the mean charge to mean weight ratio ye.
    output is the pressure (in erg/cm**3), the pressure derivative with density,
    the energy (in erg/g), and the energy derivative with density.
    """
    helmobj = hh.helmeos(dens, temp, abar, zbar)
    if (dens < 0.0):
        # print 'bad pass in routine fergas'
        return 0.0, 0.0, 0.0, 0.0
    return helmobj.ptot, helmobj.dpdd, helmobj.etot, helmobj.dedd


def invert_helm(den, pres, abar=4.0, zbar=2.0, temp=1e7):
    """
    given the pressure, ye, and a guess for the density,
    find the density and dpdd
    """
    # local variables
    eostol = 1.0e-6
    fpmin  = 1.0e-14
    # save the initial guess
    denold = den
    # newton loop
    for i in range(100):
        pres1, dpdd, ener, dedd = helm(den, abar=abar, zbar=zbar, temp=temp)
        f  = pres1/pres - 1.0
        df = dpdd/pres
        if (df == 0.0):
            return denold, 0.0
        ratio  = f/df
        dennew = den - ratio
        z      = abs((dennew - den)/den)
        den    = dennew
        if (z < eostol or abs(ratio) <= fpmin):
            break
    # call it one more time with the converged value of the density
    pres1, dpdd, ener, dedd = helm(den, abar=abar, zbar=zbar, temp=temp)
    return den, dpdd
