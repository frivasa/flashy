"""Polytropic profiles, Helmholtz temperature based on composition.
"""
# Giant Hammer
from scipy.integrate import solve_ivp
from ..utils import msol, rsol, G, c, gas_constant, np
from ..nuclear import convXmass2Abun
from ..datahaul.plainText import dataMatrix
_lown = 1.0e-10


def buildPolytropic(rstart, pc, rhoc, species=['he4'], xmass=[1.0],
                    gamma=4.0/3, rmax=rsol, mind=1e-5):
    """Generates a polytropic profile satisfying a given central pressure, density,
    and heat capacity ratio (index).
    Adds Temperatures based on composition through a Helmholtz free energy EoS.

    Args:
        rstart(float): starting radius for the solution.
        pc(float): starting pressure.
        rhoc(float): starting density.
        species(float): nuclides to include in composition (Helmholtz).
        xmass(float): mass fractions for each species (Helmholtz).
        gamma(float): polytrope index (1+1/n = gamma).
        rmax(float): furthest radius to solve the system.
        mind(float): curoff density for generated profile.

    Returns:
        dataMatrix: profile object with properties as attributes.

    """
    y0 = [rhoc*4.0*np.pi*np.power(rstart, 3.0)/3.0, pc]
    pcons = [pc/np.power(rhoc, gamma), gamma]
    pheidippides = solve_ivp(fun=lambda t, y: polyHydro(t, y, pcons),
                             method='LSODA',
                             jac=lambda t, y: jac(t, y, pcons),
                             t_span=(rstart, rmax), y0=y0, events=athens)
    rs, ms, ps = pheidippides.t, pheidippides.y[0], pheidippides.y[1]
    ds = np.array([polyDens(p, pcons) for p in ps])

    # trim off the edges
    ncut = np.where(ds > mind)[0][-1]
    rs, ms, ps, ds = rs[:ncut], ms[:ncut], ps[:ncut], ds[:ncut]

    ts = [polyTemp(p, d, species, xmass) for (p, d) in zip(ps, ds)]
    mult = len(rs)
    keys = ['radius', 'dens', 'pres', 'temp']
    dblock = np.column_stack([rs, ds, ps, ts])
    for i, x in enumerate(xmass):
        keys.append(species[i])
        dblock = np.column_stack((dblock, [x]*mult))
    return dataMatrix([keys, dblock])


def polyTemp(pres, rho, species, xmass):
    """returns the ideal gas temperature subject to a composition."""
    yis, abar, zbar = convXmass2Abun(species, xmass)
    molarmass = np.sum(yis)*1e-3  # turn to SI kg/mole
    # turn dyne to N and cm to m for gas_constant SI -> 1e-7
    return molarmass*pres/rho/gas_constant*1e-7


def polyDens(pres, pcons):
    arg = pres/pcons[0]
    pw = 1.0/float(pcons[1])
    dens = np.power(arg, pw)
    if dens < _lown:
        return 0.0
    else:
        return dens


def athens(t, y):
    athens.terminal = True
    if y[1] < _lown:
        return 0.0
    else:
        return 1.0


def jac(x, y, pcons):
    """Jacobian for polyH."""
    den = polyDens(y[1], pcons)
    mdm = 0
    mdp = 0
    pdm = -G*den/x/x
    pdp = 0
    return np.array([[mdm, mdp], [pdm, pdp]])


def polyHydro(x, y, pcons):
    """Returns the RHS for hydro+mass conservation subject to a
    Polytropic EoS
    """
    con1 = 4.0e0 * np.pi
    c2 = c*c
    # map the input vector
    massg = y[0]
    pres = y[1]

    # Poly dens
    den = polyDens(pres, pcons)
    if den < _lown:
        dydx = [0.0, 0.0]
    else:
        # d(massg)/dr
        dydx = [con1 * x*x * den]
        # d(press)/dr
        cor = (1.0 + pres/(den*c2)) *\
              (1.0 + (con1*pres*x**3)/(massg*c2)) /\
              (1.0 - (2.0*G*massg)/(x*c2))
        # cor = 1.0
        dydx.append(-G*massg/x/x*den*cor)
    return dydx
