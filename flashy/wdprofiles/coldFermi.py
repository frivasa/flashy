"""
Fermi eos wd profiles, adapted from Timmes public cowd code.
EoS is a partially relativistic completely degenerate (cold) fermi gas
    (\eta\beta ~ 1 and \eta=\infty)
    (Cox & Giuli section 24.6c)
20180430: Helmholtz Eos added for temperature. Composition can also change.
"""
from ..utils import msol, rsol, G, c, h, m_e, Avogadro, np, kb
from ..nuclear import convXmass2Abun
from ..datahaul.plainText import dataMatrix
from ..datahaul.helmholtz import getTemps
from ..post import nonRelFermi
# Giant Hammer
from scipy.integrate import solve_ivp
# Tools for micromanagement
# from scipy.integrate import RK45
# mercury = RK45(derv, rads[0], [ms[0], ps[0]], rads[2], max_step=100)


def buildFermiHelmhotz(denc, xmass, species, pdens=0, hack=False):
    """Solves an IVP for a completely degenerate
    Fermi gas under hydrostatic equilibrium,
    then uses a Helmholtz EoS to assign a temperature
    value to the profile (this is somewhat
    wrong but better than guessing or putting an average).
    WARN: for high densities (over 1e8), helmholtz temp breaks
    (reaches minimum 1e4 K and 'bounces').
    WARN: pdens maxes around 600 points, beyond
    that the helmholtz package breaks.

    Args:
        denc(float): initial central density.
        xmass(float list): mass fractions for the used species.
        species(str list): nuclide code for each species.
        pdens(int): force a set number of points in
        the output (WARNING: runs ivp twice).
        hack(bool): avoid helmholtz temperature and
        return a 'weighted' fermi temperature.

    Returns:
        (dataMatrix): tabbable profile object.

    """
    ymass, abar, zbar = convXmass2Abun(species, xmass)
    ye = zbar/abar
    r, m, d, p = buildFermiProfile(denc, ye=ye, pdens=pdens)
    tH = getTemps([d[0]], [p[0]], [xmass], species)
    if hack:
        t = [nonRelFermi(de)/kb/tH for de in d]
    else:
        t = getTemps(d, p, len(d)*[xmass], species)
    keys = ['radius', 'dens', 'pres', 'temp']
    datablock = np.column_stack([r, d, p, t])
    mult = len(r)
    for i, x in enumerate(xmass):
        keys.append(species[i])
        datablock = np.column_stack((datablock, [x]*mult))
    return dataMatrix([keys, datablock])


def buildFermiProfile(denc, ye=0.5, start=1e4,
                      stop=1e10, pdens=400, debug=False):
    """Solves an IVP for a completely degenerate
    Fermi gas under hydrostatic equilibrium.

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
    y0 = setBC(denc, ye=ye, start=start)
    # run to find the edge
    pheidippides = solve_ivp(fun=lambda t, y: derv(t, y, ye=ye),
                             method='BDF', jac=jac,
                             t_span=(start, stop), y0=y0)  # events=athens)
    # run again to get the desired number of points and focus on the edge
    core = np.logspace(np.log10(start),
                       np.log10(0.9*pheidippides.t[-1]), pdens)
    # near the outermost 10%, use pdens points
    edge = np.linspace(0.91*pheidippides.t[-1], 1.01*pheidippides.t[-1], pdens)
    rads = np.append(core, edge)
    stsize = edge[-1]-edge[-2]
    pheidippides = solve_ivp(fun=lambda t, y: derv(t, y, ye=ye), jac=jac,
                             method='BDF',
                             max_step=stsize*0.5,
                             t_span=(start, stop), y0=y0,
                             t_eval=rads)
    rs, ms, ps = pheidippides.t, pheidippides.y[0], pheidippides.y[1]
    ds = [invert_fergas(denc, p, ye)[0] for p in ps]
    if debug:
        print(pheidippides['message'])
        print('{:e}'.format(pheidippides['t'][-1]))
        print('{:e}'.format(pheidippides['y'][0][-1]))
        return pheidippides, ds
    else:
        return rs, ms, ds, ps


def setBC(dens, ye=0.5, start=1e4):
    presc, dpresdd, ener, denerdd = fergas(dens, ye)
    # intial conditions: first cell's mass, and central pressure
    con1 = 4.0e0 * np.pi
    ms0 = con1 * start**3 * dens
    ps0 = presc - 0.5e0 * con1 * G * start**2 * dens**2
    return [ms0, ps0]


def jac(x, y, denc=1e9):
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_fergas(y[0]/vol, y[1], 0.5)
    mdm = 0
    mdp = 0
    pdm = -G*den/x/x
    pdp = 0
    return np.array([[mdm, mdp], [pdm, pdp]])


def derv(x, y, denc=1e9, ye=0.5, genrel=True):
    # this routine sets up the continuity and hydrostatic equilibrium ode's.
    # x is the radial coordinate, y(1) is the gravitational mass,
    # y(2) is the pressure
    con1 = 4.0e0 * np.pi
    c2 = c*c
    # map the input vector
    massg = y[0]
    pres = y[1]

    # cold ideal fermi gas
    # guess through mean density
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_fergas(y[0]/vol, pres, ye)
    # here is d(massg)/dr
    dydx = [con1 * x*x * den]
    # here is d(press)/dr
    if genrel:
        cor = (1.0 + pres/(den*c2)) *\
              (1.0 + (con1*pres*x**3)/(massg*c2)) /\
              (1.0 - (2.0*G*massg)/(x*c2))
    else:
        cor = 1.0
    dydx.append(-G*massg/x/x*den*cor)
    return dydx


def fergas(den, ye):
    """the eos for a partial relativistic completely degenerate (cold) fermi gas
    (\eta\beta ~ 1 and \eta=\infty)
    (Cox & Giuli section 24.6c)

    input is the density den and the mean charge to mean weight ratio ye.
    output is the pressure (in erg/cm**3),
    the pressure derivative with density,
    the energy (in erg/g), and the energy derivative with density.
    """
    hbar = 0.5*h/np.pi
    lamb = hbar/(m_e*c)
    lam3 = lamb*lamb*lamb
    xcon = 3.0e0 * np.pi*np.pi * lam3 * Avogadro
    pcon = m_e*c*c/(lam3*8.0e0*np.pi*np.pi)

    if (den < 0.0):
        # print 'bad pass in routine fergas'
        return 0.0, 0.0, 0.0, 0.0
    deni = 1.0e0/den
    x = np.power(xcon*den*ye, 1/3.0)
    dxdd = x*deni/3.0
    x2 = x*x
    x3 = x2*x
    # the pressure in erg/cm**3
    # note: fac3 is a way of writing arc-sinh(x)
    fac1 = np.sqrt(1.0e0 + x2)
    fac2 = 2.0*x2/3.0 - 1.0
    fac3 = np.log(x + fac1)
    pres = pcon * (x*fac1*fac2 + fac3)
    # pressure derivative with density
    dfac1 = x/fac1
    dfac2 = 4.0*x/3.0
    dfac3 = (1.0 + dfac1)/(x + fac1)
    dpdx = pcon * (fac1*fac2 + x*dfac1*fac2 + x*fac1*dfac2 + dfac3)
    dpdd = dpdx * dxdd
    # the internal energy in erg/cm**3
    gac1 = 8.0 * x3/3.0
    gac2 = fac1 - 1.0
    ener = pcon*gac1*gac2 - pres
    # energy derivative with density
    dgac1 = 8.0 * x2
    dgac2 = dfac1
    dedx = pcon*(dgac1*gac2 + gac1*dgac2) - dpdx
    dedd = dedx * dxdd
    # convert the energies into erg/g
    ener = ener * deni
    dedd = (dedd - ener) * deni
    return pres, dpdd, ener, dedd


def invert_fergas(den, pres, ye):
    """given the pressure, ye, and a guess for the density,
    find the density and dpdd
    """
    # local variables
    eostol = 1.0e-6
    fpmin = 1.0e-14
    # save the initial guess
    denold = den
    # newton loop
    for i in range(100):
        pres1, dpdd, ener, dedd = fergas(den, ye)
        f = pres1/pres - 1.0
        df = dpdd/pres
        if (df == 0.0):
            return denold, 0.0
        ratio = f/df
        dennew = den - ratio
        z = abs((dennew - den)/den)
        den = dennew
        if (z < eostol or abs(ratio) <= fpmin):
            break
    # call it one more time with the converged value of the density
    pres1, dpdd, ener, dedd = fergas(den, ye)
    return den, dpdd
