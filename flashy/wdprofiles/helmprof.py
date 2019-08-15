"""Direct hydrostatic eq subject to helmholtz eos.
"""
from ..utils import msol, np
import helmholtz as hh

_targetmass = msol
_densmin = 1e-11  # (2019/03/28)
_presmin = 2.05e+13


# hydro, jacobian, eos calls and bc
def setBC(dens, temp=1e7, abar=4.0, zbar=2.0, start=1e4):
    """starting point for solver."""
    helmobj = hh.helmeos(dens, temp, abar, zbar)
    # intial conditions: first cell's mass, and central pressure
    con1 = 4.0e0 * np.pi
    ms0 = con1 * start**3 * dens
    ps0 = helmobj.ptot[0] - 0.5e0 * con1 * G * start**2 * dens**2
    return [ms0, ps0]


def jac(x, y, abar=4.0, zbar=2.0, temp=1e7):
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_helm(y[1], deng=y[0]/vol, abar=abar,
                            zbar=zbar, temp=temp)
    mdm = 0
    mdp = 0
    pdm = -G*den/x/x
    pdp = 0
    return np.array([[mdm, mdp], [pdm, pdp]])


def derv(x, y, abar=4.0, zbar=2.0, temp=1e7, genrel=True):
    # this routine sets up the continuity and hydrostatic equilibrium ode's.
    # x is the radial coordinate, y(1) is the gravitational mass,
    # y(2) is the pressure
    con1 = 4.0e0*np.pi
    c2 = c*c
    # map the input vector
    massg = y[0]
    pres = y[1]

    # guess through mean density
    vol = 4*np.pi*np.power(x, 3.0)/3.0
    den, dpdd = invert_helm(pres, deng=y[0]/vol, abar=abar,
                            zbar=zbar, temp=temp)

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


def helm(dens, abar=4.0, zbar=2.0, temp=1e7):
    helmobj = hh.helmeos(dens, temp, abar, zbar)
    return helmobj.ptot[0], helmobj.dpd[0], helmobj.etot[0], helmobj.ded[0]


def invert_helm(pres, deng=4e6, abar=4.0, zbar=2.0, temp=1e7):
    """find the density and dpdd from a guess and a pressure."""
    # local variables
    steptol = 1.0e-6
    maxiter = 50
    deni = deng
    for i in range(maxiter):
        helmobj = hh.helmeos(deni, temp, abar, zbar)
        presi, dpd = helmobj.ptot[0], helmobj.dpd[0]
        z = abs((pres-presi)/pres)
        if z < steptol:
            return deni, dpd
        else:
            f = presi/pres - 1.0
            df = dpd/pres
            ratio = f/df
            deni = deni - ratio
    print("Reached max iterations ({:d})".format(maxiter))
    return deni, dpd


# stopping condition
def athens(t, y):
    athens.terminal = True
    athens.direction = 1.0
    if abs(y[0]-_targetmass) < 1e-6:
        print('reached 1msun')
        return 0.0
    else:
        return 1.0
