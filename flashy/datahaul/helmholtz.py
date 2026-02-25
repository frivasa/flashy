import helmholtz
import flashy.datahaul.hdf5yt as hdf5yt
import numpy as np
from flashy.post import getRayleighVelocities
from flashy.nuclear import conv_xmass2abun


# external starkiller/Helmholtz package wrappers
def getTemps(rhos, pres, xmasses, species, returnObj=False):
    trojan = wrapVector(rhos, pres, xmasses, species)
    sol = helmholtz.helmeos_DP(*trojan)
    if returnObj:
        return sol
    else:
        return sol.temp


def getPres(rhos, temps, xmasses, species, returnObj=False):
    trojan = wrapVector(rhos, temps, xmasses, species)
    sol = helmholtz.helmeos(*trojan)
    if returnObj:
        return sol
    else:
        return sol.ptot


def getExtHelmCJ(fname, geom='spherical'):
    fi, ai, fo, ao, cjdat = buildHelmTrojan(fname, offset=1, geom=geom)
    # crank the abomination inwards
    # ri = adiabat(fi, ai, 0.0)
    ri = cj_cond(fi, ai)
    if isinstance(ri, int):
        vcji = [0.0]
    else:
        denom = 1.0/fi[1][0]
        vcji = denom*np.sqrt((ri.ptot-fi[0][0])/(1.0/fi[1][0]-1.0/ri.den))
    # rcj, pcj = res.den, res.ptot
    # crank the abomination outwards
    # ro = adiabat(fo, ao, 0.0)
    ro = cj_cond(fo, ao)
    if isinstance(ro, int):
        vcjo = [0.0]
    else:
        denom = 1.0/fo[1][0]
        vcjo = denom*np.sqrt((ro.ptot-fo[0][0])/(1.0/fo[1][0]-1.0/ro.den))
    # rcj, pcj = res.den, res.ptot
    # cjdat are my cj values ( xin, xout, cjin, cjout, time, xmatch )
    return vcji[0], vcjo[0], cjdat


def wrapVector(rho, var, xmass, species):
    """
    Wrapper for J.Schwab's Helmholtz python module, which is
    itself also one for Timmes' Helmholtz EoS.
    Checks rho for length and comparison so that it can be
    used for single points.

    Args:
        rho(float): query density.
        var(float): auxilliary thermodynamic variable (T, P, S or E(int)).
        xmass(float): query mass fractions.
        species(str): query nuclide list.

    Returns:
        (trojan): [[rho/s], [var/s], [abar/s], [zbar/s]]

    """
    try:
        pnts = len(rho)
    except TypeError:
        _, abar, zbar = conv_xmass2abun(species, xmass)
        return [[rho], [var], [abar], [zbar]]
    allp = []
    for p in range(pnts):
        _, abar, zbar = conv_xmass2abun(species, xmass[p])
        allp.append((rho[p], var[p], abar, zbar))
    vr, vv, va, vz = zip(*allp)
    return [vr, vv, va, vz]


def cj_cond(fuel, ash):
    """Calculate Cj velocity from Helmholtz Eos calculation
    fits density for a set temperature

    Args:
        fuel (trojan): [[pressure, eint], [rho], [temp], [abar], [zbar]]
        ash (trojan): [[pressure, eint], [rho], [temp], [abar], [zbar]]

    """
    tol = 1e-8
    itmax = 30
    comprho = 1e10
    for it in range(itmax):
        # res = adiabat(fuel, ash, q)
        res = helmholtz.helmeos(*ash[1:])
        res.den = fuel[1][0]*(1.0+(res.ptot-fuel[0][0])/(res.gam1*res.ptot))
        # print res.den, type(res.den[0])
        dden = res.den[0] - comprho

        if abs(dden) < tol*res.den:
            return res
        elif res == -1:
            return -1
        else:
            ash[1][0] = res.den
            comprho = res.den
    return -1


def adiabat(fuel, ash, q):
    """Hack-feeding data to helmholtz:
    fuel/ash = [[pressure, eint], [rho], [temp], [abar], [zbar]]

    fits a temperature to a hugoniot curve
    """
    # q value -- we need the change in molar fractions
    # call ener_gener_rate(eos_state_ash
    # % xn(:)/aion(:) - eos_state_fuel % xn(:)/aion(:), q_burn)
    tol = 1e-8
    itmax = 30
    for it in range(itmax):
        res = helmholtz.helmeos(*ash[1:])
        aux = (1.0/fuel[1][0] - 1.0/res.den)
        f = fuel[0][1] + q - res.etot + 0.5*(fuel[0][0]+res.ptot)*aux
        dfdT = -res.det + 0.5*res.dpt*aux
        dT = -f/dfdT
        if abs(dT) < tol*res.temp:
            return res
        else:
            ash[2][0] = ash[2][0] + dT[0]
    return -1


def buildHelmTrojan(fname, offset=1, geom='spherical'):
    """Frankensteinian bridge between flash checkpoints and
    J.Schwab's Helmholtz python module.
    Joined inward/outward and spewing cj data to avoid
    calling yt more than once.

    Args:
        fname(str): filename.
        offset(int): zone offset from shock.

    """
    props = ['dens', 'temp', 'pres', 'eint']
    nprops = len(props)
    data, species = hdf5yt.getLineout(fname, fields=props,
                                      species=True, geom=geom)
    nspecs = len(species)
    xin, xout, cjin, cjout, time, xmatch = getRayleighVelocities(fname)
    # get fuel and ash for outward shock
    # xin/xout == ray @ len([x for x in ray['r'][rsort] if x<xout])
    inw = len([x for x in data[0] if x < xout])+offset
    ouw = inw - 2*offset
    inv, ouv = [], []
    for i in range(1, nprops+1):
        inv.append(data[i][inw])
        ouv.append(data[i][ouw])
    xmin, xmou = [], []
    for i in range(nprops+1, nprops+1+nspecs):
        xmin.append(data[i][inw])
        xmou.append(data[i][ouw])
    _, abar, zbar = conv_xmass2abun(species, xmin)
    # [pressure, eint], [rho], [temp], [abar], [zbar]]
    fuelo = [[inv[-2], inv[-1]], [inv[0]], [inv[1]], [abar], [zbar]]
    _, abar, zbar = conv_xmass2abun(species, xmou)
    asho = [[ouv[-2], ouv[-1]], [ouv[0]], [ouv[1]], [abar], [zbar]]

    # get fuel and ash for inward shock
    inw = len([x for x in data[0] if x < xin]) - offset
    ouw = inw + 2*offset
    inv, ouv = [], []
    for i in range(1, nprops+1):
        inv.append(data[i][inw])
        ouv.append(data[i][ouw])
    xmin, xmou = [], []
    for i in range(nprops+1, nprops+1+nspecs):
        xmin.append(data[i][inw])
        xmou.append(data[i][ouw])
    _, abar, zbar = conv_xmass2abun(species, xmin)
    fueli = [[inv[-2], inv[-1]], [inv[0]], [inv[1]], [abar], [zbar]]
    _, abar, zbar = conv_xmass2abun(species, xmou)
    ashi = [[ouv[-2], ouv[-1]], [ouv[0]], [ouv[1]], [abar], [zbar]]
    return fueli, ashi, fuelo, asho, [xin, xout, cjin, cjout, time, xmatch]
