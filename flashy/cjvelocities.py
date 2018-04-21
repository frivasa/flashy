"""cj velocity main module"""
# external package. leaving here for now...
from flashy.IOutils import sys
sys.path.append('../03.helmholtz/')
import helmholtz

import flashy.utils as ut
import flashy.datahaul.hdf5yt as reader
from flashy.nuclear import convXmass2Abun

def getVelocities(fname):
    """Returns positions of the shock, and both in an outer cj 
    velocities for a file, plus the starting x_match.
    (xin, xout, cjin, cjout, float(ray.ds.current_time), ray.ds.parameters['x_match'])
    """
    fields = ['sound_speed', 'density', 'pressure']
    data, _ = reader.profilePackage(fname, fields=fields, species=False)
    time, params, _, _ = reader.getMeta(fname)
    rad, cs, dens, pres = data[0], data[1], data[2], data[3]
    shockin, shockout = ut.locateShock(rad, cs, params['x_match'], vvv=False)
    xin, xout = rad[shockin], rad[shockout]
    cjin = ut.roughCJ(dens, pres, shockin)
    cjout = ut.roughCJ(dens, pres, shockout)
    return xin, xout, cjin, cjout, time, params['x_match']

# external starkiller package wrappers
def buildHelmTrojan(fname, offset=1):
    """Frankensteinian bridge between flash checkpoints and 
    J.Schwab's Helmholtz python module. 
    Joined inward/outward and spewing cj data to avoid calling 
    getRay more than once since it is a weak reference.
    
    Args:
        fname(str): filename.
        offset(int): zone offset from shock.
    
    """
    
    props = ['dens', 'temp', 'pres', 'eint']
    nprops = len(props)
    data, species = reader.profilePackage(fname, fields=props, species=True)
    nspecs = len(species)
    xin, xout, cjin, cjout, time, xmatch = getVelocities(fname)
    # get fuel and ash for outward shock
    # xin/xout == ray @ len([x for x in ray['r'][rsort] if x<xout])
    inw = len([x for x in data[0] if x<xout])+offset
    ouw = inw - 2*offset
    inv, ouv = [], []
    for i in range(nprops):
        inv.append(data[i][inw])
        ouv.append(data[i][ouw])
    xmin, xmou = [], []
    for i in range(nprops, nspecs):
        xmin.append(data[i][inw])
        xmou.append(data[i][ouw])
    _, abar, zbar = convXmass2Abun(species, xmin)
    fuelo = [[inv[-2], inv[-1]], [inv[0]], [inv[1]], [abar], [zbar]]
    _, abar, zbar = convXmass2Abun(species, xmou)
    asho = [[ouv[-2], ouv[-1]], [ouv[0]], [ouv[1]], [abar], [zbar]]
    
    # get fuel and ash for inward shock
    inw = len([x for x in data[0] if x<xin]) - offset
    ouw = inw + 2*offset
    inv, ouv = [], []
    for i in range(nprops):
        inv.append(data[i][inw])
        ouv.append(data[i][ouw])
    xmin, xmou = [], []
    for i in range(nprops, nspecs):
        xmin.append(data[i][inw])
        xmou.append(data[i][ouw])
    _, abar, zbar = convXmass2Abun(species, xmin)
    fueli = [[inv[-2], inv[-1]], [inv[0]], [inv[1]], [abar], [zbar]]
    _, abar, zbar = convXmass2Abun(species, xmou)
    ashi = [[ouv[-2], ouv[-1]], [ouv[0]], [ouv[1]], [abar], [zbar]]
    return fueli, ashi, fuelo, asho, [xin, xout, cjin, cjout, time, xmatch]


def getCJ(fname):
    fi, ai, fo, ao, cjdat = fp.buildHelmTrojan(fname, offset=1)
    # crank the abomination inwards
    ri = cj_cond(fi, ai, 0.0)
    if isinstance(ri, int):
        vcji = [0.0]
    else:
        vcji = 1.0/fi[1][0]*fp.np.sqrt((ri.ptot-fi[0][0])/(1.0/fi[1][0]-1.0/ri.den))
    #rcj, pcj = res.den, res.ptot
    # crank the abomination outwards
    ro = cj_cond(fo, ao, 0.0)
    if isinstance(ro, int):
        vcjo = [0.0]
    else:
        vcjo = 1.0/fo[1][0]*fp.np.sqrt((ro.ptot-fo[0][0])/(1.0/fo[1][0]-1.0/ro.den))
    #rcj, pcj = res.den, res.ptot
    # cjdat are my cj values ( xin, xout, cjin, cjout, time, xmatch )
    return vcji[0], vcjo[0], cjdat


def cj_cond(fuel, ash, q):
    """Hack-feeding data to helmholtz:
    fuel/ash = [[pressure, eint], [rho], [temp], [abar], [zbar]]
    """
    tol = 1e-8
    itmax = 40
    comprho = 1e10
    for it in range(itmax):
        res = helmholtz.helmeos(*ash[1:])
        res.den = fuel[1][0]*(1.0+(res.ptot-fuel[0][0])/(res.gam1*res.ptot))
        dden = res.den - comprho
        
        if abs(dden)<tol*res.den:
            return res
        elif res==-1:
            return -1
        else:
            ash[1][0] = res.den
            comprho = res.den
    return -1


def adiabat(fuel, ash, q):
    """Hack-feeding data to helmholtz:
    fuel/ash = [[pressure, eint], [rho], [temp], [abar], [zbar]]
    """
    tol = 1e-8
    itmax = 30
    for it in range(itmax):
        res = helmholtz.helmeos(*ash[1:])
        aux = (1.0/fuel[1][0] - 1.0/res.den)
        f = fuel[0][1] + q - res.etot + 0.5*(fuel[0][0]+res.ptot)*aux
        dfdT = -res.det + 0.5*res.dpt*aux
        dT = -f/dfdT
        if abs(dT)<tol*res.temp:
            return res
        else:
            ash[2][0] = ash[2][0] + dT[0]
    return -1