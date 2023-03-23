"""post processing for otp data."""
from flashy.datahaul.hdf5yt import getLineout, wedge2d, wedge3d
from flashy.datahaul.hdfdirect import directMeta
from flashy.datahaul.ytfields import _alphas
import flashy.utils as ut
from .utils import h, m_e, np, c, kb, Avogadro
from scipy.optimize import newton, minimize
import flashy.datahaul.parData as pd
import flashy.paraMan as pman
from .IOutils import os


def get2Dtaus(bview, fname, wedges=5):
    """calculates soundspeed timescale and burning timescale
    for every cell in the domain.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file name.
        wedges(int): parallel extraction slices.

    Returns:
        (np.array): [tauC, tauE]

    """
    # get data from file in wedges and join everything into 'dat'
    kwargs = {'fname': os.path.abspath(fname), 'wedges': wedges,
              'fields': ['dx', 'gamc', 'pres', 'dens', 'eint', 'enuc']}
    res = pman.throwHammer(bview, wedges, pd.par_wedge2d, **kwargs)
    dat = pd.glue2dWedges(res.get())
    # work on data
    soundspeeds = np.sqrt(dat[1]*dat[2]/dat[3])
    tauC = dat[0]/soundspeeds
    tauE = dat[4]/dat[5]
    filter = np.logical_not(np.isnan(tauE))
    tauE = tauE[filter]
    tauC = tauC[filter]
    # output as a plot friendly zip
    return dat[0][filter], tauC, tauE


def par_radialSpeeds(wedgenum, wedges=5, fname='', geom='cartesian',
                     dimension=2, ref='x', avoid=0.0, antipode=False):
    """wedge-parallelizable radial speed extraction."""
    # slight offset to avoid division by zero.
    if avoid:
        cuts = np.linspace(179.9 - avoid, 0.1 + avoid, wedges + 1) - 90
    else:
        cuts = np.linspace(179.9, 0.1, wedges + 1) - 90
    wedges = list(zip(cuts, cuts[1:]))
    start, stop = wedges[wedgenum]
    if stop*start > 0.0:
        stop = -stop
    else:
        start = abs(start)
        stop = abs(stop)
    print(start, stop)
    return radialSpeeds(fname, elevation=start, depth=stop,
                        antipode=antipode, geom=geom,
                        dimension=dimension, ref=ref)


def radialSpeeds(fname, elevation=5, depth=5,
                 geom='cartesian', dimension=2, ref='x', antipode=False):
    """radius vs radial speed for a wedge centered at the equator.
    """
    if dimension == 2:
        rawd = wedge2d(fname, elevation, depth,
                       fields=['x', 'y', 'z', 'velx', 'vely', 'velz'])
        rs, vrs = [], []
        for x, y, z, vx, vy, vz in zip(*rawd):
            vec = np.array([x, y, z])
            vel = np.array([vx, vy, vz])
            r = np.sqrt(vec.dot(vec))
            v = np.sqrt(vel.dot(vel))
            rxv = np.cross(vec, vel)
            normrxv = np.sqrt(rxv.dot(rxv))
            angle = np.arcsin(normrxv/r/v)
            vrad = v*np.cos(angle)
            rs.append(r)
            vrs.append(vrad)
    elif dimension == 3:
        rawd = wedge3d(fname, elevation, depth,
                       fields=['spherical_radius',
                               'velocity_spherical_radius'],
                       reference=ref, antipode=antipode)
        rs, vrs = rawd
    else:
        rawd, _ = getLineout(fname, fields=['velx'], geom=geom, species=False)
        rs, vrs = rawd
    return rs, vrs


def par_speedHisto(wedgenum, wedges=5, fname='', geom='cartesian',
                   resolution=4e7, dimension=2, ref='x',
                   antipode=False, cylvel=False, species=_alphas):
    """wedge-parallelizable speedHisto."""
    # slight offset to avoid division by zero.
    cuts = np.linspace(179.9, 0.1, wedges+1)-90
    wedges = list(zip(cuts, cuts[1:]))
    start, stop = wedges[wedgenum]
    if stop*start > 0.0:
        stop = -stop
    else:
        start = abs(start)
        stop = abs(stop)
    print(start, stop)
    return speedHisto(fname, resolution=resolution, velrange=[1e9, 5e9],
                      elevation=start, depth=stop, geom=geom,
                      dimension=dimension, ref=ref, cylvel=cylvel,
                      antipode=antipode, species=species)


def speedHisto2d(fname, radius, sphereCapR=1e9, ref='y',
                 resolution=4e7, velrange=[1e9, 5e9],
                 geom='cylindrical', species=_alphas):
    """Calculate speeds within a cylinder oriented towards an axis 
    and return masses.
    histogram bin range is fixed
    # recommended resolution: 100 km/s (Fink, 2010)
    # default resolution: 400 km/s
    # default max velocity: a sixth of c.

    Args:
        fname(str): filename.
        radius(float): half diameter of borehole cylinder.
        sphereCapR(float): max radius to consider (sphere cutoff).
        ref(str): reference axis.
        resolution(float): bin size in cm/s.
        geom(str): specify geometry for 1d file.
        fields(str list): specify which fields of data to get.
        species(str list): specify which species data to get.

    Returns:
        np.array list, np.array: Masses for species queries, bin limits.

    """
    rawd, ispecies = wedge2d(fname, radius, sphereCapR=sphereCapR,
                             reference=ref) # , fields=fields)
    offset = len(rawd)-len(ispecies)
    pick = { 'x': 0, 'y': 1, 'z': 2 }[ref]
    speeds = rawd[pick][:]/1e5
    print('using {} = position '.format(ref), pick)
    print('cells picked:', len(speeds))
    print('getting mass of each species per cell....')
    celltpls = []
    for i in range(len(rawd[0])):
        # multiply each X_i by the cell_masses
        masses = [rawd[j][i]*rawd[1][i]
                  for j in range(offset, len(ispecies)+offset)]
        celltpls.append((speeds[i], masses))
    
    print('sorting by speed...')
    sortedcells = sorted(celltpls)
    speeds, massgrid = zip(*sortedcells)
    # massgrid = [[Mh1_0, Mhe4_0,...], [Mh1_1, Mhe4_1,...]]
    # respecting increasing speed order

    # get ranges for histogram bins
    print("data range {:E} {:E}".format(np.max(speeds), np.min(speeds)))
    vmin, vmax = velrange
    binnum = int((vmax - vmin)/resolution)
    print("post.speedHisto: Velocity range: {:e} {:e}".format(vmin, vmax))
    print("post.speedHisto: Bins: {}".format(binnum))
    try:
        counts, bins = np.histogram(speeds, bins=int(binnum))
    except ValueError:
        print('post.speedHisto: zero speed range.')
        return None

    # sort by species
    species = [s.strip() for s in species]
    buckets = []
    for s in species:
        buckets.append(np.zeros(len(bins)))
    psp = dict(zip(species, range(len(species))))
    for i in range(len(ispecies)):
        if ispecies[i] in species:  # main species filter
            weights = [0]
            start = 0
            for c in counts:
                # since the data is already sorted by speed
                # counts follows the data in each bin
                totalspmass = sum([m[i] for m in massgrid[start:start+c]])
                # force the summed data value as histo "weight"
                weights.append(totalspmass)
                start += c
            weights = np.array(weights)
            # finally sum the array of masses to the corresponding histo
            buckets[psp[ispecies[i]]] += weights
        else:
            continue
    return buckets, bins


def speedHisto3d(fname, radius, sphereCapR=1e9, ref='z',
                 resolution=4e7, velrange=[1e9, 5e9],
                 geom='cartesian', dimension=3,
                 fields=['velx', 'vely', 'velz', 'cell_mass'], 
                 species=_alphas):
    """Calculate speeds within a cylinder oriented towards an axis 
    and return masses.
    histogram bin range is fixed
    # recommended resolution: 100 km/s (Fink, 2010)
    # default resolution: 400 km/s
    # default max velocity: a sixth of c.

    Args:
        fname(str): filename.
        radius(float): half diameter of borehole cylinder.
        sphereCapR(float): max radius to consider (sphere cutoff).
        ref(str): reference axis.
        resolution(float): bin size in cm/s.
        geom(str): specify geometry for 1d file.
        dimension(int): specify file dimension.
        fields(str list): specify which fields of data to get.
        species(str list): specify which species data to get.

    Returns:
        np.array list, np.array: Masses for species queries, bin limits.

    """
    # fields = fields + species
    rawd, ispecies = wedge3d(fname, radius, sphereCapR=sphereCapR,
                             reference=ref)
    offset = len(rawd)-len(ispecies)
    pick = { 'x': 0, 'y': 1, 'z': 2 }[ref]
    speeds = rawd[pick][:]/1e5  # take z velocities, do km/s
    print('using {} = position '.format(ref), pick)
    print('cells picked:', len(speeds))
    print('getting mass of each species per cell....')
    celltpls = []
    for i in range(len(rawd[0])):
        # multiply each X_i by the cell_masses
        masses = [rawd[j][i]*rawd[1][i]
                  for j in range(offset, len(ispecies)+offset)]
        celltpls.append((speeds[i], masses))
    
    print('sorting by speed...')
    sortedcells = sorted(celltpls)
    speeds, massgrid = zip(*sortedcells)
    # massgrid = [[Mh1_0, Mhe4_0,...], [Mh1_1, Mhe4_1,...]]
    # respecting increasing speed order

    # get ranges for histogram bins
    print("data range {:E} {:E}".format(np.max(speeds), np.min(speeds))) 
    vmin, vmax = velrange
    binnum = int((vmax - vmin)/resolution)
    print("post.speedHisto: Velocity range: {:e} {:e}".format(vmin, vmax))
    print("post.speedHisto: Bins: {}".format(binnum))
    try:
        counts, bins = np.histogram(speeds, bins=int(binnum))
    except ValueError:
        print('post.speedHisto: zero speed range.')
        return None

    # sort by species
    # F5+ no longer uses 'he4' but 'he4 ' (with ' ')
    # species = [s.strip() for s in species]
    buckets = []
    for s in species:
        buckets.append(np.zeros(len(bins)))
    psp = dict(zip(species, range(len(species))))
    for i in range(len(ispecies)):
        if ispecies[i] in species:  # main species filter
            weights = [0]
            start = 0
            for c in counts:
                # since the data is already sorted by speed
                # counts follows the data in each bin
                totalspmass = sum([m[i] for m in massgrid[start:start+c]])
                # force the summed data value as histo "weight"
                weights.append(totalspmass)
                start += c
            weights = np.array(weights)
            # finally sum the array of masses to the corresponding histo
            buckets[psp[ispecies[i]]] += weights
        else:
            continue
    return buckets, bins


def speedHisto_legacy(fname, resolution=4e7, velrange=[1e9, 5e9],
               elevation=5, depth=5, geom='cartesian', dimension=2,
               ref='x', antipode=False, cylvel=False,
               maxR=1.0e9, fields=['velx', 'vely', 'velz', 'cell_mass'], 
               species=_alphas):
    """Calculate speeds within a wedge and return masses.
    histogram bin range is fixed so that one can mix wedges.
    # recommended resolution: 100 km/s (Fink, 2010)
    # default resolution: 400 km/s
    # default max velocity: a sixth of c.

    Args:
        fname(str): filename.
        resolution(float): bin size in cm/s.
        velrange(float list): historgram range.
        elevation(float): equator-north pole degree.
        depth(float): equator-south pole degree.
        geom(str): specify geometry for 1d file.
        dimension(int): specify file dimension.
        ref(str): reference axis (3D only, see datahaul.hdf5yt.wedge3d).
        antipode(bool): antipodal wedge (see datahaul.hdf5yt.wedge3d).
        cylvel(bool): take velx (axial velocity in cylindrical) as speed.
        maxR(float): sphere cap for cutout (for large 3d files)
        species(str list): specify which species data to get.

    Returns:
        np.array list, np.array: Masses for species queries, bin limits.

    """
    # split into wedges process each, then read outputs...
    # wedge2d call without 'fields' calculates cylindrical volumes
    if dimension == 2:
        rawd, ispecies = wedge2d(fname, elevation, depth, cylvel=cylvel)
    elif dimension == 3:
        fields = fields + species
        rawd, ispecies = wedge3d(fname, elevation, depth, maxR=maxR,
                                 fields=fields,
                                 reference=ref, antipode=antipode)
    else:
        rawd, ispecies = getLineout(fname, fields=['velx', 'vely',
                                                  'velz', 'cell_mass'],
                                   geom=geom)
        rawd = rawd[1:]  # remove the radius column

    offset = len(rawd)-len(ispecies)
    # see datahaul.hdf5yt.wedge2d for positions of data in rawd
    # speeds = rawd[0][:]
    pick = { 'x': 0, 'y': 1, 'z': 2 }[ref]
    speeds = rawd[pick][:]  # take z velocities
    print('using {} = position '.format(ref), pick)
    print('cells picked:', len(speeds))
    print('getting mass of each species per cell....')
    celltpls = []
    for i in range(len(rawd[0])):
        # multiply each X_i by the cell_masses
        masses = [rawd[j][i]*rawd[1][i]
                  for j in range(offset, len(ispecies)+offset)]
        celltpls.append((speeds[i], masses))
    
    print('sorting by speed...')
    sortedcells = sorted(celltpls)
    speeds, massgrid = zip(*sortedcells)
    # massgrid = [[Mh1_0, Mhe4_0,...], [Mh1_1, Mhe4_1,...]]
    # respecting increasing speed order

    # get ranges for histogram bins
    vmin, vmax = velrange
    binnum = int((vmax - vmin)/resolution)
    print("post.speedHisto: Velocity range: {:e} {:e}".format(vmin, vmax))
    print("post.speedHisto: Bins: {}".format(binnum))
    try:
        counts, bins = np.histogram(speeds, bins=int(binnum))
    except ValueError:
        print('post.speedHisto: zero speed range.')
        return None

    # sort by species
    species = [s.strip() for s in species]
    buckets = []
    for s in species:
        buckets.append(np.zeros(len(bins)))
    psp = dict(zip(species, range(len(species))))
    for i in range(len(ispecies)):
        if ispecies[i] in species:  # main species filter
            weights = [0]
            start = 0
            for c in counts:
                # since the data is already sorted by speed
                # counts follows the data in each bin
                totalspmass = sum([m[i] for m in massgrid[start:start+c]])
                # force the summed data value as histo "weight"
                weights.append(totalspmass)
                start += c
            weights = np.array(weights)
            # finally sum the array of masses to the corresponding histo
            buckets[psp[ispecies[i]]] += weights
        else:
            continue
    return buckets, bins


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
    """(1D) Returns positions of the shock, and both inner and
    outer rayleigh line velocities joining
    both sides of the shocked cell.
    (xin, xout, cjin, cjout,
     float(ray.ds.current_time), ray.ds.parameters['x_match'])
    """
    # fields = ['sound_speed', 'density', 'pressure']
    # data, _ = reader.getLineout(fname, fields=fields,
    #                             species=False, geom=geom)
    # time, params, _, _, _ = reader.getMeta(fname)
    # rad, cs, dens, pres = data[0], data[1], data[2], data[3]

    fields = ['sound_speed', 'density', 'pressure']
    time, pars, _, _, _ = directMeta(fname)
    if len(direction) > (pars['dimensionality']-1):
        print("Direction doesn't match dimensionality: "
              "{}".format(pars['dimensionality']))
        return None
    data, _ = getLineout(fname, fields=fields, species=False,
                         direction=direction, geom=pars['geometry'])
    # time, params, _, _, _ = directMeta(fname)
    rad, cs, dens, pres = data[0], data[1], data[2], data[3]

    # this fails for x_match = y_match = z_match = 0.0
    linepos = ut.estimateMatch(direction, pars, vvv=False)
    # print(linepos)
    shockin, shockout = ut.locateShock(rad, cs, linepos, vvv=False)
    # shockin, shockout = ut.locateShock(rad, cs, pars['x_match'], vvv=False)

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
        (float tuple): specific volume, pressure, CJVelocity,
        matchhead position, time.

    """
    pos, dens, pres, gamc, cjest, pm, time = getShockConditions(fname,
                                                                addvar='gamc',
                                                                inward=inward,
                                                                **kwargs)
    # set bulk properties
    fv, fp, fg = 1.0/dens[-1], pres[-1], gamc[-1]
    av, ap, ag = 1.0/dens[0], pres[0], gamc[0]
    try:
        v, p, cj = newtonCJ(cjest, fv, fp, fg, av, ap, ag, width=width)
    except:
        # print('getNewt error')
        v, p, cj = 0.0, 0.0, cjest
    return v, p, cj, pos[1], pm, time


def newtonCJ(cjest, fuelv, fuelp, fgam, ashv, ashp, agam,
             width=0.8, fueleint=0, asheint=0):
    """fits CJ velocity to a pair of states by varying
    the speed of the rayleigh line.

    Args:
        cjest(float): starting estimate for velocity.
        fuelv(float): fuel state specific volume.
        fuelp(float): fuel state pressure.
        fgam(float): fuel state sp. heat ratio.
        ashv(float): ash state specific volume.
        ashp(float): ash state pressure.
        agam(float): ash state sp. heat ratio.
        width(float): specific volume gap for fitting.
        single(bool): skip HRlower

    Returns:
        (float tuple): specific volume, pressure, CJVelocity.

    """
    miniu = minimize(x0=cjest, tol=1e-14,
                     fun=lambda x: diffHRupper(x, ph=fuelp, vh=fuelv,
                                               gh1=fgam, gh2=agam,
                                               pr=ashp, vr=ashv,
                                               env=[ashv*(1.0-width),
                                                    ashv*(1.0+width)]))
    minil = minimize(x0=cjest, tol=1e-14,
                     fun=lambda x: diffHRlower(x, ph=fuelp, vh=fuelv,
                                               gh1=fgam, gh2=agam,
                                               pr=ashp, vr=ashv,
                                               env=[ashv*(1.0-width),
                                                    ashv*(1.0+width)]))
    cjposu, cjspdu = miniu.fun[0], miniu.x[0]
    cjposl, cjspdl = minil.fun[0], minil.x[0]
    cjpos = 0.5*(cjposu+cjposl)
    if asheint:
        cjpres = shockhugoniot_eint(cjpos, p1=fuelp, v1=fuelv,
                                    eint1=fueleint, eint2=asheint)
    else:
        cjpres = shockhugoniot(cjpos, p1=fuelp, v1=fuelv, g1=fgam, g2=agam)
    # nus are too small, so they break the rayleigh calculation
    # either that or the speed is unphysical, thus pathologic.
    if ashv - cjpos < 0:
        cjspd = rayleighSpeed(fuelp, fuelv, cjpres, cjpos)
    else:
        cjspd = rayleighSpeed(ashp, ashv, cjpres, cjpos)
#     print('{:.5E} {:.5E} {:.5E}'.format(cjspd, cjspdu, cjspdl))
    return cjpos, cjpres, cjspd


def getShockConditions(fname, inward=False, addvar='temp', direction=[]):
    """(1D) Returns bulk conditions at both sides of shock.
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
    time, pars, _, _, _ = directMeta(fname)
    if len(direction) > (pars['dimensionality']-1):
        print("Direction doesn't match dimensionality: "
              "{}".format(pars['dimensionality']))
        return None
    data, _ = getLineout(fname, fields=fields, species=False,
                         direction=direction, geom=pars['geometry'])
    # time, params, _, _, _ = directMeta(fname)
    rad, cs, dens, pres, var = data[0], data[1], data[2], data[3], data[4]

    # this fails for x_match = y_match = z_match = 0.0
    linepos = ut.estimateMatch(direction, pars, vvv=False)
    # print(linepos)
    shockin, shockout = ut.locateShock(rad, cs, linepos, vvv=False)
    # print('check shockConditions')
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
    xvar = [var[ind-offset], var[ind], var[ind+offset]]
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
    par2 = np.power(3.0*np.pi*np.pi, 2.0/3)
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
    par2 = np.power(3.0/8.0/np.pi, 1.0/3)
    par3 = np.power(Avogadro*dens*ye, 1.0/3)
    # print(par1, par2, par3)
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
    if denom < 0:
        print('Negative specific volume difference: '
              'dP {:.2e} dnu {:.2e}'.format(nom, denom))
        denom = 1.0
    return v2*np.sqrt(nom/denom)


def shockhugoniot(v2, p1=1e23, v1=0.02, g1=1.6666, g2=1.6666):
    """returns the huigoniot adiabat pressure corresponding to a
    given specific volume while passing through a set point (v1,p1)."""
    g1fac = g1*v1/(g1-1.0)
    g2fac = g2*v2/(g2-1.0)
    var = 0.5*(v1+v2)
    return p1*(g1fac-var)/(g2fac-var)


def shockhugoniot_eint(v2, p1=1e23, v1=0.02, eint1=1e17, eint2=1e18):
    """returns the huigoniot adiabat pressure corresponding to a
    given specific volume while passing through a set point (v1,p1).
    using eint
    eint2-eint2 = 1/2 (p1+p2)(v1+v2)
    """
    p2 = 2.0*(eint2-eint1)/(v1+v2)-p1
    return p2


def rayleigh(v2, p1=1e23, v1=0.02, speed=1e5):
    """returns the Rayleigh line pressure for a line crossing (v1,p1)."""
    sq = np.power(speed/v1, 2)
    fac = sq*(v1-v2)
    return p1 + fac


def diff(v, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666,
         pr=1e22, vr=0.02, speed=2.4e10):
    """yields the difference between the hugoniot adiabat and a rayleigh line,
    both passing through (v1, p1)"""
    hug = shockhugoniot(v, p1=ph, v1=vh, g1=gh1, g2=gh2)
    ray = rayleigh(v, p1=pr, v1=vr, speed=speed)
    return abs(ray) - abs(hug)


def diffHRupper(sp, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666, pr=1e22,
                vr=0.02, env=[0.04, 0.08]):
    """minimizes hugoniot adiabat - rayleigh line difference,
    starting near a strong detonation.
    (starting from the "top").
    """
    w = newton(func=lambda x: diff(x, ph=ph, vh=vh, gh1=gh1, gh2=gh2, pr=pr,
                                   vr=vr, speed=sp), x0=env[0], tol=1e-14)
    return w


def diffHRlower(sp, ph=1e23, vh=0.02, gh1=1.6666, gh2=1.6666, pr=1e22,
                vr=0.02, env=[0.04, 0.08]):
    """minimizes hugoniot adiabat - rayleigh line difference,
    starting near a weak detonation.
    (starting from the "bottom").
    """
    s = newton(func=lambda x: diff(x, ph=ph, vh=vh, gh1=gh1, gh2=gh2, pr=pr,
                                   vr=vr, speed=sp), x0=env[1], tol=1e-10)
    return s
