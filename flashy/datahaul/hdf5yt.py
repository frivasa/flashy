import yt
import os
import numpy as np
from yt.utilities.exceptions import YTFieldNotFound
from flashy.utils import get_bearing, by_mass, rot
from flashy.IOutils import log, setFolders
from flashy.nuclear import sort_nuclides, msol
import flashy.datahaul.ytfields as ytf
# avoid yt warnings
from yt.funcs import mylog
mylog.setLevel(50)
_radiusMult = 10


def getLineout(fname, fields=['density', 'temperature', 'pressure'],
               species=True, radius=5e11, geom='cartesian',
               direction=[], origin=[0.0, 0.0, 0.0], srcnames=True):
    """Returns a np.array with radius, dens, temp, pres and species specified.

    Args:
        fname(str): filename to extract data from.
        fields(str list): specify fields to extract (species are added).
        species(bool): toggle nuclide data.
        radius(float): reach of lineout.
        geom(str): geometry (['cartesian'], 'spherical').
        direction(list of float): angles of lineout.
            takes x/(x,z) as normals for 2D/3D.
            also sets dimensionality.
        origin(float list): override origin of ray to build profile.
        srcnames(bool): passthrough to hdf5yt.getFields.

    Returns:
        dblock (numpy array): matrix with fields as columns.
        species (list of str): names of species in the checkpoint.

    """
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    ds = yt.load(fname)
    # spherical (r, theta (polar), phi(azimuth))
    # cartesian (x, y, z)
    cname, bearing = get_bearing(direction, geom=geom)
#     print("calculated bearing: ",bearing)
#     ray = ds.ray([0.0, 0.0, 0.0], radius*bearing)
    ray = ds.ray(origin, radius*bearing)
    rs = np.argsort(ray['t'])
    dblock = ray[cname][rs].value

    for f in fields:
        dblock = np.vstack((dblock, ray[f][rs].value))
    _, sps = getFields(ds.field_list)
    if species:
        for s in sps:
            # force ('flash', s) to ensure it retireves the checkpoint data
            # instead of a yt variable.
            try:
                dblock = np.vstack((dblock,
                                    ray[('flash',
                                         "{}".format(s))][rs].value))
            except YTFieldNotFound:
                dblock = np.vstack((dblock,
                                    ray[('flash',
                                         "{:<4}".format(s))][rs].value))
    _, sps = getFields(ds.field_list, srcnames=srcnames)
    return dblock, sps


def probeDomain(fname):
    """Returns main domain properties from a file.

    Args:
        fname(str): filename path.

    Returns:
        (float, float, float list, float list):
            time, timestep, domain widths, domain edges

    """
    ds = yt.load(fname)
    t, dt = ds.parameters['time'], ds.parameters['dt']
    widths = ds.domain_width.value
    edges = ds.domain_left_edge.value, ds.domain_right_edge.value
    return t, dt, widths, edges


def getFields(flist, srcnames=True):
    """filters flash checkpoint field list, extracting species found
    in the checkpoint (including named exceptions, see source).

    Args:
        flist(tuple list): ds.derived_field_list from yt.
        srcnames(bool): return original names for species (n, d, t).

    Returns:
        fields (list of str): field names in checkpoint.
        species (list of str): nuclide codes as in checkpoint.

    """
    fields, species = [], []
    exceptions = ['n', 'p', 'd', 't']
    parsedvals = ['n1', 'p1', 'd2', 't3']
    for (t, field) in flist:
        if any(char.isdigit() for char in field):
            species.append(field)
        elif field.strip() in exceptions:
            species.append(field)
        else:
            fields.append(field)
    species = sort_nuclides(species)
    # once sorted, rename the exceptions
    if srcnames:
        for i, e in enumerate(parsedvals):
            if e in species:
                species[species.index(e)] = exceptions[i]
    return fields, species


def getMeta(fname, print_stats=False, showUNK=False):
    """returns metadata from checkpoint and yt-dereived fields.

    Args:
        fname(str): filename to check.
        print_stats(bool): print ds.print_stats (grid info).
        showUNK(bool): print ds.field_list (unk variables in chk).

    Returns:
        (float): simtime of checkpoint.
        (dict): flash.par dictionary.
        (str list): fields in checkpoint.
        (str list): species in checkpoint.
        (str list): fullpath, filename of checkpoint.

    """
    ds = yt.load(fname)
    if print_stats:
        ds.print_stats()
    if showUNK:
        print(ds.field_list)
    fields, species = getFields(ds.derived_field_list)
    filepaths = ds.fullpath, ds.parameter_filename
    return float(ds.current_time), ds.parameters, fields, species, filepaths


def getExtrema(fname, flist=['density', 'temperature', 'pressure']):
    """returns a list of tuples with extrema for given fields(flist).

    Args:
        fname(str): filename to probe.
        flist(list of str): fields to query.

    Returns:
        (list of np.arrays): list of paired extrema [min max]
        for each field queried.

    """
    ds = yt.load(fname)
    ad = ds.all_data()
    if len(flist) == 1:
        return [ad.quantities.extrema(flist).value]
    else:
        return [x.value for x in ad.quantities.extrema(flist)]


def getYields(fname, subsetR=0.0):
    """returns time, summed masses and species in a flash otp file.
    all units are msun or cgs.
    1D assumes spherical geometry
    2D assumes cylindrical geometry and trims cart_ from filename.

    Args:
        fname(str): filename to inspect.
        subsetR(float): radius cutoff for extraction.

    Returns:
        time (float) [s],
        species names: (list of str),
        masses: (float list) [msun]

    """
    if 'cart' in fname:
        stem, bud = os.path.split(fname)
        filename = os.path.join(stem, bud[5:])
    else:
        filename = fname
    ds = yt.load(filename)
    if subsetR:
        log.warning('R subset for yields: {:.2E}'.format(subsetR))
        ad = ds.sphere([0.0, 0.0, 0.0], (subsetR, 'cm'))
    else:
        ad = ds.all_data()
    _, species = getFields(ds.field_list)
    species = [s.ljust(4) for s in species]
    # 2d glitches due to dz being nonsensical
    if ds.parameters['dimensionality'] == 1:
        try:
            dx = ad['path_element_x'].value
            x = ad['x'].value
        except YTFieldNotFound:
            dx = ad['path_element_r'].value
            x = ad['r'].value
        fac = 4.0/3.0*np.pi
        if len(x) > 1:
            rdiff = np.diff(x)
            rdiff = np.append(rdiff, rdiff[-1])
            skewed = 0.5*rdiff + x
            skewed = np.insert(skewed, 0, 0.0)
            rcub = skewed**3
            vols = fac*np.diff(rcub)
        else:
            vols = np.array([fac*ad['dr'].value])
        cell_masses = vols*ad['density'].value/msol
        msunMs = []
        for sp in species:
            msunMs.append(np.sum(cell_masses*ad[sp].value))
        log.warning('Assumed spherical cells for volume')
    elif ds.parameters['dimensionality'] == 2:
        log.warning('Cylindrical volume cells')
        masses = ad.quantities.weighted_average_quantity(species, 'cell_mass')
        total = ad.quantities.total_mass()
        msunMs = [m.value*total[0].value/msol for m in masses]
    else:
        masses = ad.quantities.weighted_average_quantity(species, 'cell_mass')
        total = ad.quantities.total_mass()
        msunMs = [m.value*total[0].value/msol for m in masses]
    return ds.current_time.value, species, msunMs


def wedge2d(fname, radius, sphereCapR=1e9, reference='y',
            fields=[], cylvel=False):
    """same as wedge3d but for a cylindrical 2d domain
    Assumes cylindrical geometry to calculate accurate cell_mass
    default fields: vx, vy, vz, cell_mass and species.
    """
    ds = yt.load(fname)
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    domx, domy, _ = ds.domain_width.value
    # center is midway to edge of domain
    # normal is pointing in direction of axis
    center, normal, height = {
        'x': [(0.5*domx, 0.0, 0.0), (domx, 0.0, 0.0), 0.5*domx],
        'y': [(0.0, 0.5*domy, 0.0), (0.0, domy, 0.0), 0.5*domy]
    }[reference]
    cylinder = ds.disk(center=ds.arr(center, "code_length"),
                       normal=ds.arr(normal, "code_length"),
                       radius=(radius, "code_length"),
                       height=(height, "code_length"))
    capDistance = ds.sphere([0, 0, 0], sphereCapR)
    borehole = ds.intersection([cylinder, capDistance])
    # get fields and data from data file
    if not fields:
        _, species = getFields(ds.field_list)
        # calculate cylindrical volumes (dz and z are non-sensical in 2D)
        dx = borehole['path_element_x'].value
        dy = borehole['path_element_y'].value
        r = borehole['x'].value
        cylvol = 2.0*np.pi*dy*dx*r
        cell_masses = cylvol*borehole['density'].value
        if cylvel:
            # asking for "away from symmetry axis" velocity
            speeds = borehole['velx'].value
        else:
            # 2d, so only use x and y
#             vx = borehole['velx'].value
#             vy = borehole['velx'].value
#             speeds = np.sqrt(np.power(vx,2) + np.power(vy,2))
            speeds = borehole['velocity_magnitude'].value
        rawd = []
#         fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
#         for f in fields:
#             if len(f) < 4 and len(f) > 1:  # 'o16 ', 'n14 ', etc
#                 rawd.append(borehole[f.ljust(4)].value)
#             else:
#                 rawd.append(borehole[f].value)
        rawd.append(borehole['velx'].value)
        rawd.append(borehole['vely'].value)
        rawd.append(borehole['velz'].value)
        rawd.append(cell_masses)
        for s in species:
            rawd.append(borehole[s.ljust(4)].value)
        log.warning('Borehole modified (cylindrical masses)')
        return rawd, species  # WARNING: this might be humongous.
    else:
        rawd = []
        for f in fields:
            rawd.append(borehole[f].value)
        log.info('borehole Unmodified')
        # WARN: this assumes same amount of fields always
        # vel xyz and cell_mass
        log.warning('Assuming species start at 5th slot')
        species = fields[4:]
        return rawd, species


def wedge3d(chkp, radius, sphereCapR=1e9, reference='z', fields=[]):
    """cut a cylinder from the origin towards an axis 
    to perform velocity vs mass fraction measurements.
    default fields: vx, vy, vz, cell_mass and species

    Args:
        fname(str): file name
        radius(float): radius of cylinder (half diam of borehole)
        sphereCapR(float): maximum radius to consider (sphere).
        reference(str): pointing direction of cylinder and 
            which velocity to use.
        fields(str list): override default fields.

    Returns:
        (tuple of np.arrays): raw data sorted by field.

    """
    ds = yt.load(chkp)
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    domx, domy, domz = ds.domain_width.value
    # center is midway to edge of domain
    # normal is pointing in direction of axis
    center, normal, height = {
        'x': [(0.5*domx, 0.0, 0.0), (domx, 0.0, 0.0), 0.5*domx],
        'y': [(0.0, 0.5*domy, 0.0), (0.0, domy, 0.0), 0.5*domy],
        'z': [(0.0, 0.0, 0.5*domz), (0.0, 0.0, domz), 0.5*domz]
    }[reference]
    
    cylinder = ds.disk(center=ds.arr(center, "code_length"),
                       normal=ds.arr(normal, "code_length"),
                       radius=(radius, "code_length"),
                       height=(height, "code_length"))
    capDistance = ds.sphere([0, 0, 0], sphereCapR)
    print('pre intersection')
    borehole = ds.intersection([cylinder, capDistance])
    # get fields and data from data file
    if not fields:
        print('extracting data. (no fields specified)')
        _, species = getFields(ds.field_list)
        species = [x.ljust(4) for x in species]
        fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
        print('fields to extract:', len(fields))
        # rawd 0 1 2 3 fixed.
        # offset = 4
        # ask yt for data, as always this takes forever.
        rawd = []
        for i, f in enumerate(fields):
            if not i: 
                print(borehole[f].size)
            rawd.append(borehole[f].value)
        return rawd, species  # WARNING: this might be humongous.
    else:
        print('extracting data. fields:', fields)
        rawd = []
        for f in fields:
            rawd.append(borehole[f].value)
        # WARN: this assumes same amount of fields always
        # vel xyz and cell_mass
        log.warning('Assuming species start at 5th slot')
        species = fields[4:]
        return rawd, species


def wedge2d_disks(fname, elevation, depth, fields=[], cylvel=False):
    """cut a wedge in a 2d rectangular domain to perform velocity
    vs mass fraction measurements.
    
    Assumes cylindrical geometry to calculate accurate cell_mass
    
    default fields: vx, vy, vz, cell_mass and species

    Args:
        fname(str): file name
        elevation(float): equator-north pole wedge angle.
        depth(float): equator-south pole wedge angle.
        fields(str list): override default fields.
        cylvel(bool): return 'velx' for speed data.

    Returns:
        (tuple of np.arrays): raw data, [species list].

    """
    ds = yt.load(fname)
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    domx, domy, _ = ds.domain_width.value

    # top cylinder pivoting on the top right of the domain
    # this sets the equator-south pole angle
    alt = np.deg2rad(depth)

    diag = np.sqrt(0.25*domy*domy + domx*domx)
    diagang = np.arctan(0.5*domy/domx)

    alpha = 0.5*np.pi-alt
    N = domx/np.cos(alpha)
    delta = domx/np.cos(alpha) - 0.5*domy/np.sin(alpha)
    h = np.tan(alpha)*domx - 0.5*domy
    normal = [-domx, -0.5*domy - h, 0.0]

    if depth > 0.0:
        eps = h * np.cos(alpha)*np.cos(alpha)/np.sin(alpha)
        height = N - delta+eps
        center = [domx, 0.5*domy, 0.0]
    else:
        height = np.cos(alt)*0.5*domy
        center = [0.0, 0.5*domy, 0.0]

    cylinder1 = ds.disk(center=ds.arr(center, "code_length"),
                        normal=ds.arr(normal, "code_length"),
                        # absurd radius to grab all cells.
                        radius=(_radiusMult*domy, "code_length"),
                        height=(height, "code_length"))

    # bottom cylinder pivoting on the bottom left of the domain
    # this sets the equator-north pole angle
    alt = np.deg2rad(elevation)

    diag = np.sqrt(0.25*domy*domy + domx*domx)
    diagang = np.arctan(0.5*domy/domx)

    alpha = 0.5*np.pi-alt
    N = domx/np.cos(alpha)
    delta = domx/np.cos(alpha) - 0.5*domy/np.sin(alpha)
    h = np.tan(alpha)*domx - 0.5*domy
    normal = [-domx, 0.5*domy + h, 0.0]

    if elevation > 0.0:
        eps = h * np.cos(alpha)*np.cos(alpha)/np.sin(alpha)
        height = N - delta+eps
        center = [domx, -0.5*domy, 0.0]
    else:
        height = np.cos(alt)*0.5*domy
        center = [0.0, -0.5*domy, 0.0]

    cylinder2 = ds.disk(center=ds.arr(center, "code_length"),
                        normal=ds.arr(normal, "code_length"),
                        # absurd radius to grab all cells.
                        radius=(_radiusMult*domy, "code_length"),
                        height=(height, "code_length"))
    wedge = ds.intersection([cylinder1, cylinder2])

    # get fields and data from data file
    if not fields:
        _, species = getFields(ds.field_list)
        # calculate cylindrical volumes (dz and z are non-sensical in 2D)
        dx = wedge['path_element_x'].value
        dy = wedge['path_element_y'].value
        r = wedge['x'].value
        cylvol = 2.0*np.pi*dy*dx*r
        cell_masses = cylvol*wedge['density'].value
        if cylvel:
            # asking for "away from symmetry axis" velocity
            speeds = wedge['velx'].value
        else:
            # 2d, so only use x and y
#             vx = wedge['velx'].value
#             vy = wedge['velx'].value
#             speeds = np.sqrt(np.power(vx,2) + np.power(vy,2))
            speeds = wedge['velocity_magnitude'].value
        rawd = []
#         fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
#         for f in fields:
#             if len(f) < 4 and len(f) > 1:  # 'o16 ', 'n14 ', etc
#                 rawd.append(wedge[f.ljust(4)].value)
#             else:
#                 rawd.append(wedge[f].value)
        rawd.append(speeds)
        rawd.append(cell_masses)
        for s in species:
            rawd.append(wedge[s.ljust(4)].value)
        log.warning('Wedge modified (cylindrical masses)')
        return rawd, species  # WARNING: this might be humongous.
    else:
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        log.info('Wedge Unmodified')
        return rawd


def wedge3d_disks(chkp, elevation, depth, reference='x', fields=[], 
                  antipode=False, maxR=1e9):
    """cut a wedge in a 3d rectangular domain to perform velocity
    vs mass fraction measurements.
    default fields: vx, vy, vz, cell_mass and species

    Args:
        fname(str): file name
        elevation(float): equator-north pole wedge angle.
        depth(float): equator-south pole wedge angle.
        reference(str): equinox reference for the equator.
        fields(str list): override default fields.
        antipode(bool): invert selection, antipodal wedge.

    Returns:
        (tuple of np.arrays): raw data sorted by field.

    """
    ds = yt.load(chkp)
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    domx, domy, domz = ds.domain_width.value
    rad = 100*np.max(ds.domain_width.value)  # huge radius to grab all cells.
    inv = -1.0 if antipode else 1.0  # flip selection

    tvec, bvec = {
        'x': [(0.0, 0.5*domy, 0.5*domz), (0.0, 0.5*domy, -0.5*domz)],
        'y': [(0.5*domx, 0.0, 0.5*domz), (-0.5*domx, 0.0, 0.5*domz)],
        'z': [(0.5*domx, 0.5*domy, 0.0), (0.5*domx, -0.5*domy, 0.0)]
    }[reference]
    topv = np.array(tvec)
    botv = np.array(bvec)

    if depth*elevation < 0.0:
        # in this case we are on a quadrant
        if depth < 0.0:
            # upper right, use 2 bottom cylinders:
            # normal bot cylinder up to abs(depth)
            equ = 45 - abs(depth)
            alt = np.deg2rad(equ)
            rm = rot(-alt, reference)
            normal = rm.dot(botv)
            height = np.cos(alt)*np.sqrt(botv.dot(botv))
            center = botv*inv
            cylinder1 = ds.disk(center=ds.arr(center, "code_length"),
                                normal=ds.arr(normal, "code_length"),
                                radius=(_radiusMult*domy, "code_length"),
                                height=(height, "code_length"))
            # second cylinder up to elevation but flipped to grab the wedge
            equ = 45 - elevation
            alt = np.deg2rad(equ)
            rm = rot(-alt, reference)
            normal = rm.dot(botv)
            height = np.cos(alt)*np.sqrt(botv.dot(botv))
            center = -botv*inv
            cylinder2 = ds.disk(center=ds.arr(center, "code_length"),
                                normal=ds.arr(normal, "code_length"),
                                radius=(_radiusMult*domy, "code_length"),
                                height=(height, "code_length"))
        else:
            # lower right, use 2 top cylinders:
            # first top goes down to depth and in inverted
            equ = 45 - depth  # 45 degree offset from center-origin reference
            dep = np.deg2rad(equ)
            rm = rot(dep, reference)
            normal = rm.dot(topv)
            height = np.cos(dep)*np.sqrt(topv.dot(topv))
            center = topv*inv
            cylinder1 = ds.disk(center=ds.arr(center, "code_length"),
                                normal=ds.arr(normal, "code_length"),
                                radius=(rad, "code_length"),
                                height=(height, "code_length"))
            # second top cylinder goes down to abs(elevation)
            # 45 degree offset from center-origin reference
            equ = 45 - abs(elevation)
            dep = np.deg2rad(equ)
            rm = rot(dep, reference)
            normal = rm.dot(topv)
            height = np.cos(dep)*np.sqrt(topv.dot(topv))
            center = -topv*inv
            cylinder2 = ds.disk(center=ds.arr(center, "code_length"),
                                normal=ds.arr(normal, "code_length"),
                                radius=(rad, "code_length"),
                                height=(height, "code_length"))
    else:
        # default wedge covering the equator
        # top cylinder
        equ = 45 - depth  # 45 degree offset from center-origin reference
        dep = np.deg2rad(equ)
        rm = rot(dep, reference)
        normal = rm.dot(topv)
        height = np.cos(dep)*np.sqrt(topv.dot(topv))
        center = topv*inv
        cylinder1 = ds.disk(center=ds.arr(center, "code_length"),
                            normal=ds.arr(normal, "code_length"),
                            radius=(rad, "code_length"),
                            height=(height, "code_length"))
        # bottom cylinder
        equ = 45 - elevation  # 45 degree offset from center-origin reference
        alt = np.deg2rad(equ)
        rm = rot(-alt, reference)

        normal = rm.dot(botv)
        height = np.cos(alt)*np.sqrt(botv.dot(botv))
        center = botv*inv
        cylinder2 = ds.disk(center=ds.arr(center, "code_length"),
                            normal=ds.arr(normal, "code_length"),
                            radius=(rad, "code_length"),
                            height=(height, "code_length"))
    print('pre intersection')
    capDistance = ds.sphere([0, 0, 0], maxR)
    wedge = ds.intersection([cylinder1, cylinder2, capDistance])
    # get fields and data from data file
    if not fields:
        print('extracting data. (no fields specified)')
        _, species = getFields(ds.field_list)
        species = [x.ljust(4) for x in species]
        fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
        print(len(fields))
        # rawd 0 1 2 3 fixed.
        # offset = 4
        # ask yt for data, as always this takes forever.
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        return rawd, species  # WARNING: this might be humongous.
    else:
        print('extracting data. fields:', fields)
        rawd = []
        # WARN: this assumes same amount of fields always
        # vel xyz and cell_mass
        for f in fields:
            rawd.append(wedge[f].value)
        species = fields[4:]
        return rawd, species


def getPointData(file, posx, posy):
    """get all unk data (plus speed, soundspeed and fermiDeg) from a
    single cell in a file to track.
    otpn: pointTrack

    # get cell masses
    # log.warning('Assuming cylindrical slice')
    # dx = data['path_element_x']
    # dy = data['path_element_y']
    # r = data['x']
    # cylvol = 2.0*np.pi*dy*dx*r
    # cell_masses = cylvol*data['dens']
    # print('cell mass (calc <> data)', cell_masses, data['cell_mass'])
    # print('cell vol (calc <> data)', cylvol, data['cell_volume'])

    Args:
        file(str): filepath.
        posx(float): x-axis position.
        posy(float): y-axis position.

    """
    ds = yt.load(file)
    # read in custom fields
    cfl = ['speed', 'sound_speed', 'fermiDeg']
    for f in cfl:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    ad = ds.all_data()
    cellsize = np.min(ad['dx'].v)
    sp = ds.sphere([posx, posy, 0.0], cellsize)
    flds, sps = getFields(ds.field_list)
    sps = [s.ljust(4) for s in sps]
    qu = sp.quantities
    ffl = ['x', 'y', 'cell_volume', 'cell_mass', 'sound_speed',
           'courant_time_step', 'dynamical_time',
           'path_element_x', 'path_element_y']
    allf = ffl + flds + sps
    # there's only one cell so min/max is irrelevant, nor the field queried
    vals = qu.sample_at_max_field_values('dens', sample_fields=allf)
    vals = [x.v for x in vals[1:]]  # drop the first value (repeated)
    # write all data from the cell to file
    # use radius as an alias for time so that datamatrix works
    # change x and y to posx posy again so that dm wont take them as n or p.
    allf[0] = 'posx'
    allf[1] = 'posy'
    metatxt = ['# radius {}'.format(' '.join(allf))]
    parsedvals = ['{:.8E}'.format(float(ds.current_time))]
    parsedvals += ['{:.8E}'.format(x) for x in vals]
    metatxt.append(' '.join(parsedvals))

    dest, num, name = setFolders(file, 'pointTrack')
    with open(name + '.meta', 'w') as f:
        f.write('\n'.join(metatxt))
    print("Wrote: {}".format(name + '.meta'))


def getPointInfo(file, posx, posy, fields=[], comp=False):
    """get a list of fields and composition from a point 
    in a 2D checkpoint.
    
    Args:
        file(str): filepath.
        posx(float): x-axis position.
        posy(float): y-axis position.
        fields(str list): specify fields to get (overrides output).
        comp(bool): add composition to field list.

    """
    ds = yt.load(file)
    ad = ds.all_data()
    cellsize = np.min(ad['dx'].v)
    sp = ds.sphere([posx, posy, 0.0], cellsize)
    flds, sps = getFields(ds.field_list)
    sps = [s.ljust(4) for s in sps]
    qu = sp.quantities
    if not fields:
        comp = True
    if comp:
        allf = fields + sps
    else:
        allf = fields
    # there's only one cell so min/max is irrelevant, nor the field queried
    vals = qu.sample_at_max_field_values('dens', sample_fields=allf)
    vals = [x.v for x in vals[1:]]  # drop the first value (dens)
    # write all data from the cell to file
    # use radius as an alias for time so that datamatrix works
    # change x and y to posx posy again so that dm wont take them as n or p.
    return dict(zip(allf, vals))


def getPointInfo3D(file, pos, fields=[], comp=False):
    """get a list of fields and composition from a point 
    in a 3D checkpoint.
    
    Args:
        file(str): filepath.
        pos(float list): position triad.
        fields(str list): specify fields to get (overrides output).
        comp(bool): add composition to field list.

    """
    ds = yt.load(file)
    # ad = ds.all_data()  # this is huge, don't do it
    # cellsize = np.min(ad['dx'].v)
    # do a 100km sphere and take the smallest dx
    sp = ds.sphere(pos, 100e5)
    cellsize = sp.min('dx')
    sp = ds.sphere(pos, cellsize)
    flds, sps = getFields(ds.field_list)
    sps = [s.ljust(4) for s in sps]
    qu = sp.quantities
    if not fields:
        comp = True
    if comp:
        allf = fields + sps
    else:
        allf = fields
    # there's only one cell so min/max is irrelevant, nor the field queried
    vals = qu.sample_at_max_field_values('dens', sample_fields=allf)
    vals = [x.v for x in vals[1:]]  # drop the first value (dens)
    # write all data from the cell to file
    # use radius as an alias for time so that datamatrix works
    # change x and y to posx posy again so that dm wont take them as n or p.
    return dict(zip(allf, vals))


def bark3d(chkp, r0, r1, center, left_edge, right_edge, fields=[]):
    """cut a box intersecting two spheres, forming a rounded slab.
    default fields: vx, vy, vz, cell_mass and species

    Args:
        fname(str): file name
        r0, r1 (float): slab radii.
        center (float list): box center.
        left_edge, right_edge (float list): box edges.
        fields(str list): override default fields.
        
    Returns:
        (tuple of np.arrays): raw data sorted by field.

    """
    ds = yt.load(chkp)
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)

    # box = ds.region(center, left_edge, right_edge)
    # print("using box, not region")
    box = ds.box(left_edge, right_edge)
    sp0 = ds.sphere([0.0, 0.0, 0.0], r0)
    sp1 = ds.sphere([0.0, 0.0, 0.0], r1)
    # cut top curve from box
    outer_curve = sp1 & box
    # cut bottom curve from box
    bark = outer_curve - sp0
    # bark = ds.intersection([box, sp0, sp1])
    # get fields and data from data file
    if not fields:
        _, species = getFields(ds.field_list)
        fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
        # rawd 0 1 2 3 fixed.
        # offset = 4
        # ask yt for data, as always this takes forever.
        rawd = []
        for f in fields:
            rawd.append(bark[f].value)
        return rawd, species  # WARNING: this might be humongous.
    else:
        rawd = []
        for f in fields:
            rawd.append(bark[f].value)
        return rawd
