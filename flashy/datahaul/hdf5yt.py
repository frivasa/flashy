import yt
from flashy.utils import np, getBearing, byMass, rot
from flashy.nuclear import sortNuclides, msol
from scipy.integrate import trapz
# avoid yt warnings
from yt.funcs import mylog
mylog.setLevel(50)

# custom fields
def _speed(field, data):
    vx = data['flash', 'velx']
    vy = data['flash', 'vely']
    spd = np.sqrt(vx*vx + vy*vy)
    return spd
#yt.add_field(("flash","speed"), function=_speed, units="cm/s", take_log=False)


def getLineout(fname, fields=['density', 'temperature', 'pressure'], species=True,
               radius=5e11, geom='cartesian', direction=[], srcnames=True):
    """Returns a np.array with radius, dens, temp, pres and species specified.
    
    Args:
        fname (str): filename to extract data from.
        species (bool): toggle nuclide data.
        radius (float): reach of lineout.
        geom (str): geometry (['cartesian'], 'spherical').
        direction (list of float): angles of lineout.
            takes x/(x,z) as normals for 2D/3D. 
            also sets dimensionality.
    
    Returns:
        dblock (numpy array): matrix with fields as columns.
        species (list of str): names of species in the checkpoint.
    
    """
    ds = yt.load(fname)
    # spherical (r, theta, phi)
    # cartesian (x, y, z)
    cname, bearing = getBearing(direction, geom=geom)
    ray = ds.ray([0.0, 0.0, 0.0], radius*bearing)
    rs = np.argsort(ray['t'])
    dblock = ray[cname][rs].value
    
    for f in fields:
        dblock = np.vstack((dblock, ray[f][rs].value))
    _, sps = getFields(ds.field_list)
    if species:
        for s in sps:
            dblock = np.vstack((dblock, ray[s][rs].value))
    _, sps = getFields(ds.field_list, srcnames=srcnames)
    return dblock, sps


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
    exceptions = ['n', 'd', 't' ]
    parsedvals = ['n1', 'd2', 't3']
    for (t, field) in flist:
        if any(char.isdigit() for char in field):
            species.append(field)
        elif field.strip() in exceptions:
            species.append(field)
        else:
            fields.append(field)
    species = sortNuclides(species)
    # once sorted, rename the exceptions
    if srcnames:
        for i, e in enumerate(parsedvals):
            if e in species:
                species[species.index(e)] = exceptions[i]
    return fields, species


def getMeta(fname, print_stats=False):
    """returns metadata from checkpoint
    
    Args:
        fname(str): filename to check.
    
    Returns:
        (float): simtime of checkpoint.
        (dict): flash.par dictionary.
        (list of str): fields in checkpoint.
        (list of str): species in checkpoint.
        (list of str): fullpath, filename of checkpoint.
    
    """
    ds = yt.load(fname)
    if print_stats:
        ds.print_stats()
    fields, species = getFields(ds.derived_field_list)
    filepaths = ds.fullpath, ds.parameter_filename
    return float(ds.current_time), ds.parameters, fields, species, filepaths


def getExtrema(fname, flist=['density', 'temperature', 'pressure']):
    """returns a list of tuples with extrema for given fields(flist).
    
    Args:
        fname(str): filename to probe.
        flist(list of str): fields to query.
    
    Returns:
        (list of np.arrays): list of paired extrema [min max] for each field queried.
    
    """
    ds = yt.load(fname)
    ad = ds.all_data()
    if len(flist)==1:
        return [ad.quantities.extrema(flist).value]
    else:
        return [x.value for x in ad.quantities.extrema(flist)]


def getYields(fname):
    """returns time, summed masses and species in a flash otp file.
    all units are msun or cgs.
    
    Args:
        fname(str): filename to inspect.
    
    Returns:
        time (float) [s], 
        species names: (list of str),
        masses: (float list) [msun]

    """
    ds = yt.load(fname)
    ad = ds.all_data()
    _, species = getFields(ds.field_list)
    masses = ad.quantities.weighted_average_quantity(species, 'cell_mass')
    total = ad.quantities.total_mass()
    msunMs = [m.value*total[0].value/msol for m in masses]
    
    return ds.current_time.value, species, msunMs
    

def wedge2d(fname, elevation, depth, fields=[]):
    """cut a wedge in a 2d rectangular domain to perform velocity 
    vs mass fraction measurements.
    
    Args:
        fname(str): file name
        elevation(float): equator-north pole wedge angle.
        depth(float): equator-south pole wedge angle.
        fields(str list): override default fields.

    Returns:
        (tuple of np.arrays): raw data, [species list].
    
    """
    ds = yt.load(fname)
    domx, domy, _ = ds.domain_width.value
    
    # top cylinder pivoting on the top right of the domain
    # this sets the equator-south pole angle
    alt = np.deg2rad(depth)

    diag = np.sqrt(0.25*domy*domy + domx*domx)
    diagang = np.arctan(0.5*domy/domx)

    alpha = 0.5*np.pi-alt
    N = domx/np.cos(alpha)
    delta =  domx/np.cos(alpha) - 0.5*domy/np.sin(alpha)
    h = np.tan(alpha)*domx - 0.5*domy
    normal = [-domx, -0.5*domy - h, 0.0]

    if depth>0.0:
        eps = h * np.cos(alpha)*np.cos(alpha)/np.sin(alpha)
        height = N - delta+eps
        center = [domx, 0.5*domy, 0.0]
    else:
        height = np.cos(alt)*0.5*domy
        center = [0.0, 0.5*domy, 0.0]        

    cylinder1 = ds.disk(center=ds.arr(center, "code_length"), 
                        normal=ds.arr(normal, "code_length"),
                        radius=(10*domy, "code_length"),  # absurd radius to grab all cells.
                        height=(height, "code_length"))

    # bottom cylinder pivoting on the bottom left of the domain
    # this sets the equator-north pole angle
    alt = np.deg2rad(elevation)

    diag = np.sqrt(0.25*domy*domy + domx*domx)
    diagang = np.arctan(0.5*domy/domx)

    alpha = 0.5*np.pi-alt
    N = domx/np.cos(alpha)
    delta =  domx/np.cos(alpha) - 0.5*domy/np.sin(alpha)
    h = np.tan(alpha)*domx - 0.5*domy
    normal = [-domx, 0.5*domy + h, 0.0]

    if elevation>0.0:
        eps = h * np.cos(alpha)*np.cos(alpha)/np.sin(alpha)
        height = N - delta+eps
        center = [domx, -0.5*domy, 0.0]
    else:
        height = np.cos(alt)*0.5*domy
        center = [0.0, -0.5*domy, 0.0]

    cylinder2 = ds.disk(center=ds.arr(center, "code_length"),
                        normal=ds.arr(normal, "code_length"),
                        radius=(10*domy, "code_length"),  # absurd radius to grab all cells.
                        height=(height, "code_length"))
    wedge = ds.intersection([cylinder1, cylinder2])
    
    # get fields and data from data file
    if not fields:
        _, species = getFields(ds.field_list)
        fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
        # rawd 0 1 2 3 fixed.
        # offset = 4
        # ask yt for data, as always this takes forever.
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        return rawd, species  # WARNING: this might be humongous. 
    else:
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        return rawd


def wedge3d(chkp, elevation, depth, reference='x', fields=[], antipode=False):
    """cut a wedge in a 3d rectangular domain to perform velocity 
    vs mass fraction measurements.
    
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

    if depth*elevation<0.0:
        # in this case we are on a quadrant
        if depth<0.0:
            # upper right, use 2 bottom cylinders:
            # normal bot cylinder up to abs(depth)
            equ = 45 - abs(depth)
            alt = np.deg2rad(equ)
            rm = rot(-alt, reference)
            normal = rm.dot(botv)
            height = np.cos(alt)*np.sqrt(botv.dot(botv))
            center = botv*inv
            cylinder1 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                                radius=(10*domy, "code_length"), height=(height, "code_length"))
            # second cylinder up to elevation but flipped to grab the wedge
            equ = 45 - elevation
            alt = np.deg2rad(equ)
            rm = rot(-alt, reference)
            normal = rm.dot(botv)
            height = np.cos(alt)*np.sqrt(botv.dot(botv))
            center = -botv*inv
            cylinder2 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                                radius=(10*domy, "code_length"), height=(height, "code_length"))
        else:
            # lower right, use 2 top cylinders:
            # first top goes down to depth and in inverted
            equ = 45 - depth  # 45 degree offset from center-origin reference
            dep = np.deg2rad(equ)
            rm = rot(dep, reference)
            normal = rm.dot(topv)
            height = np.cos(dep)*np.sqrt(topv.dot(topv))
            center = topv*inv
            cylinder1 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                                radius=(rad, "code_length"), height=(height, "code_length"))
            # second top cylinder goes down to abs(elevation)
            equ = 45 - abs(elevation)  # 45 degree offset from center-origin reference
            dep = np.deg2rad(equ)
            rm = rot(dep, reference)
            normal = rm.dot(topv)
            height = np.cos(dep)*np.sqrt(topv.dot(topv))
            center = -topv*inv
            cylinder2 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                                radius=(rad, "code_length"), height=(height, "code_length"))
    else:
        # default wedge covering the equator
        # top cylinder
        equ = 45 - depth  # 45 degree offset from center-origin reference
        dep = np.deg2rad(equ)
        rm = rot(dep, reference)
        normal = rm.dot(topv)
        height = np.cos(dep)*np.sqrt(topv.dot(topv))
        center = topv*inv
        cylinder1 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                            radius=(rad, "code_length"), height=(height, "code_length"))
        # bottom cylinder
        equ = 45 - elevation  # 45 degree offset from center-origin reference
        alt = np.deg2rad(equ)
        rm = rot(-alt, reference)

        normal = rm.dot(botv)
        height = np.cos(alt)*np.sqrt(botv.dot(botv))
        center = botv*inv
        cylinder2 = ds.disk(center=ds.arr(center, "code_length"), normal=ds.arr(normal, "code_length"),
                            radius=(rad, "code_length"), height=(height, "code_length"))

    wedge = ds.intersection([cylinder1, cylinder2])
    #XXX
    print('cylinders created, getting data')
    # get fields and data from data file
    if not fields:
        _, species = getFields(ds.field_list)
        fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
        # rawd 0 1 2 3 fixed.
        # offset = 4
        # ask yt for data, as always this takes forever.
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        return rawd, species  # WARNING: this might be humongous. 
    else:
        rawd = []
        for f in fields:
            rawd.append(wedge[f].value)
        return rawd
