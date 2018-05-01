import yt
from flashy.utils import np, getBearing
from flashy.nuclear import sortNuclides
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
               radius=5e11, geom='cartesian', direction=[]):
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
    cname, bearing = getBearing(direction, geom)
    ray = ds.ray([0.0, 0.0, 0.0], radius*bearing)
    rs = np.argsort(ray['t'])

    dblock = ray[cname][rs].value
    for f in fields:
        dblock = np.vstack((dblock, ray[f][rs].value))
    _, sps = getFields(ds.field_list)
    if species:
        for s in sps:
            dblock = np.vstack((dblock, ray[s][rs].value))
    return dblock, sps


def getFields(flist):
    """filters flash checkpoint field list, extracting species found
    in the checkpoint (including named exceptions, see source).
    
    Args:
        flist (list of tuples): ds.derived_field_list from yt.
    
    Returns:
        fields (list of str): field names in checkpoint.
        species (list of str): nuclide codes as in checkpoint.
    
    """
    fields, species = [], []
    exceptions = ['n', 'd', 't', 'neut']
    for (t, field) in flist:
        if any(char.isdigit() for char in field):
            species.append(field)
        elif field.strip() in exceptions:
            species.append(field)
        else:
            fields.append(field)
    species = sortNuclides(species)
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