import h5py
from flashy.nuclear import sortNuclides
from flashy.IOutils import os, getFileList
from flashy.utils import np
from decimal import Decimal
_parameter_keys = ['integer runtime parameters', 'integer scalars',
                   'string runtime parameters', 'string scalars',
                   'logical runtime parameters', 'logical scalars',
                   'real runtime parameters', 'real scalars']
_parlen = 80


def getUNK(file, srcnames=True):
    """returns fields and species found in a hdf5 file.
    (species is any unk with a number, plus (n, d, t))

    Args:
        file(str): path to file.
        srcnames(bool): return original names (else d > H2)

    Returns:
        (str list): field names.
        (str list): nuclide codes.

    """
    finn = h5py.File(file, "r")
    infval = finn['unknown names']
    unks = [decode(v[0]) for v in infval]
    fields, species = [], []
    exceptions = ['n', 'd', 't', 'd', 't']
    parsedvals = ['n1', 'H2', 'H3', 'd2', 't3']
    for field in unks:
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
    finn.close()
    return fields, species


def getBlockLineout(file):
    bulk, species = getUNK(file, srcnames=True)
    finn = h5py.File(file, 'r')
    mfracs = []
    bulks = []
    for s in species:
        data = finn[s.ljust(4)]
        mfracs.append(data.value.item())
    for b in bulk:
        if len(b) < 4:
            data = finn[b.ljust(4)]
        else:
            data = finn[b]
        bulks.append(data.value[0][0][0][0])
    pardict = dict(finn['real scalars'])
    time = pardict['time'.ljust(80).encode()]
    finn.close()
    return time, dict(zip(bulk, bulks)), species, mfracs


def getDictLikeValues(file, title, keys):
    """fetches parameter like values from a file.
    e.g.:
    irp = 'integer runtime parameters'
    irpkeys = ['lrefine_min', 'lrefine_max', 'nblockx', 'nblocky', 'nblockz']

    rrp = 'real runtime parameters'
    rrpkeys = ['xmin', 'xmax', 'ymin', 'ymax', 'zmin', 'zmax']

    isc = 'integer scalars'
    isckeys = ['nxb', 'nyb', 'nzb', 'dimensionality']

    Args:
        file(str): filename
        title(str): dataset name
        keys(str list): keys to be fetched.

    Returns:
        (float/int/str list): values for keys queried.

    """
    finn = h5py.File(file, 'r')
    data = finn[title]
    ddict = dict(data.value)
    values = []
    for k in keys:
        bink = ('{:' + str(int(_parlen)) + '}').format(k).encode()
        values.append(ddict[bink])
    return values


def getPardict(file):
    """returns every parameter found in a hdf5 file.
    (integer, string, logical, real)x(runtime parameters, scalars)

    Args:
        file(str): path to file.

    Returns:
        (dict): parameter, value dictionary.

    """
    finn = h5py.File(file, "r")
    d = {}
    for p in _parameter_keys:
#         inf = finn[p].value
        inf = finn[p][()]  # new .value syntax
        # str and bool keys are flipped
        if 'string' in p:
            d.update(dict([(decode(k), decode(v)) for (k, v) in inf]))
        elif 'logical' in p:
            k, v = zip(*inf)
            # switch zeroes/ones to strings
            v = [{1: '.true.', 0: '.false.'}[decode(i)] for i in v]
            d.update(dict([(decode(k), decode(v)) for (k, v) in inf]))
        else:  # floats and ints
            d.update(dict([(decode(k), decode(v)) for (k, v) in inf]))
    return d


def directMeta(file):
    """probes file, returning the main properties of it.

    Args:
        file(str): path to file.

    Returns:
        (float): timestamp of the file (simtime).
        (dict): parameter, value dictionary.
        (str list): field names.
        (str list): nuclide codes.
        (str list): path to folder and file itself.

    """
    pardict = getPardict(file)
    time = pardict['time']
    fields, species = getUNK(file)
    paths = []
    paths.append(os.path.dirname(os.path.abspath(file)))
    paths.append(os.path.abspath(file))
    return time, pardict, fields, species, paths


def switchGeometry(file, output, verbose=True):
    """copies hdf5 file, changing the coordinate system name to
    cartesian for yt input.

    Args:
        file(str): input filename.
        output(str): output filename.
        verbose(bool): Report file creation.

    """
    finn = h5py.File(file, "r")
    jake = h5py.File(output, "w")
    # p2 > p3: obj.iterkeys() > obj.keys()
    # p2 > p3: hdf5 works with bytes, not str: u"" > b""
    for k in finn.keys():
        finn.copy(k, jake)
    # string scalars is usually a single entry,
    # but directly changing it didn't work
    # so iterate over all values found.
    ds = jake[u'string scalars']
    newt = np.copy(ds[...])
    for i, v in enumerate(ds):
        if b"cylindrical" in v[1]:
            newt[i][1] = v[1].replace(b"cylindrical", b"cartesian  ")
    ds[...] = newt

    ds2 = jake[u'string runtime parameters']
    newt2 = np.copy(ds2[...])
    for i, v in enumerate(ds2):
        if b"cylindrical" in v[1]:
            newt2[i][1] = v[1].replace(b"cylindrical", b"cartesian  ")
    ds2[...] = newt2

    finn.close()
    jake.close()
    if verbose:
        print("Wrote {} from {}".format(output, file))


def switchback(file, output, verbose=True):
    """copies hdf5 file, reverting switchGeometry effect
    (cartesian -> cylindrical).

    Args:
        file(str): input filename.
        output(str): output filename.
        verbose(bool): Report file creation.

    """
    finn = h5py.File(file, "r")
    jake = h5py.File(output, "w")
    # p2 > p3: obj.iterkeys() > obj.keys()
    # p2 > p3: hdf5 works with bytes, not str: u"" > b""
    for k in finn.keys():
        finn.copy(k, jake)
    # string scalars is usually a single entry,
    # but directly changing it didn't work
    # so iterate over all values found.
    ds = jake[u'string scalars']
    newt = np.copy(ds[...])
    for i, v in enumerate(ds):
        if b"cartesian   " in v[1]:
            newt[i][1] = v[1].replace(b"cartesian  ", b"cylindrical")
    ds[...] = newt

    ds2 = jake[u'string runtime parameters']
    newt2 = np.copy(ds2[...])
    for i, v in enumerate(ds2):
        if b"cartesian  " in v[1]:
            newt2[i][1] = v[1].replace(b"cartesian  ", b"cylindrical")
    ds2[...] = newt2

    finn.close()
    jake.close()
    if verbose:
        print("Wrote {} from {}".format(output, file))


def buildBlockBurnFile(folders, output):
    """read all files in selected folders
    extracting single zone properties
    and writing them to a plaintext file

    Args:
        folders(str list): folders to extract data from.
        output(str): name of output file.

    Returns:
        (None): writes to file.

    """
    files = []
    for folder in folders:
        files += getFileList(folder, glob='chk', fullpath=True)
    # use first file for header.
    t, props, species, spvalues = getBlockLineout(files[0])
    head = '#time '
    head += ' '.join(list(props.keys()))
    head += ' ' + ' '.join(species)
    with open(output, 'w') as f:
        f.write(head + '\n')
        for file in files:
            t, props, species, spvalues = getBlockLineout(file)
            f.write('{:.7E} '.format(t))
            arr = np.nan_to_num(list(props.values()))
            f.write(' '.join(['{:.6E}'.format(v) for v in arr]))
            f.write(' ')
            f.write(' '.join(['{:.6E}'.format(m) for m in spvalues]))
            f.write('\n')


def turn2cartesian(folder, prefix='all', nowitness=False, silent=False):
    """Iterates over files within a folder, switching the geometry of
    hdf5 files found to cartesian.

    Args:
        folder(str): folder path.
        prefix(str): filter string (defaults to all files in the folder).
        nowitness(bool): remove non-modified files.

    """
    if prefix == 'all':
        finns = getFileList(folder)
        finns += getFileList(folder, glob='chk')
    else:
        finns = getFileList(folder)
    finns = [f for f in finns if "cart_" not in f]
    for finn in finns:
        jake = os.path.join(folder, 'cart_'+finn)
        if os.path.exists(jake):
            if not silent:
                print("{} found. Skipping.".format(jake))
            continue
        switchGeometry(os.path.join(folder, finn), jake, verbose=True)
        if nowitness:
            os.remove(os.path.join(folder, finn))


def extractVariables(source, destination, variables=['temp']):
    """creates a new hdf5 FLASH file with a reduced set of variables.

    Args:
        source(str): input filename.
        destination(str): output filename.
        variables(str list): list of named variables to extract.

    """
    finn = h5py.File(source, 'r')
    # essential mesh data for the FLASH file structure
    struct = [
        'bflags', 'block size', 'bounding box',
        'coordinates', 'gid', 'gsurr_blks',
        'integer runtime parameters',
        'integer scalars', 'logical runtime parameters',
        'logical scalars', 'node type', 'real runtime parameters',
        'real scalars', 'refine level',  'sim info',
        'string runtime parameters', 'string scalars', 'unknown names',
        'which child'
    ]
    struct += variables
    # turn selected tags to bytes
    # newunks = np.array(variables, dtype=np.bytes_)
    newunks = np.array([[unk.encode('UTF-8')] for unk in variables])
    with h5py.File(destination, 'w') as otp:
        for key in finn.keys():
            if key in struct:
                print('Wrote:', key)
                # reduce the variables to selected ones.
                if key == 'unknown names':
                    otp.create_dataset('unknown names', data=newunks)
                else:
                    otp.copy(finn[key], dest='/{}'.format(key))
    print('Wrote:', destination)
    finn.close()


def decode(entry):
    """Binary string to utf-8 decoder.
    Does not alter floats, ints, or bools.

    Args:
        entry(binary): binary string.

    Returns:
        (str): unicode string.

    """
    try:
        val = entry.strip().decode('utf-8')
        return val
    except AttributeError:
        return entry
