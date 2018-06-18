import h5py
from flashy.nuclear import sortNuclides
from flashy.IOutils import os
from flashy.utils import np
from decimal import Decimal
_hdf5_keys = [ 'integer runtime parameters', 'integer scalars',
               'string runtime parameters', 'string scalars',
               'logical runtime parameters', 'logical scalars',
               'real runtime parameters', 'real scalars' ]

def getUNK(file, srcnames=True):
    """returns fields and species foudn in a hdf5 file.
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
    exceptions = ['n', 'd', 't' ]
    parsedvals = ['n1', 'H2', 'H3']
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
    return fields, species


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
    for p in _hdf5_keys:
        inf = finn[p].value
        # str and bool keys are flipped
        if 'string' in p:
            d.update(dict([(decode(k), decode(v)) for (v, k) in inf]))
        elif 'logical' in p:
            v, k = zip(*inf)
            k = [decode(i) for i in k]
            # switch zeroes/ones to strings
            v = [{1:'.true.', 0:'.false.'}[decode(i)] for i in v]
            d.update(dict(zip(k, v)))
        else: # floats and ints
            d.update(dict([(decode(k), decode(v)) for (k, v) in inf]))
    return d


def directMeta(file):
    """probes file, returning main properties of it.
    
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