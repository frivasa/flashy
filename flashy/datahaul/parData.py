from flashy.datahaul.hdf5yt import wedge2d
from flashy.IOutils import np


def par_wedge2d(wedgenum, wedges=5, fname='', fields=[]):
    """parallelizable wrapper for 2d wedges."""
    delta = 180.0/wedges
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
    return wedge2d(fname, start, stop, fields)


def glue2dWedges(pload):
    """rejoin data after extraction."""
    fields = len(pload[0])
    data = []
    for i in range(fields):
        field = np.array([])
        for wedge in pload:
            field = np.append(field, wedge[i])
        data.append(field)
    return data
