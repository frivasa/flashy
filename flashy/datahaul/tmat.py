from flashy.IOutils import np
import flashy.paraMan as pman
import flashy.datahaul.parData as pd
from flashy.datahaul.hdf5yt import probeDomain
_nullval = -123.4
_wedges = 5


def get2dPane(bview, file, prop, wedges=_wedges):
    """get a 2d matrix with prop values from a file,
    paralelizing the extraction of the data and 
    returning timestamp, dt and physical extents

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        file(str): relative filename.
        prop(str): UNK name to retrieve.
        wedges(int): domain cuts to distribute.

    Returns:
        (float, float, float list, np.array): time, dt, extent, data

    """
    kwargs = {'wedges': wedges,
              'fname': pman.os.path.abspath(file),
              'fields': [prop, 'x', 'y', 'dx', 'dy']}
    res = pman.throwHammer(bview, wedges, pd.par_wedge2d, **kwargs)
    d = pd.glue2dWedges(res.get())
    t, dt, widths, edges = probeDomain(file)
    ledge, redge = edges
    ext = [ledge[0], redge[0], ledge[1], redge[1]]
    nmat = buildMat(d[0], *d[1:],  widths, edges)
    return t, dt, ext, np.flip(nmat, 1).T


def buildMat(prop, x, y, dx, dy, width, edges):
    mdx = np.min(dx)
    mdy = np.min(dy)
    widths = width[:2]
    cells = [int(x) for x in widths/np.array([mdx, mdy])]
    matrix = np.full(cells, _nullval)
    xposv = x/mdx
    yposv = y/mdy
    # sum the max to get only positive values for the matrix
    if any(xposv < 0):
        xposv = xposv + max(xposv)
    if any(yposv < 0):
        yposv = yposv + max(yposv)

    xpos = [int(x) for x in xposv]
    ypos = [int(x) for x in yposv]
    # get deltas as integers
    dixv = 0.5*dx/mdx
    diyv = 0.5*dy/mdy
    dix = [int(x) for x in dixv]
    diy = [int(x) for x in diyv]
    
    for k, (i, j) in enumerate(zip(xpos, ypos)):
        # check only x since blocks are squares
        if dix[k]==0:
            matrix[i, j] = prop[k]
        else:
            matrix[i-dix[k]:i+dix[k], j-diy[k]:j+diy[k]]= prop[k]
    return matrix