"""module for handling plain text 1D profiles sorted in columns and with a typical header

There's some plotting routines here that should be taken elsewhere, plus headers are too 
strict, for example, the main key is exactly 'Radius'...
"""
import linecache
from flashy.utils import np
from flashy.IOutils import cl

def writeProf(filename, profD):
    """writes a 1D profile dict to a file."""
    points = len(profD['Radius'])
    header = " ".join(profD.keys())
    with open(filename, 'w') as f:
        f.write('# {}\n'.format(header))
        f.write('{}\n'.format(points))
        for i in range(points):
            line = []
            for k in profD.keys():
                line.append('{:15.8e}'.format(profD[k][i]))
            f.write("{}\n".format(" ".join(line)))
    print 'Wrote: {}'.format(filename)


def getProps(profD):
    """returns central dens, temp, rmin, rmax, volume and mass of a 1D profile."""
    pts = len(profD['Radius'])
    mass, vol = 0 , 0
    rads = np.insert(np.copy(profD['Radius']), 0, 0.0)
    for i in range(1, pts):
        r2, r1 = profD['Radius'][i], profD['Radius'][i-1]
        dr = r2-r1
        dvol = dr * ( 3.0*r1*r2 + dr*dr ) * 4.0*np.pi/3.0
        vol += dvol
        mass += dvol*profD['dens'][i-1]
    cd, ct = profD['dens'][0], profD['temp'][0]
    msuns = mass/(1.989e33)
    res = (r2 - profD['Radius'][0])/pts
    print "Points/resolution: {} / {:E}".format(pts, res)
    print "Central Dens: {:E}".format(cd)
    print "Central Temp: {:E}".format(ct)
    print "Radius: {:E}".format(r2)
    print "Volume: {:E}".format(vol)
    print "Mass (Msun): {:E}".format(msuns)
    return pts, res, cd, ct, r2, vol, msuns

# auxiliary functions
def getProfData(filename):
    """Returns arrays for radius, dens, temp and a matrix with abundances, 
    as an ordered dictionary.
    
    """
    data = np.genfromtxt(filename, skip_header=2, comments='#')
    header = linecache.getline(filename, 1).strip(' #\n').split()
    varsinFile = cl.OrderedDict(zip(header, range(len(header))))
    outDict = cl.OrderedDict()
    for k, v in varsinFile.items():
        var = data[:, v]
        outDict[k] = var
    return outDict
