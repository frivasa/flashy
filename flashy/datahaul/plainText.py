"""module for handling plain text 1D profiles sorted in columns and with structure:
# col names
length of data(rows int)
<data block>
# comments or anything
"""
import linecache
from flashy.utils import np, byMass
from flashy.IOutils import cl
from flashy.utils import msol, byMass


class dataMatrix(object):
    """object for interacting with columned data in a plain text file.
    format for file is strict at the start since it uses np.genfromtxt to handle 
    the data block. after row 3 any comment is allowed (# marked).
    Essential columns: radius and density.
    Anything with a number becomes a species.
    
    # col names
    length of data(rows int)
    <data block>
    # comments or anything
    """
    def __init__(self, filename, comment=''):
        if isinstance(filename, list):
            self.filekeys, self.data = filename[0], filename[1]
        else:
            self.filekeys, self.data = chopFile(filename)
        self.species = []
        self.bulkprops = []
        self.meta = {}
        skip = []
        n, num = multiIndex(['radius', 'r', 'rad'], self.filekeys)
        self.radius = self.data[:, num]
        self.bulkprops.append('radius')
        skip.append(num)
        n, num = multiIndex(['rho', 'dens', 'density'], self.filekeys)
        self.density = self.data[:, num]
        skip.append(num)
        self.bulkprops.append('density')
        
        for i, k in enumerate(self.filekeys):
            if i not in skip:
                setattr(self, k, self.data[:, i])
                if any(i.isdigit() for i in k):
                    self.species.append(k)
                else:
                    self.bulkprops.append(k)
            else:
                continue
        self.masses = byMass(self.radius, self.density)
        self.bulkprops.append('masses')
        self.meta['mass'] = self.masses[-1]
        self.meta['central dens'] = self.density[0]
        self.meta['radius'] = self.radius[-1]
        self.meta['resolution'] = (self.radius[-1] - self.radius[0])/len(self.radius)
        self.meta['points'] = len(self.radius)
        self.meta['comment'] = comment

    def printMeta(self):
        for k, v in self.meta.items():
            print '{}: {:E}'.format(k.capitalize(), v)

    def writeProf(self, output, subset=[]):
        """Write profile to file
        
        Args:
            ouput(str): otp filename.
            subset(list): write a subset of keys to file.
        
        """
        if subset:
            keys = subset
        else:
            keys = self.bulkprops + self.species
        print "Writing: {}".format(" ".join(keys))
        missing = set()
        with open(output, 'w') as f:
            header = " ".join(keys)
            f.write('# {}\n'.format(header))
            f.write('{}\n'.format(self.meta['points']))
            for i in range(self.meta['points']):
                line = []
                for k in keys:
                    try:
                        line.append('{:15.8e}'.format(getattr(self, k)[i]))
                    except AttributeError:
                        missing.add(k)
                        line.append('{:15.8e}'.format(0.0))
                f.write("{}\n".format(" ".join(line)))
            if self.meta['comment']:
                f.write("# {}\n".format(self.meta['comment']))
            f.write("# Total Mass {} Msun".format(self.meta['mass']))
        if missing:
            print 'Missing keys: {} were set to zero.'.format(' '.join(missing))
        print 'Wrote: {}'.format(output)


def multiIndex(names, keys):
    """returns the index of a key with multiple probable names in a 
    list. (assumes such key only happens once)"""
    exists = False
    for n in names:
        try:
            num = keys.index(n)
            exists = True
            break
        except ValueError:
            continue
    if exists:
        return n, num
    else:
        return n, -1


def chopFile(filename):
    """Returns header names and a data matrix from a file.
    
    Args:
        filename(str): file path.
    
    Returns:
        (list): lowercase header names.
        (np.array): data matrix of shape (coords, properties)
    
    """
    data = np.genfromtxt(filename, skip_header=2, comments='#')
    header = linecache.getline(filename, 1).strip(' #\n').split()
    return [h.lower() for h in header], data
