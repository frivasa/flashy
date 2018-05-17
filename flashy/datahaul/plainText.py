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
import os


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
    rnames = ['radius', 'r', 'rad']
    dnames = ['rho', 'dens', 'density']
    
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

    def __getattr__(self, name):
        # run name through the 'dictionary' of names
        if name in self.dnames:
            return self.__getattribute__('density')
        elif name in self.rnames:
            return self.__getattribute__('radius')
        else:
            # default call
            return self.__getattribute__(name)
    
    def printMeta(self):
        for k, v in self.meta.items():
            print('{}: {:E}'.format(k.capitalize(), v))

    def writeProf(self, output, subset=[], autotag=False, Xthresh=1e-15):
        """Write profile to file
        bug: ndarray ninja breaks formatting (np.array([float]))
        
        Args:
            ouput(str): otp filename.
            subset(list): write a subset of keys to file.
            autotag(bool): setup name based on the profile (WDs)
        
        """
        basedir = os.path.dirname(output)
        if autotag:
            rhoc = '{:2.4e}'.format(self.density[0]).replace('+','').replace('.','p')
            c12 = '{:2.2f}'.format(self.c12[0]).replace('+','').replace('.','p')
            o16 = '{:2.2f}'.format(self.o16[0]).replace('+','').replace('.','p')
            filename = 'wd_fermi_helm_{}_{}_{}.dat'.format(rhoc, c12, o16)
        else:
            filename = os.path.basename(output)
        otp = os.path.join(basedir, filename)
        if subset:
            keys = dataMatrix.checkbulkprops(subset)
        else:
            keys = self.bulkprops + self.species
        print ("Writing: {}".format(" ".join(keys)))
        missing = set()
        with open(otp, 'w') as f:
            header = " ".join(keys)
            f.write('# {}\n'.format(header))
            f.write('{}\n'.format(self.meta['points']))
            for i in range(self.meta['points']):
                line = []
                for k in keys:
                    try:
                        # clip stupid numbers such as 1e-99
                        value = getattr(self, k)[i]
                        if value < Xthresh or np.isnan(value):
                            value = 0.0
                        line.append('{:15.8e}'.format(value))
                    except AttributeError:
                        missing.add(k)
                        line.append('{:15.8e}'.format(0.0))
                f.write("{}\n".format(" ".join(line)))
            if self.meta['comment']:
                f.write("# {}\n".format(self.meta['comment']))
            f.write("# Total Mass {} Msun".format(self.meta['mass']))
        if missing:
            print ('Missing keys: {} were set to zero.'.format(' '.join(missing)))
        print ('Wrote: {}'.format(otp))
    
    def checkbulkprops(klist):
        nklist = []
        for k in klist:
            # hard bulk props are always named density and radius
            if k in dataMatrix.rnames:
                # key requested is a valid radius name
                nklist.append('radius')
            elif k in dataMatrix.dnames:
                # key requested is a valid density name
                nklist.append('density')
            else:
                # key is something else, careful with temp/temperature
                    nklist.append(k)
        return nklist
                


def joinProfiles(pori, pnew, skipP=True):
    """splices two dataMatrix objects. 
    probable bug: output profile may have a different resolution.
    
    """
    cut = pnew.radius[0]
    ncut = np.where(pori.radius<cut)[0][-1]+1
    if skipP:
        keys = ['radius', 'dens', 'temp']
        nra = np.hstack((pori.radius[:ncut], pnew.radius))
        nde = np.hstack((pori.density[:ncut], pnew.density))
        nte = np.hstack((pori.temp[:ncut], pnew.temp))
        dblock = np.column_stack([nra, nde, nte])
    else:
        keys = ['radius', 'dens', 'pres', 'temp']
        nra = np.hstack((pori.radius[:ncut], pnew.radius))
        nde = np.hstack((pori.density[:ncut], pnew.density))
        npr = np.hstack((pori.pres[:ncut], pnew.pres))
        nte = np.hstack((pori.temp[:ncut], pnew.temp))
        dblock = np.column_stack([nra, nde, npr, nte])

    species = pori.species + pnew.species
    species = set(species)
    keys = keys + list(species)
    flen = len(nra)
    for i, x in enumerate(species):
        try:
            osp = getattr(pori, x)[:ncut]
        except AttributeError:
            osp = ncut*[0.0]
            offset = 0
        try:
            nsp= getattr(pnew, x)
        except AttributeError:
            nsp = (flen-ncut)*[0.0]
        xmass = np.hstack((osp, nsp))
        dblock = np.column_stack((dblock, xmass))
    return dataMatrix([keys, dblock])

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
