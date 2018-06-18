from .utils import np
from .IOutils import os, getFileList, turn2cartesian, subp, _cdxfolder
import datetime
from .paramSetup import parameterGroup, _FLASHdefaults, _otpfolder
import flashy.plot.simplot as simp
_pancake = 'PANCAKE'  # break word

class simulation(object):
    def __init__(self, folder):
        name = folder.rstrip('/')
        # folders and paths
        self.name = name
        self.root = os.path.dirname(name)
        self.cdx = os.path.join(self.root, _cdxfolder)
        self.chk = os.path.join(self.name, _otpfolder)
        self.checkpoints = getFileList(self.chk, glob='chk', fullpath=True)
        self.plotfiles = getFileList(self.chk, glob='plt', fullpath=True)
        # parameters
        self.pargroup = parameterGroup(os.path.join(self.cdx, _FLASHdefaults))
        self.pargroup.setPars(os.path.join(self.name, 'flash.par'))
        # check for Xnet
        if os.path.exists(os.path.join(self.cdx, 'Networks')):
            netname = 'Data_{}'.format(self.pargroup.meta['network'])
            self.netpath = os.path.join(self.cdx, 'Networks', netname)
        else:
            self.netpath = ''
        # initial profile
        self.profile = os.path.join(self.cdx, self.pargroup.defaults.initialWDFile['value'])
        # log file (timings, headers, and timesteps)
        self.steps, self.header, self.timings = readLog(os.path.join(name, 'run.log'))
        if self.timings:
            self.rundelts = [timingParser(t.split('\n')) for t in self.timings]
            self.runtime = sum(self.rundelts, datetime.timedelta())
        # stats (.dat output)
        data = readStats(os.path.join(name, 'stats.dat'))
        for att in data.dtype.names:
            setattr(self, att, data[att][:len(self.steps)])
        # qsub output parsing (.o files)
        glob = os.path.basename(name)+'.o'
        self.otpfiles = getFileList(folder, glob=glob, fullpath=True)
        if self.otpfiles:
            for f in self.otpfiles:
                self.addOtp(f)

    def getBlocks(self):
        """returns the blocks used for each timestep."""
        return [t.blocks for t in self.steps]
    
    def getTsteps(self):
        """get all timesteps in the run."""
        return [t.dt for t in self.steps]

    def getTimes(self):
        """get all timesteps in the run."""
        return [t.t for t in self.steps]
    
    def standardizeGeometry(self):
        """convert hdf5 files from cylindrical to cartesian."""
        if self.pargroup.params.geometry.strip('"')=='cylindrical':
            turn2cartesian(self.chk, prefix='all', nowitness=True)
            self.checkpoints = getFileList(self.chk, glob='chk', fullpath=True)
            self.plotfiles = getFileList(self.chk, glob='plt')
        else:
            print('Geometry already cartesian or spherical. Skipping conversion.')
    
    def getSlowCoord(self, direction='x'):
        """returns slowest position for each timestep in a given axis."""
        times, dirs = [], []
        for t in self.steps:
            try:
                sdir = getattr(t, 'slow{}'.format(direction))
                times.append(t.t)
                dirs.append(sdir)
            except AttributeError:
                continue
        return dirs
    
    def getAvgTstep(self):
        ts = self.getTsteps()
        low = np.min(ts)
        high = np.max(ts)
        avg = np.mean(ts)
        return low, high, avg
    
    def addOtp(self, filename):
        """read in a qsub output file. getting slowest coordinates for each timestep."""
        ns, slowps = readOtp(filename)
        if not ns:
            print('Qsub stopped in {}. Deleting.'.format(filename))
            os.remove(filename)
        for (n, p) in zip(ns,slowps):
            try:
                setattr(self.steps[n], 'slowx', p[0])
                setattr(self.steps[n], 'slowy', p[1])
                setattr(self.steps[n], 'slowz', p[2])
            except IndexError:
                continue

    def quickLook(self):
        """print a group of general information about the run."""
        nodes, info = self.pargroup.probeSimulation(verbose=False)
        info[-1] = 'Nodes used: {}'.format(nodes)
        info.append('Total runtime (hh:mm:ss): {}'.format(str(self.runtime)))
        info.append('Timesteps (s): min {:e} max {:e} mean {:e}'.format(*self.getAvgTstep()))
        info.append('IRL timestep (s): {:e}'.format(self.runtime.seconds/len(self.time)))
        limblks = np.max(self.getBlocks())
        info.append('Max blocks used (per node): {} ({:.2f})'.format(limblks, limblks/nodes))
        info.append('Checkpoints: {}'.format(len(self.checkpoints)))
        info.append('Plotfiles: {}'.format(len(self.plotfiles)))
        return '\n'.join(info)

    def time2step(self, time):
        """returns step number correspoding to a given time."""
        for i, t in enumerate(self.time):
            if time>t:
                continue
            else:
                return self.steps[i].n, self.steps[i].t
        return self.steps[-1].n-1, self.steps[-1].t


def readStats(filename):
    """reads in data from a stats.dat file."""
    data = np.genfromtxt(filename, names=True)
    return data

    
def readLog(filename):
    """reads a FLASH log file filtering everything except step lines and header/timers"""
    steps = []
    hlines = []
    with open(filename, 'r') as f:
        for i, l in enumerate(f):
            if l.startswith(' ['):
                if 'GRID' in l and 'tot blks requested' in l:
                    blks = int(l.split(':')[-1])
                elif 'step' in l and 'GRID' in l:
                    continue
                elif 'checkpoint' in l:
                    continue
                elif 'step' in l:
                    tst = tstep()
                    setattr(tst, 'blocks', blks)
                    tstamp, params = chopLogline(l.strip())
                    setattr(tst, 'tstamp', tstamp)
                    for (k, v) in params:
                        setattr(tst, k, v)
                    steps.append(tst)
            else:
                hlines.append(l.strip('\n'))
    header, timings = splitHeader(hlines)
    timings = ['\n'.join(t) for t in timings]
    return steps, '\n'.join(header), timings


def splitHeader(hlines):
    """separates the header from timings in a list of .log lines."""
    perfs, heads, runs = [], [], []
    for i, l in enumerate(hlines):
        if 'FLASH log' in l:
            heads.append(i)
        elif 'perf_summary' in l:
            perfs.append(i)
    header = hlines[heads[0]:perfs[0]]
    timings = []
    heads.append(len(hlines))
    for i, j in enumerate(perfs[::2]):
        timings.append(hlines[j:heads[i+1]])
    return header, timings


def readOtp(filename):
    """reads an output file (.oNUMBER) extracting the slowest coordinate for 
    each step."""
    with open(filename, 'r') as f:
        ns, slowp = [], []
        for i, l in enumerate(f):
            if _pancake in l:
                return [], []
            if '|' in l:
                if 'x' in l:
                    continue
                else:
                    n = l.split()[0]
                    a, b = l.index('('), l.index(')')
                    ns.append(int(n))
                    slowp.append([float(x) for x in l[a+1:b].split(',')])
        return ns, slowp


def timingParser(tlines):
    """for now this returns the timespan of a timing. 
    expand to make stats from run."""
    start, stop = datetime.datetime(1, 1, 1), datetime.datetime(1, 1, 1)
    for t in tlines:
        if 'beginning : ' in t:
            _, _, tstamp = t.partition(':')
            start = datetime.datetime.strptime(tstamp.strip(), '%m-%d-%Y  %H:%M:%S')
        if 'ending :' in t:
            _, _, tstamp = t.partition(':')
            stop = datetime.datetime.strptime(tstamp.strip(), '%m-%d-%Y  %H:%M:%S')
    delt = stop-start
    return delt


class tstep(object):
    """storing object for step info."""
    def __setattr__(self, att, val):
        super().__setattr__(att, val)

    def items(self):
        keys = self.__dict__
        values = [getattr(self, att) for att in keys]
        return zip(keys, values)
    

def unstick(phrase):
    """future: use str.partition(separator) = pre, sep, post"""
    where = phrase.find('=')
    if where==-1:
        return -1, -1  # pass on the 'error'
    return phrase[:where], phrase[where+1:]


def chopLogline(string):
    """splits a step line from a .log file."""
    br1, br2 = string.find('['), string.find(']')  # assuming index if leftmost ocurrence
    if br1==-1:
        return -1, -1
    stamp = datetime.datetime.strptime(string[br1+1:br2].strip(), '%m-%d-%Y  %H:%M:%S.%f')
    remnant = string[br2:]
    params = []
    for p in remnant.split():
        k, v = unstick(p)
        if k==-1:
            continue # not a value
        else:
#             try:
            params.append((k.strip(), float(v)))
#             except ValueError:
#                 print('bad log line: {}'.format(v))
#                 params.append((k.strip(), 1.0e-20))
    return stamp, params
        