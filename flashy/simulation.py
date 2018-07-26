from .utils import np
from .IOutils import os, getFileList, subp, _cdxfolder
import datetime
from .paramSetup import parameterGroup, _FLASHdefaults, _otpfolder
import flashy.plot.simplot as simp
from flashy.datahaul.hdfdirect import turn2cartesian
_pancake = 'PANCAKE'  # break word

class simulation(object):
    def __init__(self, folder):
        name = folder.rstrip('/')
        # folders and paths
        self.name = name
        self.root = os.path.dirname(name)
        self.cdx = os.path.join(self.root, _cdxfolder)
        self.chk = os.path.join(self.name, _otpfolder)
        self.checkpoints = getFileList(self.chk, glob='cart_flash_hdf5_chk', fullpath=True)
        if not self.checkpoints:
            self.checkpoints = getFileList(self.chk, glob='chk', fullpath=True)
        self.plotfiles = getFileList(self.chk, glob='plt', fullpath=True)
        # parameters
        self.pargroup = parameterGroup(os.path.join(self.cdx, _FLASHdefaults))
        self.pargroup.setPars(os.path.join(self.name, 'flash.par'))
        # check for Xnet
        if os.path.exists(os.path.join(self.cdx, 'Networks')):
            netname = 'Data_{}'.format(self.pargroup.meta['network'].split('_')[-1])
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
        else:
            self.runtime, self.rundelts = None, None
        # stats (stats.dat output)
        data = readStats(os.path.join(name, 'stats.dat'))
        for att in data.dtype.names:
            setattr(self, att, data[att][:len(self.steps)])
        # qsub output parsing (.o files)
        glob = os.path.basename(name) + '.o'
        self.otpfiles = getFileList(folder, glob=glob, fullpath=True)
        if self.otpfiles:
            for f in self.otpfiles:
                self.addOtp(f)
        # cj output parsing (specific .dat files)
        glob = 'shockDetect_'
        self.CJ = getFileList(folder, glob=glob, fullpath=True)
        # yield files
        yieldf = os.path.join(folder, 'spec')
        try:
            glob = '.yield'
            self.yields = getFileList(yieldf, glob=glob, fullpath=True)
            glob = '.decayed'
            self.decayedYields = getFileList(yieldf, glob=glob, fullpath=True)
        except FileNotFoundError:
            self.yields = []
            self.decayedyields = []

    def getStepProp(self, which):
        """get a property of all steps in a run.

        Args:
            which(str): property to extract.

        Returns:
            (list): step ordered property list.

        """
        return [getattr(t, which) for t in self.steps]

    def standardizeGeometry(self):
        """convert hdf5 files from cylindrical to cartesian."""
        if self.pargroup.params.geometry.strip('"')=='cylindrical':
            turn2cartesian(self.chk, prefix='all', nowitness=False)
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
        ts = self.getStepProp('dt')
        low = np.min(ts)
        high = np.max(ts)
        avg = np.mean(ts)
        return low, high, avg

    def addOtp(self, filename):
        """read in a qsub output file. getting slowest coordinates for each timestep."""
        ns, slowps, dthy, dtbu = readOtp(filename)
        if not ns:
            print('Qsub stopped in {} or No steps achieved.'.format(filename))
            #os.remove(filename)
        for (n, p, dth, dtb) in zip(ns, slowps, dthy, dtbu):
            try:
                setattr(self.steps[n], 'slowx', p[0])
                setattr(self.steps[n], 'slowy', p[1])
                setattr(self.steps[n], 'slowz', p[2])
                setattr(self.steps[n], 'dthydro', dth)
                setattr(self.steps[n], 'dtburn', dtb)
            except IndexError:
                continue

    def quickLook(self):
        """print a group of general information about the run."""
        nodes, info = self.pargroup.probeSimulation(verbose=False)
        # nodes from probe may change if a forced PE count is used,
        # that meta is lost so get PEs from log file (3rd row).
        nodes = int(self.header.split('\n')[2].split()[-1])
        info[-1] = 'Nodes used: {}'.format(nodes)
        info.append('Achieved timesteps: {}'.format(len(self.time)))
        info.append('Tstep sizes (s): min {:e} max {:e} mean {:e}'.format(*self.getAvgTstep()))
        if self.timings:
            info.append('Total runtime (hh:mm:ss): {}'.format(str(self.runtime)))
            info.append('IRL timestep (s): {:e}'.format(self.runtime.total_seconds()/len(self.time)))
        limblks = np.max(self.getStepProp('blocks'))
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

    def mergeTimings(self):
        """Averages timing information."""
        props = ['secs', 'calls', 'percent', 'depth']
        if not self.timings:
            return {}
        #elif len(self.timings)==1:
        else:  # just use the first one.
            tsteps, ttime, md = readTiming(self.timings[0])
            return tsteps, ttime, ttime/tsteps, md
        # some keys change in timings. skipping merging for now
        tims = []
        for timing in self.timings:
            tims.append(readTiming(timing))
        tsteps = sum([a[0] for a in tims])
        ttime = sum([a[1] for a in tims])
        fmd = tims[0][2]
        avg = len(tims)
        print(fmd.keys())
        for md in [a[2] for a in tims[1:]]:
            print (md.keys())
            for k in md.keys():
                print(k)
                fmd[k]['secs'] += md[k]['secs']
                fmd[k]['calls'] += md[k]['calls']
                fmd[k]['percent'] += md[k]['percent']
                for j in md[k].keys():
                    if j not in props:
                        fmd[k][j]['secs'] += md[k][j]['secs']
                        fmd[k][j]['calls'] += md[k][j]['calls']
                        fmd[k][j]['percent'] += md[k][j]['percent']
        for k in fmd.keys():
            fmd[k]['secs'] /= avg
            fmd[k]['calls'] /= avg
            fmd[k]['percent'] /= avg
            for j in fmd[k].keys():
                if j not in props:
                    fmd[k][j]['secs'] /= avg
                    fmd[k][j]['calls'] /= avg
                    fmd[k][j]['percent'] /= avg
        return tsteps, ttime, ttime/tsteps, fmd


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
    try:
        header, timings = splitHeader(hlines)
        timings = ['\n'.join(t) for t in timings]
    except IndexError:
        header = hlines
        timings = None
    return steps, '\n'.join(header), timings


def splitHeader(hlines):
    """separates the header from timings in a list of .log lines."""
    perfs, heads, runs = [], [], []
    for i, l in enumerate(hlines):
        if 'FLASH log' in l:
            heads.append(i)
        # average timers may fail sometimes so you get only percentages.
        elif 'perf_summary: code performance summary for process' in l:
            perfs.append(i)
    header = hlines[heads[0]:perfs[0]]
    timings = []
    heads.append(len(hlines))
    for i, j in enumerate(perfs):
        timings.append(hlines[j:heads[i+1]])
    return header, timings


def readOtp(filename):
    """reads an output file (.oNUMBER) extracting the slowest coordinate for
    each step."""
    with open(filename, 'r') as f:
        ns, slowp, dthydro, dtburn = [], [], [], []
        for i, l in enumerate(f):
            if _pancake in l:
                return [], [], [], []
            if '|' in l:
                if 'x' in l:
                    continue
                else:
                    n = l.split()[0]
                    a, b = l.index('('), l.index(')')
                    ns.append(int(n)-1)
                    slowp.append([float(x) for x in l[a+1:b].split(',')])
                    dth, dtb = l.split()[-2:]
                    dthydro.append(float(dth))
                    dtburn.append(float(dtb))
        return ns, slowp, dthydro, dtburn


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


def readTiming(timingString):
    lim = []
    depth, indent = 0, 0
    for i, l in enumerate(timingString.split('\n')):
    #     depth = len(l) - len(l.lstrip(' '))
        if 'seconds in monitoring period' in l:
            tdelt = float(l.split()[-1])
        elif 'number of evolution steps' in l:
            steps = int(l.split()[-1])
        elif '---' in l:
            start = i
        elif '====' in l:
            stop = i
            break
        else:
            indent = len(l) - len(l.lstrip(' '))
            if indent>depth and indent<20:  # first lines have large indents
                depth = indent
    md = {}
    for i, l in enumerate(timingString.split('\n')[start+1:stop]):
        specialsplt = [a for a in l.split('  ') if a!='']
        delt = len(l) - len(l.lstrip(' '))
        mp, secs, calls, avgt, perc = specialsplt
        mp = mp.strip().replace(' ', '_')
        if delt==1:
            md[mp] = {}
            md[mp]['secs'] = float(secs)
            md[mp]['calls'] = int(calls)
            md[mp]['percent'] = float(perc)
            stem = mp
        else:
            md[stem][mp]={}
            md[stem][mp]['secs'] = float(secs)
            md[stem][mp]['calls'] = int(calls)
            md[stem][mp]['percent'] = float(perc)
            md[stem][mp]['depth'] = delt
    return steps, tdelt, md

