from .utils import np
from .IOutils import os, getFileList, rename, _cdxfolder
from .IOutils import blockGenerator, getLines, getBlock, getRepeatedBlocks
import datetime
import tempfile
from .paramSetup import (parameterGroup, _FLASHdefaults,
                         _otpfolder, _statsfile, _logfile)
from flashy.datahaul.hdfdirect import turn2cartesian
_DEFpar = 'flash.par'  # default par name


class simulation(object):
    def __init__(self, folder):
        name = folder.rstrip('/')
        self.name = os.path.basename(name)
        self.cdx = os.path.join(os.path.dirname(name), _cdxfolder)
        self.chk = os.path.join(name, _otpfolder)
        # parameters (user set and defs)
        self.pargroup = parameterGroup(os.path.join(self.cdx, _FLASHdefaults))
        self.pargroup.setPars(os.path.join(name, _DEFpar))
        prof = self.pargroup.defaults.initialWDFile['value']
        self.profile = os.path.join(self.cdx, prof)
        # backwards compatibility: change names of otp
        if not getFileList(folder, glob=_logfile) == [_logfile]:
            rename(folder, '.log', _logfile)
            rename(folder, '.dat', _statsfile)
        self.pargroup.meta.update(buildMetaFromLog(os.path.join(folder,
                                                                _logfile)))
        self.netpath = os.path.join(self.cdx, 'Networks',
                                    self.pargroup.meta['network'])
        if self.pargroup.meta['geometry'] == 'cylindrical':
            self.standardizeGeometry()
            self.checkpoints = getFileList(self.chk,
                                           glob='cart_flash_hdf5_chk',
                                           fullpath=True)
            self.plotfiles = getFileList(self.chk,
                                         glob='cart_flash_hdf5_plt_',
                                         fullpath=True)
        else:
            self.checkpoints = getFileList(self.chk, glob='chk', fullpath=True)
            self.plotfiles = getFileList(self.chk, glob='plt', fullpath=True)
        try:
            self.steps, self.refs, self.chkp, self.header, self.timings = \
                readLogAndStats(os.path.join(name, _logfile),
                                os.path.join(name, _statsfile))
        except Exception as e:
            self.steps, self.header, self.timings = [], [], []
        self.time = self.getStepProp('t')
        addIRLdeltas(self.steps, outlier=5.0)
        deltas = self.getStepProp('irldelta')
        self.runtime = sum(deltas, datetime.timedelta(0))
        self.irlstep = np.mean(deltas[:-1])
        # qsub output parsing (.o files)
        glob = os.path.basename(name) + '.o'
        self.otpfiles = getFileList(folder, glob=glob, fullpath=True)
        # treat files as overwriting addenda to the timesteps
        if self.otpfiles:
            for f in self.otpfiles:
                self.addOtp(f)
        glob = 'shockDetect_'  # cj output parsing (specific .dat files)
        self.CJ = getFileList(folder, glob=glob, fullpath=True)
        yieldf = os.path.join(folder, 'spec')   # yield files
        try:
            glob = '.yield'
            self.yields = getFileList(yieldf, glob=glob, fullpath=True)
            glob = '.decayed'
            self.decayedYields = getFileList(yieldf, glob=glob, fullpath=True)
        except FileNotFoundError:
            self.yields = []
            self.decayedyields = []

    def getStepProp(self, which, range=[0, 0], src='steps'):
        """get a property of all steps in a run.

        Args:
            which(str): property to extract.
            range(int list): pick steps within range.
            src(str): source of data, either 'steps' or 'refs'

        Returns:
            (list): step ordered property list.

        """
        if sum(range) != 0:
            cut = slice(*range)
            return np.array([getattr(t, which)
                             for t in getattr(self, src)[cut]])
        else:
            return np.array([getattr(t, which)
                             for t in getattr(self, src)])

    def standardizeGeometry(self, verbose=False):
        """convert hdf5 files from cylindrical to cartesian."""
        if self.pargroup.params.geometry.strip('"') == 'cylindrical':
            turn2cartesian(self.chk, prefix='all',
                           nowitness=False, silent=True)
        elif verbose:
            print('Geometry already cartesian or spherical. '
                  'Skipping conversion.')

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
        otpls = getLines(filename, ') |  ')
        # FLASH magically changes de dts shown in the otp
        # noticeable by printing a new header line
        dtnames = otpls[0].split('|')[-1]
        dtnames = dtnames.split()
        dtdict = dict(zip(dtnames, [0.0]*len(dtnames)))
        stags = ['slowx', 'slowy', 'slowz']
        curdtn = dtnames
        nums = np.array([step.n for step in self.steps])
        for otpl in otpls[1:]:
            try:
                num, dtv, slowc = otpBreakdown(otpl)
            except ValueError:
                curdtn = otpl.split('|')[-1]
                curdtn = curdtn.split()
                continue
            valdict = dict(zip(curdtn, dtv))
            dtdict.update(valdict)
            loc = np.where((nums-num) == 0)[0]
            if loc.size > 0:
                loc = loc[0]
            else:
                # this should be the last step which has no
                # corresponding log step
                continue
            target = self.steps[loc]
            for t, v in zip(stags, slowc):
                setattr(target, t, v)
            for t, v in dtdict.items():
                setattr(target, t, v)
        # finally set first log step to zero (log-otp offset)
        for t in dtnames+stags:
            setattr(self.steps[0], t, 0.0)

    def quickLook(self, refsimtime=0.1, refstep=100, rsets=6, rounding=8):
        """print a group of general information about the run.

        Args:
            refsimtime(float): simtime to backtrace steps and walltime.
            refstep(int): step to backtrace simtime and walltime.
            rsets(int): assumed MPI Ranks/node based on machine.
            rounding(int): precision for max walltime.

        Returns:
            (str): string block of information about the run.

        """
        nodes, info = self.pargroup.probeSimulation(verbose=False)
        # nodes from probe may change if a forced PE count is used,
        # that meta is lost so get PEs from log file (3rd row).
        info[-1] = ''  # erase node recommendation from probeSim
        nodes = int(self.header.split('\n')[2].split()[-1])
        info.insert(len(info), 'Nodes used: {:>19}'.format(int(nodes/rsets)))
        info.insert(len(info), 'MPI Ranks: {:>20}'.format(nodes))
        info.append('Achieved timesteps: {:>11}'.format(len(self.steps)))
        tmax = self.getStepProp('t')[-1] + self.getStepProp('dt')[-1]
        info.append('Achieved simtime: {:>13}'.format(round(tmax, rounding)))
        info.append('Tstep sizes (s) [min/max(mean)]: '
                    '{:e}/{:e} ({:>10e})'.format(*self.getAvgTstep()))
        info.append('Total runtime (hh:mm:ss): {}'.format(str(self.runtime)))
        info.append('IRL timestep (hh:mm:ss): {}'.format(self.irlstep))
        limblks = np.max(self.getStepProp('totblocks', src='refs'))
        info.append('Max blocks used (per node): '
                    '{} ({:.2f})'.format(limblks, limblks/nodes))
        info.append('Checkpoints/Plotfiles: '
                    '{}/{}'.format(len(self.checkpoints), len(self.plotfiles)))
        # simulation figure of merit
        # info += self.getTfom(refsimtime)
        # info += self.getNfom(refstep)
        return '\n'.join(info)

    def getTfom(self, refsimt, tol=1e-3):
        """get a figure of merit for the simulation:
        for a given simtime, get irl step and steps achieved.
        simtime is fixed, so:
        ++ higher mean tstep
        ++ more steps achieved
        -- larger min max variation
        """
        N, simtime = self.time2step(refsimt)
        simdelta = simtime-refsimt
        if abs(simdelta) > tol:
            return ["\n Reference simulation time not reached."]
            # return 0, int(N), datetime.timedelta(0), datetime.timedelta(0)
        tstamps = [step.timestamp for step in self.steps[:int(N)]]
        start = tstamps[0]
        tstamps.insert(0, start)
        deltas = np.array([b-a for (a, b) in zip(tstamps[:-1], tstamps[1:])])
        totalwall = sum(deltas, datetime.timedelta(0))
        irlstep = np.mean(deltas)
        ts = self.getStepProp('dt')
        breadth = ts[:int(N)]
        factor = 1e5
        avg = np.mean(breadth)  # order of 1e-5
        tmin = np.min(breadth)
        tmax = np.max(breadth)
        delta = (tmin-tmax)
        components = [delta*factor, avg*factor, N, factor*(simtime-refsimt)]
        otp = []
        otp.append('\nStats for {:.4f}s Simtime:'.format(refsimt))
        otp.append('Reference steps achieved: {}'.format(int(N)))
        otp.append('IRL timestep (to {}): {} '.format(int(N), irlstep))
        otp.append('Walltime to reference: {}'.format(totalwall))
        otp.append('Figure of Merit: {}'.format(int(sum(components))))
        return otp

    def getNfom(self, refstep):
        try:
            simtime = self.steps[refstep].t
        except IndexError:
            return ["\n Reference simulation step not reached."]
        tstamps = [step.timestamp for step in self.steps[:refstep]]
        start = tstamps[0]
        tstamps.insert(0, start)
        deltas = np.array([b-a for (a, b) in zip(tstamps[:-1], tstamps[1:])])
        totalwall = sum(deltas, datetime.timedelta(0))
        irlstep = np.mean(deltas)
        ts = self.getStepProp('dt')
        breadth = ts[:refstep]
        factor = 1e5
        avg = np.mean(breadth)  # order of 1e-5
        tmin = np.min(breadth)
        tmax = np.max(breadth)
        delta = (tmin-tmax)
        components = [delta*factor, avg*factor, refstep]
        otp = []
        otp.append('\nStats for {} Steps:'.format(refstep))
        otp.append('Simtime reached: {:.4f}'.format(simtime))
        otp.append('IRL timestep (to {:.4f}): {} '.format(simtime, irlstep))
        otp.append('Walltime to reference: {}'.format(totalwall))
        otp.append('Figure of Merit: {}'.format(int(sum(components))))
        return otp

    def time2step(self, time):
        """returns step number correspoding to a given time."""
        for i, t in enumerate(self.time):
            if time > t:
                continue
            else:
                return self.steps[i].n, self.steps[i].t
        return self.steps[-1].n-1, self.steps[-1].t


def readLogAndStats(logfile, statsfile):
    """reads both log and stats file,
    correlating the data to each simulation step.
    """
    # read log data
    rsteps, skips, refs, chkps, header, timings = readLog(logfile)

    # read stats data and remove pre-restart values
    data = np.genfromtxt(statsfile, names=True)  # this skips # lines
    for att in data.dtype.names:
        # remove repeated steps detected in run.log from the stats data rows
        purgedAtt = np.delete(data[att], skips)
        for step, p in zip(rsteps, purgedAtt):
            setattr(step, att, p)
    return rsteps, refs, chkps, header, timings


def readLog(logfile):
    """read a flash log file, building step objects,
    grepping the first header, listing available timings
    and associating refinements to the steps.
    The method used here takes into account for
    repeated timesteps (which happpen after restarts)
    but not for repeated refinements which is likely for failing
    runs but not dramatic for the delta in blocks.

    Args:
        logfile(str): filepath to *.log file from flash run.

    Return:
        (sim.steps): every numbered step in the run (no repeated steps).
        (str): first flash header from the log file.
        (str list): end of run timings when available.
        (int list): rows of repeated steps(this is used to clean stats.dat).

    """
    # get first header from title to the first logged step
    header = getBlock(logfile, 'FLASH log file', ' [ ', skip=0)
    header = '\n'.join(header)
    # get all timings found in the log
    timings = getRepeatedBlocks(logfile,
                                'perf_summary: code performance summary for',
                                'LOGFILE_END: FLASH run complete.')
    timings = ['\n'.join(t) for t in timings]
    # log steps
    steps = [stepBreakdown(sline) for sline in getLines(logfile, 'step: ')]
    realsteps, removedrows = clearRestarts(steps)
    # refinement data
    refinementLines = getLines(logfile, '[GRID amr_refine_derefine]')
    refdata = [refBreakdown(rfb) for rfb in blockGenerator(refinementLines)]
    tstamps = [step.timestamp for step in steps]
    zero = datetime.timedelta(0)
    tags = ['minblocks', 'maxblocks', 'totblocks',
            'leaf_minblocks', 'leaf_maxblocks', 'leaf_totblocks']
    # find the simsteps where the refinement took place
    corr = []
    reftime = datetime.datetime.now()
    numreft = np.array([(reftime - rlvl[0]).total_seconds()
                        for rlvl in refdata])
    numstept = np.array([(reftime - t).total_seconds() for t in tstamps])
    for n in numreft:
        delts = numstept - n
        pos = np.where(delts > 0.0)[0]
        if pos.size > 0:
            corr.append(pos[-1])
        else:
            corr.append(0)
    # group ref-steps with the data from the corresponding simtime steps
    refs = []
    for i, c in enumerate(corr):
        ref = tstep()
        setattr(ref, 'n', steps[c].n)
        setattr(ref, 't', steps[c].t)
        setattr(ref, 'timestamp', steps[c].timestamp)
        for t, v in zip(tags, refdata[i][1:]):
            setattr(ref, t, v)
        refs.append(ref)
    # checkpoint/plotfile timestamps
    chkLs = getLines(logfile, "[IO_write")
    chkdata = [chkBreakdown(rfb) for rfb in blockGenerator(chkLs, step=4)]
    chkp = []
    for tst, (typ, num) in chkdata:
        chst = tstep()
        setattr(chst, 'timestamp', tst)
        setattr(chst, 'filetype', typ)
        setattr(chst, 'number', num)
        chkp.append(chst)
    return realsteps, removedrows, refs, chkp, header, timings


def addIRLdeltas(steps, outlier=3.0):
    """calculate walltime deltas between timesteps.
    Args:
        steps(sim.step list): simulation log steps.
        outlier(float): outlier delta for restarts (in minutes).

    Returns:
        (timedelta, timedelta): total walltime, mean walltime step

    """
    tstamps = [step.timestamp for step in steps]
    deltas = [b-a for (a, b) in zip(tstamps[:-1], tstamps[1:])]
    deltas.append(datetime.timedelta(0))
    for st, d in zip(steps, deltas):
        if d > datetime.timedelta(minutes=outlier):
            setattr(st, 'irldelta', datetime.timedelta(0))
        else:
            setattr(st, 'irldelta', d)


class tstep(object):
    """storing object for step info."""
    def __setattr__(self, att, val):
        super().__setattr__(att, val)

    def items(self):
        keys = self.__dict__
        values = [getattr(self, att) for att in keys]
        return zip(keys, values)


def buildMetaFromLog(log):
    """Looks up relevant run information from a
    log file to update the simulation object.
    """
    keys = ['network', 'geometry', 'cells', 'maxblocks', 'dimension']
    units = getBlock(log, 'FLASH Units used:', '=========')[1:]
    speciesblk = getBlock(log, 'Species Constituents', '=========')[1:]
    if any(['XNet' in u for u in units]):
        if len(speciesblk) < 15:
            netname = 'alpha'
        else:
            netname = 'SN{}'.format(len(speciesblk))
        network = 'Data_{}'.format(netname)
    else:
        network = 'Other'
    lines = getLines(log, 'geometry')
    for geom in ['cylindrical', 'spherical', 'cartesian']:
        if any([geom in line for line in lines]):
            break
    lines = getLines(log, 'zones')
    cells = []
    for l in lines[:3]:
        cells.append(int(l.split()[-1]))
    lines = getLines(log, 'Dimensionality')
    dimension = int(lines[0].split()[-1])
    lines = getLines(log, 'Max Number of Blocks/Proc')
    maxblocks = float(lines[0].split()[-1])
    return dict(zip(keys, [network, geom, cells, maxblocks, dimension]))


def readTiming(timingString, avgProc=False):
    """return timing data from a specified timing block.
    timing block has 2 subblocks, first one is for a single
    proc while the second is an avg from all procs.

    Args:
        timingString(str list): timing block to analize.
        avgProc(bool): return "second" timing subblock
        (avg over all procs).

    Returns:
        (list list): data per line with depths at pos 1.
        struct (time, steps, colnames, (unit, depth, [colvalues]))

    """
    # simulate a file
    tp = tempfile.NamedTemporaryFile(mode='w')
    tp.write(timingString)
    tp.flush()
    # chop lines
    time = getLines(tp.name, 'seconds in monitoring period')
    time = float(time[0].split(':')[-1])
    steps = getLines(tp.name, 'number of evolution steps')
    steps = int(steps[0].split(':')[-1])
    tblocks = getRepeatedBlocks(tp.name, 'accounting unit', '======')
    tblock = tblocks[-1] if avgProc else tblocks[0]
    # get data
    div = 35  # split for name and values
    colkeys = tblock[0].split()
    keys1 = colkeys[::2]
    keys2 = colkeys[1::2]
    colnames = [a+" "+b for a, b in zip(keys1, keys2)]
    data = []
    for line in tblock[2:]:
        unit, vals = line[:div], line[div:]
        ldep = len(unit) - len(unit.lstrip(' '))
        data.append((unit, ldep-1, [float(v) for v in vals.split()]))
    return time, steps, colnames, data


def clearRestarts(steps):
    """Remove repeated steps from restarting failed runs.

    Args:
        (sim.step list): list of simulation.step objects.

    Returns:
        (sim.step list): cleansed steps.
        (int list): removed step indices.

    """
    numbers = np.array([step.n for step in steps])  # get step numbers
    checkdelt = np.diff(numbers - len(numbers))  # find restarts
    restartpos = np.where(checkdelt != 1)[0]  # get restart positions
    repsteps = abs(checkdelt[restartpos])+1  # get wasted steps
    remsteps = []
    for pos, reps in zip(restartpos, repsteps):
        for i in range(int(reps)):
            remsteps.append(pos-i)
    realsteps = [st for i, st in enumerate(steps) if i not in remsteps]
    return realsteps, remsteps


def chkBreakdown(refBlock):
    """specific breakdown of a checkpoint/plotfile write log entry."""
    # get tstamp and type from first line
    tstampstr, _, rest = refBlock[0].partition(']')
    tstamp = datetime.datetime.strptime(tstampstr, '[ %m-%d-%Y %H:%M:%S.%f ')
    _, _, ftype = rest.partition('=')
    # get number from second line
    number = int(refBlock[1][-4:])
    return tstamp, (ftype, number)


def refBreakdown(refBlock):
    """detailed breakdown of the refinement block from a flash log.
    e.g.:
    '[ 06-17-2019  18:15:06.134 ]
        [GRID amr_refine_derefine]: initiating refinement'
    '[ 06-17-2019  18:15:06.144 ] [GRID amr_refine_derefine]:
        redist. phase.  tot blks requested: 10484'
    '[GRID amr_refine_derefine]
        min blks 39    max blks 46    tot blks 10484'
    '[GRID amr_refine_derefine] min leaf blks 30
        max leaf blks 38    tot leaf blks 7895'
    '[ 06-17-2019  18:15:06.185 ]
        [GRID amr_refine_derefine]: refinement complete'

    Arg:
        (str list): 5 line refinement log stamps.

    Returns:
        (datetime.datetime, 6 floats): timestamp,
        min/max/tot blocks base and leaf.

    """
    blockstencil = [2, 5, 8]
    leafbstencil = [3, 7, 11]
    # get tstamp from first line
    tstampstr, _, _ = refBlock[0].partition(']')
    tstamp = datetime.datetime.strptime(tstampstr, '[ %m-%d-%Y %H:%M:%S.%f ')
    # block information from lines 3-4
    _, _, blockstr = refBlock[2].partition(']')
    min, max, tot = [int(x) for i, x in enumerate(blockstr.split())
                     if i in blockstencil]
    _, _, blockstr = refBlock[3].partition(']')
    lmin, lmax, ltot = [int(x) for i, x in enumerate(blockstr.split())
                        if i in leafbstencil]
    return tstamp, min, max, tot, lmin, lmax, ltot


def stepBreakdown(sline):
    """detailed breakdown a step from a flash log.
    e.g.:
    '[ 06-17-2019  18:15:07.507 ] step: n=46142 t=1.314454E+00 dt=1.000000E-10'

    Args:
        (str): log line

    Returns:
        (simulation.step): step object with relevant data.

    """
    tags = ['timestamp', 'n', 't', 'dt']
    tstampstr, _, props = sline.partition(']')
    values = [datetime.datetime.strptime(tstampstr.strip('['),
                                         ' %m-%d-%Y %H:%M:%S.%f ')]
    pstr = props.replace('=', ' ')
    values += [float(x) for x in pstr.split() if x[0].isdigit()]
    simstep = tstep()
    for t, v in zip(tags, values):
        setattr(simstep, t, v)
    return simstep


def otpBreakdown(otpline):
    """detailed breakdown of a bash output file from flash.
    e.g.:
    '  14 5.9196E-07 1.2839E-07  '
    '( 5.243E+05,  4.168E+08,   0.00    ) |  2.966E-04 1.158E-05'
    Note that there's an offset in the ouput with respect to the log
    since the first timestep is not written to std out.

    Args:
        (str): step from the bash output.

    Returns:
        (int): log step number (otp number + 1).
        (float list): timestep values found.
        (float list): slowest coordinates (x, y, z).

    """
    stepinfo, _, extradts = otpline.partition('|')
    extradts = extradts.split()
    stepprops, _, pos = stepinfo.strip(' )').partition('(')
    # t and dt already in log so discard
    otpn, _, _ = stepprops.split()
    x, y, z = pos.split(',')
    otpn = int(otpn)+1  # log vs otp offset
    slowc = [float(i.replace(',', '')) for i in [x, y, z]]
    # take the first dt and float it
    hydrodt = float(extradts[0])
    try:
        otherdts = [float(exdt) for exdt in extradts[1:]]
    except ValueError:
        otherdts = [-1.0]*len(extradts[1:])
    dtvals = [hydrodt]+otherdts
    return otpn, dtvals, slowc
