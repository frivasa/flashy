from .utils import np, pd
from .IOutils import os, log, getFileList, rename, _cdxfolder, _metaname
from .IOutils import blockGenerator, getLines, getBlock, getRepeatedBlocks
import datetime
import tempfile
from .paramSetup import (parameterGroup, _FLASHdefaults,
                         _otpfolder, _statsfile, _logfile)
from flashy.datahaul.hdfdirect import turn2cartesian
_DEFpar = 'flash.par'  # default par name
log.name = __name__


class simulation(object):
    def __init__(self, folder):
        name = folder.rstrip('/')
        self.name = os.path.basename(name)
        self.cdx = os.path.join(os.path.dirname(name), _cdxfolder)
        self.chk = os.path.join(name, _otpfolder)

        # parameters (user set and defs)
        self.pargroup = parameterGroup(os.path.join(self.cdx, _FLASHdefaults))
        self.pargroup.setPars(os.path.join(name, _DEFpar))
        try:
            prof = self.pargroup.defaults.initialWDFile['value']
        except (TypeError, AttributeError) as e:
            prof = 'no profile, different kind of simulation.'
            
        self.profile = os.path.join(self.cdx, prof)

        # backwards compatibility: change names of otp
        if not getFileList(folder, glob=_logfile) == [_logfile]:
            rename(folder, '.log', _logfile)
            rename(folder, '.dat', _statsfile)
        try:
            inp = os.path.join(folder,_logfile)
            self.logfile = inp
            self.pargroup.meta.update(buildMetaFromLog(inp))
        except FileNotFoundError:
            self.logfile = ''
            pass

        # network and checkpoints
        self.netpath = os.path.join(self.cdx, 'Networks',
                                    self.pargroup.meta['network'])
        if self.pargroup.meta['geometry'] == 'cylindrical':
            self.standardize_geometry()
            self.checkpoints = getFileList(self.chk,
                                           glob='cart_flash_hdf5_chk',
                                           fullpath=True)
            self.plotfiles = getFileList(self.chk,
                                         glob='cart_flash_hdf5_plt_',
                                         fullpath=True)
            self.partfiles = getFileList(self.chk,
                                         glob='cart_flash_hdf5_part_',
                                         fullpath=True)
        else:
            self.checkpoints = getFileList(self.chk, glob='chk', fullpath=True)
            self.plotfiles = getFileList(self.chk, glob='plt', fullpath=True)
            self.partfiles = getFileList(self.chk, glob='part', fullpath=True)
        # add ouput files if any
        glob = os.path.basename(name) + '.o'
        self.otpfiles = getFileList(folder, glob=glob, fullpath=True)
        glob = os.path.basename(name) + '.e'
        self.errfiles = getFileList(folder, glob=glob, fullpath=True)

        # build step data from FLASH output (log and stats)
        try:
            log.debug('Reading metadata.')
            self.steps, self.chkp, self.header, self.timings = \
                read_log_and_stats(os.path.join(name, _logfile),
                                os.path.join(name, _statsfile))
            log.debug("Setting up time deltas")
            self.time = self.get_step_prop('t')
            deltas = getIRLdeltas(self.steps)
            self.steps = self.steps.assign(irldelta=pd.Series(deltas))
            self.runtime = sum(deltas, datetime.timedelta(0))
            self.irlstep = np.mean(deltas[:-1])
            log.debug("Adding otp files to timesteps")
            # Add qsub output and error (.o, .e files)
            # treat files as overwriting addenda to the timesteps
            if self.otpfiles:
                for i, f in enumerate(self.otpfiles):
                    try:  # skip any non-typical otp
                        if not i:
                            self.addOtp(f, init=True)
                        else:
                            self.addOtp(f)
                        log.debug("Added {}".format(os.path.basename(f)))
                    except IndexError:
                        log.debug("Skipped {} (IndexError)".format(os.path.basename(f)))
                        continue
                    except KeyError:
                        log.debug("Skipped {} (KeyError)".format(os.path.basename(f)))
                        continue
        except FileNotFoundError:
            log.critical('No steps set. (log file not found)', exc_info=False)
            self.steps, self.header, self.timings = [], [], []
        except Exception as e:
            log.debug('Skipped reading metadata.')
            log.critical('No steps set.', exc_info=True)
            self.steps, self.header, self.timings = [], [], []

        # CJ file (if any)
        glob = 'shockDetect_'  # cj output parsing (specific .dat files)
        self.CJ = getFileList(folder, glob=glob, fullpath=True)

        # yields (if any)
        yieldf = os.path.join(folder, 'spec')
        try:
            glob = '.yield'
            self.yields = getFileList(yieldf, glob=glob, fullpath=True)
            glob = '.decayed'
            self.decayedYields = getFileList(yieldf, glob=glob, fullpath=True)
        except FileNotFoundError:
            self.yields = []
            self.decayedyields = []

    def get_step_prop(self, which, range=[0, 0], src='steps'):
        """get a property of all steps in a run.

        Args:
            which(str): property to extract.
            range(int list): pick steps within range.
            src(str): source of data: 'steps', 'refs' or 'chkp'

        Returns:
            (list): step ordered property list.

        """
        if sum(range) != 0:
            cut = slice(*range)
            return np.array(self.steps[which])[cut]
        else:
            return np.array(self.steps[which])

    def print_step_props(self):
        """print all keys available to the steps."""
        selfnamelist = ['steps', 'refs', 'chkp']
        for sname in selfnamelist:
            try:
                _all = getattr(self, sname)
                print(sname, len(_all))
                print(_all.keys())
            except AttributeError:
                print('no {} in simulation.object'.format(sname))

    def standardize_geometry(self, verbose=False):
        """convert hdf5 files from cylindrical to cartesian."""
        if self.pargroup.params.geometry.strip('"') == 'cylindrical':
            turn2cartesian(self.chk, prefix='all',
                           nowitness=False, silent=True)
        elif verbose:
            print('Geometry already cartesian or spherical. '
                  'Skipping conversion.')

    def get_slow_coord(self, direction='x'):
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

    def get_avg_tstep(self):
        ts = self.get_step_prop('dt')
        low = np.min(ts)
        high = np.max(ts)
        avg = np.mean(ts)
        return low, high, avg

    def addOtp(self, filename, init=False):
        otpls = getLines(filename, ') |  ')
        dtnames = otpls[0].split('|')[-1]
        dtnames = dtnames.split()
        stags = ['slowx', 'slowy', 'slowz']
        tags = stags + dtnames
        ntags = len(tags) + 3
        
        all = otpls[1:]
        lines = len(all)
        block = ' '.join(all)
        chars = '()|,'
        for c in chars:
            block = block.replace(c, '')
        rdat = np.fromstring(block, sep=' ')
        rdat = rdat.reshape(lines, ntags)

        if init:
            for t in tags:
                stub = [0]*len(self.steps.index)
                self.steps = self.steps.assign(**{t: pd.Series(stub).values})
        for i, t in enumerate(tags):
            # log - otp offset is 1 step
            self.steps.loc[rdat[:, 0]-1, t] = rdat[:, 3 + i]


    def quick_look(self, retlist=False, refsimtime=0.1,
                  refstep=100, rsets=42, rounding=8):
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
        try:
            nodes = int(self.header.split('\n')[2].split()[-1])
        except AttributeError:  # workaround for different simuilation types
            return "simulation.quick_look: quick_look unavailable for this simulation."
        info.insert(len(info), 'Nodes used: {:>19}'.format(int(nodes/rsets)))
        info.insert(len(info), 'MPI Ranks: {:>20}'.format(nodes))
        info.append('Achieved timesteps: {:>11}'.format(len(self.steps)))
        tmax = self.get_step_prop('t')[-1] + self.get_step_prop('dt')[-1]
        info.append('Achieved simtime: {:>13}'.format(round(tmax, rounding)))
        info.append('Tstep sizes (s) [min/max(mean)]: '
                    '{:e}/{:e} ({:>10e})'.format(*self.get_avg_tstep()))
        info.append('Total runtime (hh:mm:ss): {}'.format(str(self.runtime)))
        info.append('IRL timestep (hh:mm:ss): {}'.format(self.irlstep))
        limblks = np.max(self.get_step_prop('totblocks', src='refs'))
        lowb = np.max(self.get_step_prop('minblocks', src='refs'))
        topb = np.max(self.get_step_prop('maxblocks', src='refs'))
        info.append('Max blocks used (min/max/per node): '
                    '{} ({}/{}/{:.2f})'.format(limblks, lowb, topb, limblks/(nodes/rsets)))
        info.append('Checkpoints/Plotfiles: '
                    '{}/{}'.format(len(self.checkpoints), len(self.plotfiles)))
        # finally read meta (if any)
        try:
            cwd = os.path.dirname(self.chk)
            with open(os.path.join(cwd, _metaname)) as f:
                metalines = [l.strip('\n') for l in f]
            info += metalines
        except OSError:
            info.append('No {} for this run.'.format(_metaname))
        if retlist:
            return info
        else:
            return '\n'.join(info)

    def get_t_fom(self, refsimt, tol=1e-3):
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
        ts = self.get_step_prop('dt')
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

    def get_n_fom(self, refstep):
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
        ts = self.get_step_prop('dt')
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


def read_log_and_stats(logfile, statsfile):
    """reads both log and stats file,
    correlating the data to each simulation step.
    """
    # read log data
    rsteps, chkps, header, timings = read_log(logfile)
    log.debug("finished reading {}".format(logfile))
    log.debug("Steps from f.log: {}".format(len(rsteps)))
    # read stats data and remove pre-restart values
    try:
        data = np.genfromtxt(statsfile, names=True)
    except OSError:
        print("no intergral quantities file.")
        log.debug("no s.dat file: skipping comparison")
        return rsteps, chkps, header, timings
    _, rep = np.unique(data['time'], axis=0, return_index=True)
    data = data[rep]
    log.debug("Steps from s.dat: {}".format(data.shape[0]))
    log.debug("Joining .log and .dat")
    diff = len(rsteps) - data.shape[0]
    if diff >= 1:
        # more steps in log than in stats, fill with zeros
        log.debug(".log > .dat, filling with zeros.")
        for name in data.dtype.names:
            stub = np.append(data[name], [0.0]*diff)
            rsteps = rsteps.assign(**{name: pd.Series(stub).values})
    else:
        # more steps in stats than in log
        # this is OK when diff is exactly 1:
        # stats writes summary final step which is not in log
        log.warning(".log < .dat, triming .dat.")
        for name in data.dtype.names:
            stub = data[name][:len(rsteps)]
            rsteps = rsteps.assign(**{name: pd.Series(stub).values})
    return rsteps, chkps, header, timings


def read_log(logfile):
    """read a flash log file grepping the first header,
    listing available timings
    and associating refinements to steps.
    The method used here takes into account for
    repeated timesteps (which happpen after restarts)
    but not for repeated refinements which is likely for failing
    runs but not dramatic for the delta in blocks.

    Args:
        logfile(str): filepath to *.log file from flash run.

    Return:
        (pd.DataFrame): every numbered step in the run (no repeated steps).
        (int list): removed step indices.
        (sim.steps): refinement step list.
        (sim.steps): checkpoint step list.
        (str): first header found in log file.
        (str list): all timing blocks found in log file.

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
    all = getLines(logfile, 'step: ')
    lines = len(all)
    chnk = [b.split(']')[0] for b in all]
    tstamps = []
    for b in chnk:
        c = b.strip('[')
        ts = datetime.datetime.strptime(c, ' %m-%d-%Y %H:%M:%S.%f ')
        tstamps.append(ts)
    chnk = [b.split('step: ')[1] for b in all]
    block = ' '.join(chnk)
    chars = 'n=td'
    for c in chars:
        block = block.replace(c, '')
    rdat = np.fromstring(block, sep=' ')
    rdat = rdat.reshape(lines, 3)
    steps = pd.DataFrame(data=rdat,
                         index=pd.RangeIndex(lines),
                         columns=['n', 't', 'dt'])
    steps = steps.assign(timestamp=pd.Series(tstamps).values)
    log.debug(".read_log: initial read steps: {}".format(len(steps)))
    steps, removedrows = clearRestarts(steps)
    log.debug(".read_log: cleared steps: {}".format(len(steps)))
    # add mpi ranks and the number of each restart
    cues, mpis = restartBreakdown(logfile)
    for i, s in enumerate(cues):
        # first submit, change all steps to catch first submit range blindly
        if not i:
            stub = [1]*len(steps.index)
            steps = steps.assign(submit_number=pd.Series(stub).values)
            stub = [mpis[0]]*len(steps.index)
            steps = steps.assign(mpi_ranks=pd.Series(stub).values)
            continue
        mask = steps['timestamp'] > s
        steps.loc[mask, 'submit_number'] = i + 1
        steps.loc[mask, 'mpi_ranks'] = mpis[i]
    log.debug('.read_log: adding refinement')
    refinementLines = getLines(logfile, '[GRID amr_refine_derefine]')
    # particles maim the output to 2 lines instead of 5
    # this can be picked by an extra refinement complete line at end
    subs = refinementLines[:10]
    stpchk = [l for l in subs if 'refinement complete' in l]
    if len(stpchk) > 0:
        stp = 5
    else:
        stp = 2
    tags = ['minblocks', 'maxblocks', 'totblocks']
    nsteps = len(steps.index)
    for i, rfb in enumerate(blockGenerator(refinementLines, step=stp)):
        tstamp, vals = refBreakdown(rfb)
        if not i:
            ddict = {}
            for j, t in enumerate(tags):
                ddict[t] = [vals[j]]*nsteps
        mask = steps['timestamp'] > tstamp
        try:
            loc = np.where(mask==True)[0][0]
        except IndexError:
            otag = '.read_log: refBreakdown mismatch loc {} vs {} steps'
            log.debug(otag.format(loc, nsteps))
        for j, t in enumerate(tags):
            ddict[t][loc:] = [vals[j]]*(nsteps-loc)
    for t in tags:
        steps = steps.assign(**{t: pd.Series(ddict[t]).values})
    log.debug('.read_log: adding checkpoint write-times')
    # checkpoint/plotfile timestamps
    # this varies between F versions and dimensionality.
    # assuming 2 lines with tag per checkpoint
    chkLs = getLines(logfile, "[IO_writeCheckpoint")
    chkLs += getLines(logfile, "[IO_writePlotfile")
    chkdata = [chkBreakdown(rfb) for rfb in blockGenerator(chkLs, step=2)]
    chkp = []
    for tst, (typ, num) in chkdata:
        chst = tstep()
        setattr(chst, 'timestamp', tst)
        setattr(chst, 'filetype', typ)
        setattr(chst, 'number', num)
        chkp.append(chst)
    return steps, chkp, header, timings


def getIRLdeltas(steps):
    """calculate walltime deltas between timesteps
    accounting for resubmits by counting submits.

    Args:
        steps(sim.step list): simulation log steps.

    Returns:
        (timedelta, timedelta): total walltime, mean walltime step

    """
    tstamps = steps['timestamp']
    submits = steps['submit_number']
    deltas = np.array([b-a for (a, b) in zip(tstamps[:-1], tstamps[1:])])
    subdelts = np.diff(submits)
    subdelts = np.append(subdelts, 0)
    ind = np.where(subdelts == 1)[0]
    mask = np.zeros(len(deltas), dtype=bool)
    mask[ind] = True
    # stub dt is the mean delta for non-issue steps
    stubdt = np.mean(deltas[~mask])
    deltas[mask] = stubdt
    deltas = np.append(deltas, stubdt)
    return deltas


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
        (pd.DataFrame): list of simulation.step objects.

    Returns:
        (pd.DataFrame): cleansed steps.
        (int list): removed step indices.

    """
    numbers = np.array(steps.index)
    totnstp = len(steps)
    checkdelt = np.diff(numbers - totnstp)  # find restarts
    skipmask = np.append(checkdelt, 1)
    restartpos = np.where(checkdelt != 1)[0]  # get restart positions
    repsteps = abs(checkdelt[restartpos])+1  # get wasted steps

    if repsteps.size > 0:
        remsteps = []
        for pos, reps in zip(restartpos, repsteps):
            rr = int(reps)
            skipmask[pos-rr+1:pos+1] = 5  # any value greater than one
            for i in range(rr):
                remsteps.append(pos-i)
        realmask = np.where(skipmask == 1)[0]
        steps = steps.drop(realmask)
        steps.index = pd.RangeIndex(len(steps.index))
        return steps, remsteps
    else:
        return steps, []


def restartBreakdown(logfile):
    """read starting line in each run/restart
    header for MPI ranks and starting time.

    Args:
        logfile(str): path to log file

    Returns:
        (datetime list): starting times.
        (int list): mpi ranks.

    """
    lines = getLines(logfile, 'FLASH log file')
    startfmt = 'FLASH log file:  %m-%d-%Y  %H:%M:%S.%f    '
    mpistrs = getLines(logfile, 'Number of MPI tasks')
    cues, mpis = [], []
    for input, mpistr in zip(lines, mpistrs):
        cue = datetime.datetime.strptime(input.split('Run')[0], startfmt)
        mpinum = int(mpistr.split()[-1])
        cues.append(cue)
        mpis.append(mpinum)
    return cues, mpis


def chkBreakdown(refBlock):
    """specific breakdown of a checkpoint/plotfile write log entry.
    e.g.:
    '[ 05-29-2020  06:28:27.525 ] [IO_writeCheckpoint] open: type=checkpoint
        name=../17ctm_p4248_90_ms32_t3.0/chk/flash_hdf5_chk_0132
     [ 05-29-2020  06:28:31.937 ] [IO_writeCheckpoint] close: type=checkpoint
        name=../17ctm_p4248_90_ms32_t3.0/chk/flash_hdf5_chk_0132'

    particles add a line:
    '[ 10-19-2020  10:58:01.960 ] [IO_writeCheckpoint] open: type=checkpoint
        name=../00ctm_p4248_90_ms32_t3.0/chk/flash_hdf5_chk_0011
     [ 10-19-2020  10:58:02.210 ] [IO_writeParticles]: done called Particles_updateAttributes()
     [ 10-19-2020  10:58:02.216 ] [IO_writeCheckpoint] close: type=checkpoint
         name=../00ctm_p4248_90_ms32_t3.0/chk/flash_hdf5_chk_0011

    """
    # get tstamp and type from first line
    tstampstr, _, rest = refBlock[0].partition(']')
    tstamp = datetime.datetime.strptime(tstampstr, '[ %m-%d-%Y %H:%M:%S.%f ')
    _, _, ftype = rest.partition('=')
    # get number from first line too due to changing text block size
    try:
        number = int(refBlock[0][-4:])
    except ValueError:  # sometimes name of file is not in first line 
        number = int(refBlock[1][-4:])
    return tstamp, (ftype, number)


def refBreakdown(refBlock):
    blockstencil = [2, 5, 8]
    # get tstamp from first line
    try:
        tstampstr, _, _ = refBlock[0].partition(']')
    except IndexError:
        log.debug("ref block timestamp fail: \n{}".format('\n'.join(refBlock)))
        tstampstr = '[ 01-01-2000  20:20:20.200 '
    tstamp = datetime.datetime.strptime(tstampstr, '[ %m-%d-%Y %H:%M:%S.%f ')
    # block information from line 2 or 3 (depending on input from log)
    if len(refBlock) > 2:
        blkline = 2
    else:
        blkline = 1
    try:
        _, _, blockstr = refBlock[blkline].partition(']')
        min, max, tot = [int(x) for i, x in enumerate(blockstr.split())
                         if i in blockstencil]
    # if anything fails set to a fixed number
    except IndexError:
        log.debug("ref block numbers fail: \n{}".format('\n'.join(refBlock)))
        min, max, tot = 2, 2, 2
    except ValueError:
        log.debug("ref block numbers fail: \n{}".format('\n'.join(refBlock)))
        min, max, tot = 2, 2, 2
    return tstamp, [min, max, tot]


def refBreakdown_old(refBlock):
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
        (datetime.datetime): timestamp
        (6 floats): min/max/tot blocks base and leaf.

    """
    blockstencil = [2, 5, 8]
    leafbstencil = [3, 7, 11]
    # get tstamp from first line
    try:
        tstampstr, _, _ = refBlock[0].partition(']')
    except IndexError:
        log.debug("ref block timestamp fail: \n{}".format('\n'.join(refBlock)))
        tstampstr = '[ 01-01-2000  20:20:20.200 '
    tstamp = datetime.datetime.strptime(tstampstr, '[ %m-%d-%Y %H:%M:%S.%f ')
    # block information from lines 3-4
    try:
        _, _, blockstr = refBlock[2].partition(']')
        min, max, tot = [int(x) for i, x in enumerate(blockstr.split())
                         if i in blockstencil]
        _, _, blockstr = refBlock[3].partition(']')
        lmin, lmax, ltot = [int(x) for i, x in enumerate(blockstr.split())
                            if i in leafbstencil]
    except IndexError:
        log.debug("ref block numbers fail: \n{}".format('\n'.join(refBlock)))
        min, max, tot, lmin, lmax, ltot = 1, 2, 3, 4, 5, 6
    except ValueError:
        log.debug("ref block numbers fail: \n{}".format('\n'.join(refBlock)))
        min, max, tot, lmin, lmax, ltot = 1, 2, 3, 4, 5, 6
    return tstamp, [min, max, tot, lmin, lmax, ltot]


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
