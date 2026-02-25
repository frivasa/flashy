import os
import sys
import shutil
import logging
from .utils import np
from subprocess import PIPE, Popen
# import imageio
import itertools

_cdxfolder = "cdx"
_otpfolder = "chk"
_logfile = 'f.log'
_statsfile = 's.dat'
_cdxpfol = "cdx/Profiles"
_metaname = "meta.txt"
_deftree = 'FLASH4'
_proj = 'gen006'  # proj where tree is hosted
_juptree = 'http://localhost:1988/tree/'
# summit has POWER arch
# if os.getenv('HOSTTYPE', 'powerpc64le') == 'powerpc64le':
_FLASH_DIR = "/gpfs/alpine/{}/proj-shared/"\
             "frivas/00.code/{}".format(_proj, _deftree)
_AUX_DIR = "/gpfs/alpine/{}/proj-shared/"\
           "frivas".format(_proj)
_SCRIPT_DIR = "/gpfs/alpine/{}/proj-shared/"\
              "frivas/07.miscellaneous/bash/".format(_proj)
# restore maxint number
# (removed in p3 due to arbitrary int length, but FORTRAN doesn't know.)
_maxint = 2147483647
log = logging.getLogger(__name__)
handler = logging.StreamHandler(stream=sys.stdout)
logstr = "[%(levelname)-8s] %(name)s: %(message)s"
formatter = logging.Formatter(logstr)
handler.setFormatter(formatter)
log.addHandler(handler)
log.setLevel(20)
# (50, 10, -10): critical, error, warn, info, debug


def setupFLASH(module, runfolder='', kwargs={'threadBlockList': 'true'},
               nbs=[16, 16, 16], geometry='cylindrical',
               maxbl=500, codetree='FLASHOR'):
    """prints adequate ./setup at _FLASH_DIR with given parameters,
    writing the code to runfolder. (FLASH setup script runs only on py 2.X).

    Arguments:
        module(str): name of Simulation folder to setup.
        runfolder(str): run folder (creates _cdx within it).
        kwargs(dict): keyword arguments to setup
        nbs(int tuple): cells per block for setup per dimension.
        geometry(str): cartesian, spherical, cylindrical(default).
        maxbl(int): maximum blocks per proc. elem.
        codetree(string): force flash ver. (e.g. FLASHOR > FLASH5).

    """
    if not runfolder:
        # theres only ap13, ap19 or xnet so treat by cases
        if '+xnet' in kwargs.keys():
            net = kwargs['xnetData'].split('_')[-1]
        elif '+a13' in kwargs.keys():
            net = 'ap13'
        elif '-with-unit' in kwargs.keys():
            argu = kwargs['-with-unit']
            net = 'ap{}'.format(int(argu[-2:]))
        # first 'three' without a 't' based on # speakers is
        # Mandarin -> German -> Japanese
        dnum = {1: 'one', 2: 'two', 3: 'san'}[len(nbs)]
        cells = 'x'.join([str(x) for x in nbs])
        name = '{}.{}.{}cells.{}maxb'.format(net, geometry, cells, maxbl)
        # TODO: remove/patch relative paths
        relpath = 'fruns.{}/{}/{}'.format(dnum, module, name)
        runfolder = os.path.join(_AUX_DIR, relpath)
    else:
        name = runfolder
    destination = os.path.join(runfolder, _cdxfolder)
    path = os.path.abspath(destination)
    if not os.path.exists(destination):
        os.makedirs(destination)
    else:
        print('Emptying {}'.format(destination))
        shutil.rmtree(destination)
        os.makedirs(destination)

    if len(nbs) == 1:
        dimstr = "-1d -nxb={}".format(*nbs)
    elif len(nbs) == 2:
        dimstr = "-2d -nxb={} -nyb={}".format(*nbs)
    else:
        dimstr = "-3d -nxb={} -nyb={} -nzb={}".format(*nbs)

    cstub = 'cd {} \n ./setup {} -auto '\
            '-objdir="{}" {} -geometry={} -maxblocks={} '
    comm = cstub.format(_FLASH_DIR, module, path,
                        dimstr, geometry, maxbl)
    kwstr = ''
    for (k, v) in kwargs.items():
        if v == '':
            kwstr += ' {}'.format(k)
        else:
            kwstr += ' {}={}'.format(k, v)
    # kwstr = ' '.join(['{}={}'.format(k,v) for (k, v) in kwargs.items()])
    comm = comm + kwstr
    print('\nGenerated run name: {}\n'.format(name))
    comm = comm.replace(_deftree, codetree)
    print(comm)
    # p = Popen(['/bin/bash'], stdin=PIPE, stdout=PIPE)
    # out, err = p.communicate(input=comm.encode())
    # exitcode = p.returncode
    # copy parameter varying script
    try:
        iterpath = os.path.join(_AUX_DIR, '07.miscellaneous/bash/iterator')
        print('cp {} {}/.'.format(iterpath, destination))
    except:
        print('bash iterator not found, skipping.')

    print("************************************************************"
          "************************************************************"
          "\nFLASH uses a preprocessing bash script to set an env_var so "
          "popen doesn't work. Please run the above commands by hand.\n"
          "************************************************************"
          "************************************************************")
    return comm  # , exitcode, out, err


def mergeRuns(runlist, newname):
    """builds a newname runfolder with aggregate properties from
    each run in runlist.

    Args:
        runlist(str list): paths of runs to join.
        newname(str): pathname of aggregated run.

    """

    sortinp = sorted([i.strip('/') for i in runlist])
    runfolder = os.path.abspath(newname)
    print("Creating:", runfolder)
    os.makedirs(os.path.join(runfolder, _otpfolder))
    # first move all checkpoints to the new chk folder
    dst = os.path.join(os.path.abspath(runfolder), _otpfolder)
    for r in sortinp:
        src = os.path.join(os.path.abspath(r), _otpfolder)
        for f in os.listdir(src):
            shutil.copy2(os.path.join(src, f), dst)
    # next copy all non-appendable output
    dst = os.path.abspath(runfolder)
    for r in sortinp:
        src = os.path.abspath(r)
        srcf = os.listdir(src)
        files = [f for f in srcf if os.path.isfile(os.path.join(src, f))]
        nonappfiles = [f for f in files if f != _logfile and f != _statsfile]
        for f in nonappfiles:
            if f == 'flash.par':
                parname = os.path.basename(r)
                shutil.copy2(os.path.join(src, f),
                             os.path.join(dst, "{}.par".format(parname)))
            else:
                shutil.copy2(os.path.join(src, f), dst)
    # finally append log and stat file
    dst = os.path.abspath(runfolder)
    refstats = os.path.join(dst, _statsfile)
    reflog = os.path.join(dst, _logfile)
    for r in sortinp:
        src = os.path.abspath(r)
        stats = os.path.join(src, _statsfile)
        log = os.path.join(src, _logfile)
        os.system("cat " + stats + " >> " + refstats)
        os.system("cat " + log + " >> " + reflog)


def getFileList(folder, glob='plt', fullpath=False):
    """Returns a filename list subject to a prefix 'glob'.

    Args:
        folder(str): look-in path.
        prefix(str): filter string for files
        fullpath(bool): return absolute path for each file.

    Returns:
        (str list)

    """
    names = sorted(os.listdir(folder))
    fnames = [x for x in names if glob in x]
    if fullpath:
        return [os.path.join(os.path.abspath(folder), x) for x in fnames]
    else:
        return fnames


def rename(path, glob, newname):
    """rename file subject to glob to newname."""
    templog = getFileList(path, glob=glob, fullpath=True)
    if templog:
        os.rename(templog[0],
                  os.path.join(os.path.split(templog[0])[0], newname))
        print('IOutils.rename: changed {} to {}'.format(templog, newname))


def getLines(filename, keyword):
    """returns all lines matching a keyword within them"""
    lines = []
    with open(filename, 'r') as f:
        for l in f:
            if keyword in l:
                lines.append(l.strip())
    return lines


def getBlock(filename, initkey, endkey, skip=0):
    """returns all lines within keywords, clearing empty ones.

    Args:
        filename(str): input file.
        initkey(str): initial block keyword.
        endkey(str): ending block keyword.
        skip(int): skip the initkey skip times before extracting.

    Returns:
        (str list): block extracted without blank lines.

    """
    lines = []
    save = False
    par = 0
    with open(filename, 'r') as f:
        for l in f:
            if initkey in l:
                par += 1
                if par > skip:
                    save = True
            elif save and endkey in l:
                break
            if save and len(l.strip(' \n')) > 0:
                lines.append(l.strip('\n'))
    return list(filter(None, lines))


def getRepeatedBlocks(filename, initkey, endkey):
    """return every key bracketed block of text from a file."""
    blocks = []
    while(True):
        block = getBlock(filename, initkey, endkey, skip=len(blocks))
        if block:
            blocks.append(block)
        else:
            break
    return blocks


def blockGenerator(lines, step=5):
    """generator for line lists of 'step' size."""
    for i in range(0, len(lines), step):
        yield lines[i:i+step]


def pairGen(objlist):
    """ turns a list into subsequent pairs between the objects """
    a, b = itertools.tee(objlist)
    next(b, None)
    return zip(a, b)


def fortParse(arg, dec=True):
    """returns a parsed variable from a parameter (bool,
    str, or number)

    Args:
        arg(str): parameter value.
        dec(bool): add "" to strings for printing.

    Returns:
        str: decorated argument for fortran parsing.

    """
    booleans = ['false', 'true']
    query = arg.strip('."\' ').lower()
    try:
        val = np.float(query.replace('d', 'E'))
        # FORT does not handle big ints nor 'X.XEXX ints.
        if int(val) == val:
            if abs(val) > _maxint:
                return '{:E}'.format(val)
            else:
                return int(val)
        else:
            return '{:E}'.format(val)
        return '{:E}'.format(val)
    except ValueError:
        if query in booleans:
            return '.{}.'.format(query)
        else:
            if dec:
                return '"{}"'.format(arg.strip('"\' '))
            else:
                return arg.strip('"\' ')


def execute(outpath):
    """qsubs the sumbit.pbs at outpath

    Args:
        outpath (str): runfolder

    Returns:
        (tuple): STDOUT, STDERR, ERRCODE

    """
    command = 'qsub submit.pbs'
    p = Popen(command.split(), cwd=os.path.abspath(outpath),
              stdin=PIPE, stdout=PIPE, stderr=PIPE)
    r, e = p.communicate()
    exitcode = p.returncode
    return r, e, exitcode


def cpList(files, src, dst):
    """Copies a file list between folders.

    Args:
        files(str list): list of filenames.
        src(str): source folder.
        dst(str): destination folder.

    """
    for f in files:
        shutil.copy('/'.join([src, f]), '/'.join([dst, f]))


def getMachineSubJob(machine, proj, time, nodes, ompth, otpf, **kwargs):
    """returns scheduler code and extension for submit file according
    to the specified queued machine.

    Args:
        machine(str): machine name.
        proj(str): project allocation code.
        time(str): walltime request.
        nodes(int): nodes to request.
        ompth(int): omp thread number.
        otpf(str): output filename.
        **kwargs: additional keywords for ad-hoc changes.

    """
    if machine == 'summit':
        return summitBatch(proj, time, nodes, otpf, **kwargs)
    else:
        return titanBatch(proj, time, nodes, ompth, otpf, **kwargs)


def titanBatch(proj, time, nodes, ompth, otpf, **kwargs):
    """Titan/PBS submit maker.
    builds a submit.pbs with a typical header, specifying walltime and nodes.
    rhea: mpirun --map-by ppr:N:node:pe=Th or -n
    debug: --display-map / --report-bindings
    Rhea max: 48 hours on 16 nodes (2x8 core p/node: -np 256)
    Titan: <125 nodes 2h, <312 nodes 6h...

    Args:
        proj(str): project allocation code.
        time(str): walltime request.
        nodes(int): nodes to request.
        ompth(int): omp thread number.
        otpf(str): output filename.
        **kwargs: additional keywords for ad-hoc changes.

    """

    schedlist = []
    schedlist.append('#!/bin/bash')
    schedlist.append('#PBS -V')  # pass env vars to nodes
    schedlist.append('#PBS -l gres=atlas1')
    schedlist.append('#PBS -A {}'.format(proj))
    schedlist.append('#PBS -l walltime={},nodes={}'.format(time, nodes))
    schedlist.append('#PBS -N {}'.format(os.path.basename(otpf)))
    schedlist.append('#PBS -j oe')  # join err and otp
    schedlist.append('#PBS -o {}'.format(otpf))
    if kwargs.get('j1', False):
        nodes *= 2
        ompth = 8
        launcher = 'aprun -n{} -d{} -j1 ./flash4 &'.format(nodes, ompth)
    else:
        launcher = 'aprun -n{} -d{} ./flash4 &'.format(nodes, ompth)
    if 'mail' in kwargs:
        schedlist.append('#PBS -M {}'.format(kwargs['mail']))
        schedlist.append('#PBS -m {}'.format(kwargs['abe']))
    schedlist.append('export OMP_NUM_THREADS={}'.format(int(ompth)))
    schedlist.append('export CRAY_CUDA_MPS=1')
    return launcher, 'qsub', '.pbs', schedlist


def summitBatch(proj, time, nodes, otpf,
                defs = {'RSpN': 2, 'mpipRS': 21}, **kwargs):
    """Summit/BSUB submit maker.
    builds a submit.lsf with a typical header, specifying walltime and nodes.

    Args:
        proj(str): project allocation code.
        time(str): walltime request.
        nodes(int): nodes to request.
        smt(int): omp/smt thread number.
        otpf(str): output filename.
        defs(dict): RSpN (RS per node) and mpipRS (ranks per RS: -a #).
        **kwargs: additional keywords for ad-hoc changes.

    """
    Ncores = 42
    Ngpu = 6
    cpRS = int(Ncores/defs['RSpN'])
    cpmpi = int(cpRS/defs['mpipRS'])
    gpupRS = int(Ngpu/defs['RSpN'])
    RS = nodes*defs['RSpN']
    schedlist = []
    # translate to BSUB
    schedlist.append('#!/bin/bash')
    schedlist.append('#BSUB -env all')  # pass env vars to nodes
    # schedlist.append('#BSUB -l gres=atlas1')
    schedlist.append('#BSUB -P {}'.format(proj))
    schedlist.append('#BSUB -nnodes {}'.format(nodes))
    schedlist.append('#BSUB -W {}'.format(time[:-3]))  # hh:mm without :ss
    schedlist.append('#BSUB -J {}.lsf'.format(os.path.basename(otpf)))
    if 'smt' in kwargs:
        ompth = cpRS/defs['mpipRS']*kwargs['smt']
        schedlist.append('#BSUB -alloc_flags "gpumps smt{}"'.format(smt))
    else:
        ompth = cpRS/defs['mpipRS']*4
        schedlist.append('#BSUB -alloc_flags "gpumps smt4"')
    schedlist.append('#BSUB -N')
    # schedlist.append('#BSUB -outdir {}'.format(otpf))  # LSB_OUTDIR
    schedlist.append('#BSUB -o {}/{}.o%J'.format(otpf, os.path.basename(otpf)))
    schedlist.append('#BSUB -e {}/{}.e%J'.format(otpf, os.path.basename(otpf)))
    if 'mail' in kwargs:
        schedlist.append('#BSUB -u {}'.format(kwargs['mail']))
    else:
        schedlist.append('#BSUB -u rivas.aguilera@gmail.com')
    # count nodes from bash
    # NNODES=$(($(cat $LSB_DJOB_HOSTFILE | uniq | wc -l)-1))
    # 2 node 7 core + GPU per task 4 threads
    # launcher = 'stdbuf -o0 jsrun '
    # '--exit-on-error -n12 -g1 -a1 -c7 -bpacked:7 ./flash4 &'
    launcher = 'jsrun '
    launcher += '-n{} -r{} -a{} '.format(RS, defs['RSpN'], defs['mpipRS'])
    launcher += '-g{} -c{} -bpacked:{} '.format(gpupRS, cpRS, cpmpi)
    launcher += '-d packed ./flash4 &'
    # -n = "Resource set"[rs] as subdivisions of a node
    # -a  MPI ranks/tasks per rs
    # -c cpus per rs (physical, non-threaded)
    # -g gpus per rs
    # -b bind withing an rs [none, rs or packed number]
    # -r rs per host=node
    # -l latency priority (cpu-gpu gpu-cpu)
    # -d launch distribution (task starting order)
    # pre 2021 08 20: hdf5/1.8.18 or 1.8.22
    if kwargs['compiler']=='gcc':
        schedlist.append('module load gcc cuda essl netlib-lapack hdf5')
        schedlist.append('module use /sw/summit/unifyfs/modulefiles')
        schedlist.append('module load unifyfs/0.9.2-mpi-mount')
        schedlist.append('export UNIFYFS_LOGIO_SPILL_SIZE=10000000000')
        schedlist.append('export UNIFYFS_LOG_VERBOSITY=5')
    else:
        schedlist.append('module load pgi cuda essl netlib-lapack hdf5')
    # write specific romio_ hints
    romiofile = writeRHints(nodes)
    schedlist.append('export ROMIO_HINTS={}'.format(romiofile))
    schedlist.append('export OMP_NUM_THREADS={}'.format(int(ompth)))
    schedlist.append('export OMP_SCHEDULE="dynamic"')
    schedlist.append('export OMP_STACKSIZE="1G"')  # def is 512
    schedlist.append('export '
                     'OMPI_LD_PRELOAD_POSTPEND='
                     '${OLCF_SPECTRUM_MPI_ROOT}/lib/libmpitrace.so')
    
    # warning this makes it extremely slow and faulty
    # profiler: export OMPI_LD_PRELOAD_POSTPEND=
    # /ccs/home/walkup/mpitrace/spectrum_mpi/libhpmprof.so
    return launcher, 'bsub', '.lsf', schedlist


def writeRHints(nodes, buffer=33554432, multi=4):
    """write multithread writing bash file."""
    name = 'romio_h{}'.format(multi*nodes)
    romiofile = os.path.join(_SCRIPT_DIR, name)
    with open(romiofile, 'w') as o:
        o.write('romio_cb_write enable\n')
        o.write('romio_ds_write disable\n')
        o.write('cb_buffer_size {:d}\n'.format(int(buffer)))
        o.write('cb_nodes {:d}\n'.format(multi*nodes))
    return romiofile


def writeSchedulerScript(subfile, code, schedcode):
    """writes submit file.

    Args:
        subfile(str): submit filename (also set as jobname)
        code(str list): ordered commands for the file.
        schedcode(str list): scheduler directives.

    """
    with open(subfile, 'w') as o:
        o.write("\n".join(schedcode))
        o.write("\n")
        o.write("\n".join(code))
        o.write("\n")


def getWalltime(nodes, machine='summit'):
    """get max walltime based on nodes requested."""
    time = globals()[machine](nodes)
    return time


def summit(nodes):
    """lsf uses hh:mm not hh:mm:ss but this is
    handled by the lsf submit writer.
    """
    if nodes < 46:
        return '02:00:00'
    elif nodes < 92:
        return '06:00:00'
    elif nodes < 922:
        return '12:00:00'
    else:
        return '24:00:00'


def titan(nodes):
    """return walltime for machine: titan"""
    if nodes < 125:
        return '02:00:00'
    elif nodes < 312:
        return '06:00:00'
    elif nodes < 3749:
        return '12:00:00'
    else:
        return '24:00:00'


def probeFile(file, showrows=3, onlyhead=True):
    """Shows 'showrows' lines from the start, midfile and
    ending of a plaintext file

    Args:
        file(str): file path.
        showrows(int): rows to show from each section.
        onlyhead(bool): only print the top of the file
        (equivalent to head -n showrows file).

    """
    with open(file, 'r') as f:
        lines = f.readlines()
    l = len(lines)
    if l < 10:
        print("".join(lines))
    else:
        l2 = int(l/2.)
        pNumbered(lines[0:showrows], 0)
        if not onlyhead:
            pNumbered(lines[l2:l2+showrows], l2)
            pNumbered(lines[-showrows:], l-showrows)


def pNumbered(rlist, offset):
    """prints numbered lines from a list.

    Args:
        rlist(str list): lines to print
        offset(int): first line number.

    """
    for i, l in enumerate(rlist):
        print('{}: {}'.format(i+offset, l.strip('\n')))


def emptyFileTree(stemfolder):
    """Empties stemfolder."""
    path = os.path.abspath(root)
    shutil.rmtree(path)
    os.makedirs(path)


def setFolders(fpath, filetag):
    """Builds folders according to a given filepath.
    Goes back a folder and creates a new one named 'filetag'
    "/path/to/file/_chkfolder/checkpoint_0001"
    > num = 0001
    > dest = /path/to/file/filetag/
    > name = dest+filetag_num

    Args:
        fpath(str): filepath.
        filetag(str): prefix for output file.

    Returns:
        (str): destination path of the file.
        (str): file suffix number.
        (str): output filename with no extension.

    """
    num = fpath[-4:]  # checkpoint number 'flash_hdf5_chk_0001'
    basedest = os.path.dirname(os.path.dirname(fpath))
    dest = os.path.join(basedest, filetag)
    name = os.path.join(dest, '{}_{}'.format(filetag, num))
    os.makedirs(dest, exist_ok=True)  # bless you, p3
    return dest, int(num), name

def argBlock(fields, space, columns):
    """splits a list into an evenly spaced string grid."""
    grid = (space, columns)
    remnant = len(fields)%grid[1]
    fit = iter(fields)
    line = "{{:{:d}}}".format(grid[0])*grid[1]
    blk = []
    for tpl in zip(*[fit]*grid[1]):
        blk.append(line.format(*tpl))
    blk.append((f"{{:{grid[0]}}}"*remnant).format(*fields[-remnant:]))
    return "\n".join(blk)
