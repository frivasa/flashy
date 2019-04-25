import os
import sys
import shutil
import collections as cl
from .utils import np
from subprocess import PIPE, Popen
import subprocess as subp
import operator
from PIL import Image
import imageio

_cdxfolder = "cdx"
_otpfolder = "chk"
_logfile = 'f.log'
_statsfile = 's.dat'
_juptree = 'http://localhost:1988/tree/'
if os.getenv('HOSTTYPE', 'powerpc64le')=='powerpc64le':  # summit has POWER arch
    _FLASH_DIR = "/gpfs/alpine/csc198/proj-shared/frivas/00.code/FLASHOR"  # summit testing
    _AUX_DIR = "/gpfs/alpine/csc198/proj-shared/frivas"
    _SCRIPT_DIR = "/gpfs/alpine/csc198/proj-shared/frivas/07.miscellaneous/bash/"
else:
    _FLASH_DIR = "/lustre/atlas/proj-shared/csc198/frivas/00.code/FLASHOR"
    _AUX_DIR = "/lustre/atlas/proj-shared/csc198/frivas/"
    _SCRIPT_DIR = "/lustre/atlas/proj-shared/csc198/frivas/07.miscellaneous/bash/"
_maxint = 2147483647  # this was removed in p3 due to arbitrary int length, but FORTIE doesn't know...

def setupFLASH(module, runfolder='', kwargs={'threadBlockList':'true'}, nbs=[16, 16, 16],
               geometry='cylindrical', maxbl=500):
    """calls ./setup at _FLASH_DIR with given parameters,
    writing the code to runfolder. (FLASH setup script runs only on py 2.X).

    Arguments:
        module(str): name of Simulation folder to setup.
        runfolder(str): run folder (creates _cdx within it).
        kwargs(dict): keyword arguments to setup
        nbs(int tuple): cells per block for setup per dimension.
        geometry(str): cartesian, spherical, cylindrical(default).
        maxbl(int): maximum blocks per proc. elem.
        manual(bool): skip terminal call, print command.

    """
    if not runfolder:
        net = 'ap13'
        if 'xnet' in kwargs:
            if kwargs['xnet']==True:
                net = kwargs['xnetData'].split('_')[-1]
        dnum = {1:'one', 2:'two', 3:'three'}[len(nbs)]
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
    
    if len(nbs)==1:
        dimstr = "-1d -nxb={}".format(*nbs)
    elif len(nbs)==2:
        dimstr = "-2d -nxb={} -nyb={}".format(*nbs)
    else:
        dimstr = "-3d -nxb={} -nyb={} -nzb={}".format(*nbs)
    
    cstub = 'cd {} \n ./setup {} -auto '\
            '-objdir="{}" {} -geometry={} -maxblocks={} '
    comm = cstub.format(_FLASH_DIR, module, path,
                        dimstr, geometry, maxbl)
    kwstr = ''
    for (k, v) in kwargs.items():
        if v=='':
            kwstr+=' {}'.format(k)
        else:
            kwstr+=' {}={}'.format(k,v)
    #kwstr = ' '.join(['{}={}'.format(k,v) for (k, v) in kwargs.items()])
    comm = comm + kwstr
    print('\nGenerated run name: {}\n'.format(name))
    print(comm)
    #p = Popen(['/bin/bash'], stdin=PIPE, stdout=PIPE)
    #out, err = p.communicate(input=comm.encode())
    #exitcode = p.returncode
    
    # copy parameter varying script
    try:
        iterpath = os.path.join(_AUX_DIR, '07.miscellaneous/bash/iterator')
#         shutil.copy2(iterpath, destination)
#         print('copied iterator')
        print('cp {} {}/.'.format(iterpath, destination))
    except:
        print('bash iterator not found, skipping.')
    
    print("************************************************************"\
          "************************************************************"
          "\nFLASH uses a preprocessing bash script to set an env_var so "\
          "popen doesn't work. Please run the above commands by hand.\n"\
          "************************************************************"\
          "************************************************************")
    return comm  #, exitcode, out, err


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
            shutil.copy2(os.path.join(src,f), dst)
    # next copy all non-appendable output
    dst = os.path.abspath(runfolder)
    for r in sortinp:
        src = os.path.abspath(r)
        files = [f for f in os.listdir(src) if os.path.isfile(os.path.join(src, f))]
        nonappfiles = [f for f in files if f!=_logfile and f!=_statsfile]
        for f in nonappfiles:
            if f=='flash.par':
                parname = os.path.basename(r)
                shutil.copy2(os.path.join(src, f), os.path.join(dst, "{}.par".format(parname)))
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
    templog = getFileList(path, glob=glob, fullpath=True)
    if templog:
        os.rename(templog[0], os.path.join(os.path.split(templog[0])[0], newname))
        print('IOutils.rename: changed {} to {}'.format(templog, newname))


def makeGIF(srcfolder, speed=0.2):
    """Join all png images within a folder in an 
    animated .gif. Outputs at srcfolder/../
    # reduce size via webm conversion.

    Args:
        srcfolder (str): folder path
        speed (float): seconds between frames
    
    """
    finns = sorted([x for x in os.listdir(srcfolder) if '.png' in x])
    outfolder = os.path.dirname(srcfolder)
    # maim the first file to get a name for the gif
    outname = os.path.join(os.path.dirname(outfolder), '{}.gif'.format(finns[0][:-9]))
    jakes = []
    for finn in finns:
        sys.stdout.write(finn + " ")
        jakes.append(imageio.imread(os.path.join(srcfolder,finn)))
    imageio.mimsave(outname, jakes, format='gif', duration=speed)
    print("\n\tSaved: {}".format(outname))


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
        val = np.float(query.replace('d','E'))
        if int(val)==val: # FORT does not handle big ints nor 'X.XEXX ints.
            if abs(val)>_maxint:
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
        shutil.copy('/'.join([src,f]), '/'.join([dst,f]))


def getMachineSubJob(machine, proj, time, nodes, ompth, ddt, otpf, **kwargs):
    """returns scheduler code and extension for submit file according 
    to the specified queued machine.
    
    Args:
        machine(str): machine name.
        proj(str): project allocation code.
        time(str): walltime request.
        nodes(int): nodes to request.
        ompth(int): omp thread number.
        ddt(bool): enable arm forge connection.
        otpf(str): output filename.
        **kwargs: additional keywords for ad-hoc changes.
    
    """
    if machine=='summit':
        return summitBatch(proj, time, nodes, ddt, otpf, **kwargs)
    else:
        return titanBatch(proj, time, nodes, ompth, ddt, otpf, **kwargs)


def titanBatch(proj, time, nodes, ompth, ddt, otpf, **kwargs):
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
        ddt(bool): enable arm forge connection.
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
        nodes*=2
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


def summitBatch(proj, time, nodes, ddt, otpf, **kwargs):
    """Summit/BSUB submit maker.
    builds a submit.lsf with a typical header, specifying walltime and nodes.
    
    Args:
        proj(str): project allocation code.
        time(str): walltime request.
        nodes(int): nodes to request.
        smt(int): omp/smt thread number.
        ddt(bool): enable arm forge connection.
        otpf(str): output filename.
        **kwargs: additional keywords for ad-hoc changes.
    
    """
    Ncores = 42
    Ngpu = 6
    defs = {'RSpN': 6, 'mpipRS': 1}
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
    schedlist.append('#BSUB -J {}'.format(os.path.basename(otpf)))
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
    # count nodes from bash NNODES=$(($(cat $LSB_DJOB_HOSTFILE | uniq | wc -l)-1))
    # 2 node 7 core + GPU per task 4 threads
    #launcher = 'stdbuff -o0 jsrun --exit-on-error -n12 -g1 -a1 -c7 -bpacked:7 ./flash4 &'
    launcher = 'stdbuf -o0 jsrun '  # --exit_on_error > Must specify either --np or (--tasks_per_rs and --nrs).
    launcher += '-n{} -r{} -a{} -g{} -c{} -bpacked:{} '.format(RS, defs['RSpN'], defs['mpipRS'], gpupRS, cpRS, cpmpi) 
    launcher += '-d packed ./flash4 &' 
    #-n = "Resource set"[rs] as subdivisions of a node
    #-a  MPI ranks/tasks per rs
    #-c cpus per rs (physical, non-threaded)
    #-g gpus per rs
    #-b bind withing an rs [none, rs or packed number]
    #-r rs per host=node
    #-l latency priority (cpu-gpu gpu-cpu)
    #-d launch distribution (task starting order)

    # hdf5 parallel to serial hack
    schedlist.append('module load pgi cuda essl netlib-lapack hdf5/1.8.18')
    # write specific romio_ hints
    romiofile = writeRHints(nodes)
    schedlist.append('export ROMIO_HINTS={}'.format(romiofile))
    schedlist.append('export OMP_NUM_THREADS={}'.format(int(ompth)))
    schedlist.append('export OMP_SCHEDULE="dynamic"')
    schedlist.append('export OMP_STACKSIZE="256M"')
    return launcher, 'bsub', '.lsf', schedlist


def writeRHints(nodes, buffer=33554432, multi=4):
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
    """lsf uses hh:mm not hh:mm:ss but this is handled by the lsf submit writer."""
    if nodes<46:
        return '02:00:00'
    elif nodes<92:
        return '06:00:00'
    elif nodes<922:
        return '12:00:00'
    else:
        return '24:00:00'


def titan(nodes):
    if nodes<125:
        return '02:00:00'
    elif nodes<312:
        return '06:00:00'
    elif nodes<3749:
        return '12:00:00'
    else:
        return '24:00:00'


def probeFile(file, showrows=3, onlyhead=True):
    """Shows 'showrows' lines from the start, midfile and
    ending of a plaintext file
    
    Args:
        file(str): file path.
        showrows(int): rows to show from each section.
        onlyhead(bool): only print the top of the file (equivalent to 
            head -n showrows file).
    
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
    