import os
import sys
import shutil
import h5py
import collections as cl
from flashy.utils import np
from subprocess import PIPE, Popen
from PIL import Image
import imageio
_cdxfolder = "cdx"
_otpfolder = "otp"
_FLASH_DIR = "/lustre/atlas/proj-shared/csc198/frivas/00.code/FLASHOR"
_AUX_DIR = "/lustre/atlas/proj-shared/csc198/frivas/"


def turn2cartesian(folder, prefix='all', nowitness=False):
    """Iterates over files within a folder, switching the geometry of 
    hdf5 files found to cartesian.
    
    Args:
        folder(str): folder path.
        prefix(str): filter string (defaults to all files in the folder).
        nowitness(bool): remove non-modified files.
    
    """
    
    if prefix=='all':
        finns = getFileList(folder)
        finns += getFileList(folder, prefix='chk')
    else:
        finns = getFileList(folder)
    finns = [f for f in finns if "cart_" not in f]
    for finn in finns:
        jake = os.path.join(folder,'cart_'+finn)
        if os.path.exists(jake):
            print("{} found. Skipping.".format(jake))
            continue
        switchGeometry(os.path.join(folder,finn), jake, verbose=True)
        if nowitness:
            os.remove(os.path.join(folder,finn))


def switchGeometry(file, output, verbose=True):
    """copies hdf5 file, changing the coordinate system name to 
    cartesian for yt input.
    
    Args:
        file(str): input filename.
        output(str): output filename.
        verbose(bool): Report file creation.
    
    """
    finn = h5py.File(file, "r")
    jake = h5py.File(output, "w")
    # p2 > p3: obj.iterkeys() > obj.keys()
    # p2 > p3: hdf5 works with bytes, not str: u"" > b""
    for k in finn.keys():
        finn.copy(k, jake)
    ds = jake[u'string scalars']
    newt = np.copy(ds[...])
    newt[0][0] = ds[0][0].replace(b"cylindrical", b"cartesian  ")
    ds[...] = newt
    ds2 = jake[u'string runtime parameters']
    newt2 = np.copy(ds2[...])
    for i, v in enumerate(ds2):
        if b"cylindrical" in v[0]:
            newt2[i][0] = v[0].replace(b"cylindrical", b"cartesian  ")
    ds2[...] = newt2
    finn.close()
    jake.close()
    if verbose:
        print("Wrote {} from {}".format(output, file))


def setupFLASH(module, runfolder='', kwargs={'threadBlockList':'true'}, nbs=[16, 16, 16],
               geometry='cylindrical', maxbl=500, debug=False):
    """calls ./setup at _FLASH_DIR with given parameters,
    writing the code to runfolder. (FLASH setup script runs only on py 2.X).

    Arguments:
        module(str): name of Simulation folder to setup.
        runfolder(str): run folder (creates _cdx, and _otp within it).
        kwargs(dict): keyword arguments to setup
        nbs(int tuple): cells per block for setup per dimension.
        geometry(str): cartesian, spherical, cylindrical(default).
        maxbl(int): maximum blocks per proc. elem.
        debug(bool): show terminal output.

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
    try:
        os.makedirs(os.path.join(runfolder, _otpfolder))
    except:
        pass
    
    if len(nbs)==1:
        dimstr = "-1d -nxb={}".format(*nbs)
    elif len(nbs)==2:
        dimstr = "-2d -nxb={} -nyb={}".format(*nbs)
    else:
        dimstr = "-3d -nxb={} -nyb={} -nzb={}".format(*nbs)
    
    cstub = 'python2 {}/bin/setup.py {} -auto '\
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
    print(comm)
    p = Popen(['/bin/bash'], stdin=PIPE, stdout=PIPE)
    out, err = p.communicate(input=comm.encode())
    exitcode = p.returncode
    if debug:
        print(out.decode())
        print(err)
    print('generated run name {}'.format(name))
    
    # copy parameter varying script
    try:
        iterpath = os.path.join(_AUX_DIR, '07.miscellaneous/bash/iterator')
        shutil.copy2(iterpath, destination)
        print('copied iterator')
    except:
        print('bash iterator not found, skipping.')

    return comm, exitcode


def getFileList(folder, prefix='plt', fullpath=False):
    """Returns a filename list subject to a prefix 'glob'.
    
    Args:
        folder(str): look-in path.
        prefix(str): filter string for files
        fullpath(bool): return absolute path for each file.
        
    Returns:
        (str list)
    
    
    """
    names = sorted(os.listdir(folder))
    fnames = [x for x in names if prefix in x]
    if fullpath:
        return [os.path.join(os.path.abspath(folder), x) for x in fnames]
    else:
        return fnames


def cpFLASHrun(runfolder, newrunfol):
    """copy the cdx folder to a new runfolder"""
    src = os.path.join(os.path.abspath(runfolder), _cdxfolder)
    dst = os.path.join(os.path.abspath(newrunfol), _cdxfolder)
    if os.path.exists(dst):
        shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copytree(src, dst)
    os.makedirs(os.path.join(newrunfol, _otpfolder))


def makeGIF(runfolder, prefix='', subf='', speed=0.2):
    """Join all png images within a folder in an animated .gif

    Args:
        runfolder (str): folder path
        prefix (str): prefix for files
        subf (str): subfolder for files
        speed (float): seconds between frames
    
    """
    if not subf:
        prepath = runfolder
    else:
        prepath = os.path.join(runfolder, "png", subf)
    finns = sorted([i for i in os.listdir(prepath) if prefix in i])
    finns = [x for x in finns if '.png' in x]
    if not prefix:
        prefix = subf
    # image resize
    #    for finn in finns:
    #        name, ext = os.path.splitext(finn)
    #        im = Image.open(finn)
    #        im.thumbnail(size)
    #        im.save(name + ".resized.png", "png")
    jakes = []
    for finn in finns:
        sys.stdout.write(finn + " ")
        jakes.append(imageio.imread(os.path.join(prepath,finn)))
    if prefix:
        expname = "{}/{}.gif".format(runfolder, prefix)
    else:
        expname = "{}/joined.gif".format(runfolder, prefix)
    imageio.mimsave(expname, jakes, format='gif', duration=speed)
    print("\n\tSaved: {}".format(expname))


def fortParse(arg, dec=True):
    """returns a parsed variable from a parameter (bool,
    str, or number)

    Args:
        arg(str): parameter value.
        dec(bool): add "" to strings for printing.

    Returns:
        str: decorated argument for fortran parsing.

    """
    try:
        val = np.float(arg.replace('d','E'))
        return arg
    except ValueError:
        if '.true.' in arg.lower():
            return arg.strip()
        elif '.false.' in arg.lower():
            return arg.strip()
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


def writePBSscript(subfile, code, pbsins=[],
                   time='12:00:00', nodes=1252, ompth=16, proj='', mail='', abe='abe'):
    """PBS submit system file cooker.
    builds a submit.pbs with a typical header, specifying walltime and nodes,
    then adding slines of code below. Exports OMP_NUM_THREADS=ompth
    titan: aprun (-j1) -n 1 -d 16
    debug: -D (int)
    rhea: mpirun --map-by ppr:N:node:pe=Th or -n
    debug: --display-map / --report-bindings
    Rhea max: 48 hours on 16 nodes (2x8 core p/node: -np 256)
    Titan: <125 nodes 2h, <312 nodes 6h...
    
    Args:
        subfile(str): submit filename (also set as jobname)
        code(str list): commands to insert in the file
        pbsins(str list): extra PBS directives.
        time(str): walltime request.
        nodes(int): nodes to request.
        ompth(int): omp thread number.
        proj(str): project code.
        mail(str): notification e-mail.
    
    """
    subHeader = []
    subHeader.append('#!/bin/bash')
    subHeader.append('#PBS -V')  # pass env vars to nodes
    subHeader.append('#PBS -l gres=atlas1')
    subHeader.append('#PBS -A {}'.format(proj))
    subHeader.append('#PBS -l walltime={},nodes={}'.format(time, nodes))
    subHeader.append('#PBS -N {}'.format(os.path.basename(subfile)[:-4]))
    subHeader.append('#PBS -j oe')  # join err and otp
    subScript = []
    subScript.append('echo Submitted from: $PBS_O_WORKDIR')
    subScript.append('date')
    subScript.append('export OMP_NUM_THREADS={}'.format(int(ompth)))
    subScript.append('export CRAY_CUDA_MPS=1')
    if mail:
        subHeader.append('#PBS -M {}'.format(mail))
        subHeader.append('#PBS -m {}'.format(abe))
    if pbsins:
        subHeader = subHeader + pbsins
    subScript = subScript + code
    with open(subfile, 'w') as o:
        o.write("\n".join(subHeader))
        o.write("\n")
        o.write("\n".join(subScript))
        o.write("\n")


def getTITANtime(nodes):
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
    