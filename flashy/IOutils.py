import os
import sys
import shutil
import h5py
import collections as cl
from flashy.utils import np
from subprocess import PIPE, Popen
from PIL import Image
import imageio


def turn2cartesian(folder, prefix='all', nowitness=False):
    if prefix=='all':
        finns = getFileList(folder)
        finns += getFileList(folder, prefix='chk')
    else:
        finns = getFileList(folder)
    finns = [f for f in finns if "cart_" not in f]
    for finn in finns:
        jake = os.path.join(folder,'cart_'+finn)
        if os.path.exists(jake):
            print "{} found. Skipping.".format(jake)
            continue
        switchGeometry(os.path.join(folder,finn), jake, verbose=True)
        if nowitness:
            os.remove(os.path.join(folder,finn))


def switchGeometry(file, output, verbose=True):
    finn = h5py.File(file, "r")
    jake = h5py.File(output, "w")
    for k in finn.iterkeys():
        finn.copy(k, jake)
    ds = jake[u'string scalars']
    newt = np.copy(ds[...])
    newt[0][0] = ds[0][0].replace("cylindrical", "cartesian  ")
    ds[...] = newt
    ds2 = jake[u'string runtime parameters']
    newt2 = np.copy(ds2[...])
    for i, v in enumerate(ds2):
        if "cylindrical" in v[0]:
            newt2[i][0] = v[0].replace("cylindrical", "cartesian  ")
    ds2[...] = newt2
    finn.close()
    jake.close()
    if verbose:
        print("Wrote {} from {}".format(output, file))


def setupFLASH(module, runfolder, kwargs={'threadBlockList':'true'}, nxb=4, nyb=8, nzb=0,
               geometry='cylindrical', maxbl=500, debug=False, portable=True):
    """calls ./setup at _FLASH_DIR with given parameters,
    writing the code to destination for compiling afterwards.
    (must be run on a Py2.X kernel)

    Arguments:
        module(str): name of Simulation folder to setup.
        runfolder(str): run folder (creates _cdx, and _otp within it).
        kwargs(dict): keyword arguments to setup
        n[xyz]b(int): cells per block for setup.
        geometry(str): cartesian, spherical, cylindrical(default).
        maxbl(int): maximum blocks per proc. elem.
        debug(bool): show terminal output.

    """
    destination = runfolder + _cdxfolder
    path = os.path.abspath(destination)
    if not os.path.exists(destination):
        os.makedirs(destination)
    else:
        print 'Emptying {}'.format(destination)
        shutil.rmtree(destination)
        os.makedirs(destination)
    try:
        os.makedirs(runfolder + _otpfolder)
    except:
        pass
    if not nzb:
        if not nyb:
            dimstr = "-1d -nxb={}".format(nxb)
        else:
            dimstr = "-2d -nxb={} -nyb={}".format(nxb, nyb)
    else:
        dimstr = "-3d -nxb={} -nyb={} -nzb={}".format(nxb, nyb, nzb)
    portag = ''
    if portable:
        portag = '-portable'
    cstub = 'cd {} && ./setup {} -auto {} '\
            '-objdir="{}" {} -geometry={} -maxblocks={} '
    comm = cstub.format(_FLASH_DIR, module, portag, path,
                        dimstr, geometry, maxbl)
    kwstr = ''
    for (k, v) in kwargs.items():
        if '+' in k:
            kwstr+=' {}'.format(k)
        else:
            kwstr+=' {}={}'.format(k,v)
    #kwstr = ' '.join(['{}={}'.format(k,v) for (k, v) in kwargs.items()])
    comm = comm + kwstr
    print comm
    p = Popen(['/bin/bash'], stdin=PIPE, stdout=PIPE)
    out, err = p.communicate(input=comm.encode())
    exitcode = p.returncode
    if debug:
        print out
        print err
    return comm, exitcode


def compileFLASH(runfolder, resultlines=20, slines=[], procs=8):
    """calls 'make -j procs' at outpath, compiling the run
    """
    comm = 'make -j {}'.format(procs)
    if slines:
        comm = '{} && {}'.format(" && ".join(slines), comm)
    path = os.path.abspath(runfolder)
    p = Popen(['/bin/bash'], cwd=path, stdin=PIPE, stdout=PIPE)
    out, err = p.communicate(input=comm.encode())
    print err
    print out
    exitcode = p.returncode
    return path, exitcode


def getFileList(folder, prefix='plt'):
    return sorted([x for x in os.listdir(folder) if prefix in x])


def cpFLASHrun(runfolder, newrunfol):
    """copy the cdx folder to a new runfolder"""
    src = os.path.abspath(runfolder) + _cdxfolder
    dst = os.path.abspath(newrunfol) + _cdxfolder
    if os.path.exists(dst):
        shutil.rmtree(dst)
        shutil.copytree(src, dst)
    else:
        shutil.copytree(src, dst)
    os.makedirs(newrunfol + _otpfolder)


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
    print "\n\tSaved: {}".format(expname)


def fortParse(arg):
    """returns a parsed variable from a parameter (bool,
    str, or number)

    Args:
        arg (str): parameter value

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
            return '"{}"'.format(arg.strip('"\' '))


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
    for f in files:
        shutil.copy('/'.join([src,f]), '/'.join([dst,f]))


def writeSubmit(subfile, code, pbsins=[],
                time='12:00:00', nodes=1252, ompth=16, proj='', mail=''):
    """builds a submit.pbs with a typical header, specifying walltime and nodes,
    then adding slines of code below. Exports OMP_NUM_THREADS=ompth
    titan: aprun (-j1) -n 1 -d 16
    debug: -D (int)
    rhea: mpirun --map-by ppr:N:node:pe=Th or -n
    debug: --display-map / --report-bindings
    Rhea max: 48 hours on 16 nodes (2x8 core p/node: -np 256)
    Titan: <125 nodes 2h, <312 nodes 6h...
    """
    subHeader = []
    subHeader.append('#!/bin/bash')
    subHeader.append('#PBS -V') # pass env vars to nodes
    #subHeader.append('#PBS -j oe') # join err and otp
    subHeader.append('#PBS -l gres=atlas1')
    subHeader.append('#PBS -A {}'.format(proj))
    subHeader.append('#PBS -l walltime={},nodes={}'.format(time, nodes))
    subHeader.append('#PBS -N {}'.format(os.path.basename(subfile)[:-4]))
    subHeader.append('#PBS -j oe')
    subScript = []
    subScript.append('date')
    subScript.append('echo Submitted from: $PBS_O_WORKDIR')
    subScript.append('echo #####################')
    subScript.append('export OMP_NUM_THREADS={}'.format(int(ompth)))
    subScript.append('export CRAY_CUDA_MPS=1')
    if mail:
        subHeader.append('#PBS -M {}'.format(mail))
        subHeader.append('#PBS -m abe')
    if pbsins:
        subHeader = subHeader + pbsins
    subScript = subScript + code
    with open(subfile, 'w') as o:
        o.write("\n".join(subHeader))
        o.write("\n")
        o.write("\n".join(subScript))
        o.write("\n")


def probeFile(file, showrows=3):
    """Shows 'showrows' lines from the start, midfile and
    ending of a plaintext file
    """
    with open(file, 'r') as f:
        lines = f.readlines()
    l = len(lines)
    if l < 10:
        print("".join(lines))
    else:
        l2 = int(l/2.)
        print("".join(lines[0:showrows]))
        print("".join(lines[l2:l2+showrows]))
        print("".join(lines[-showrows:]))


def emptyFileTree(root):
    """Empties 'root' folder."""
    path = os.path.abspath(root)
    shutil.rmtree(path)
    os.makedirs(path)
    