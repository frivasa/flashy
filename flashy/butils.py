import os
import shutil
import numpy as np
from subprocess import PIPE, Popen
from PIL import Image
import imageio
import sys


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
    subScript.append('wait')
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
