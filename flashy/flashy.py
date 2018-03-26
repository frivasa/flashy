import pandas as pd
import h5py
import numpy as np
import os
import shutil
from subprocess import PIPE, Popen
import collections as cl
import butils as bu
_cdxfolder = "/cdx"
_otpfolder = "/otp"
_FLASH_DIR = "/lustre/atlas/proj-shared/csc198/frivas/00.code/FLASHOR"
_delchar = u'!'


class parameterGroup(object):
    def __init__(self, parfile):
        self.params = {}
        self.defaults = cl.OrderedDict()
        self.docked = cl.OrderedDict()
        if "setup_params" in parfile:
            self.setPars(parfile, defaults=True)
        else:
            self.setPars(parfile)

    def setPars(self, parfile, defaults=False):
        """sets the parameters from a .par file in the
        object
        """
        if defaults:
            self.defaults.update(readSetupParams(parfile))
        else:
            cd = makeParDict(parfile)
            self.params.update(cd)
            if self.defaults:
                self.addDocs()
                    
    def addDocs(self):
        if not self.defaults or not self.params:
            print "paramGroup: addDocs is missing params or defaults."
            return
        else:
            for k,v in self.params.items():
                self.docked[k] = {}
                self.docked[k]['value'] = v
                self.docked[k]['default'] = self.defaults[k]['default']
                self.docked[k]['comment'] = self.defaults[k]['comment']
                self.docked[k]['family'] = self.defaults[k]['family']

    def tabulateDocked(self, docked):
        """returns a pandas dataframe with defaults from setup_params.
        """
        if docked:
            A = pd.DataFrame(self.docked)
        else:
            A = pd.DataFrame(self.defaults)
        A = A.transpose()
        A.index.name = 'Parameter'
        return A

    def tabulate(self, defaults=False):
        """returns a pandas dataframe with parameters
        being used.
        """
        if defaults:
            return self.tabulateDocked(False)
        elif self.docked:
            return self.tabulateDocked(True)
        else:
            A = pd.DataFrame(self.params.items(), columns=['Parameter', 'Value'])
            return A.set_index('Parameter')

    def writeParfile(self, outfile, terse=False):
        outpath = os.path.abspath(outfile)
        if self.defaults and self.params:
            self.addDocs()
            writeDictionary(self.docked, outfile, meta=True, terse=terse)
        else:
            writeDictionary(self.params, outfile, meta=False)

    def comp(self):
        """compares self.defaults to self.params, returning only changed
        parameters"""
        shared = {}
        if not len(self.defaults.keys()):
            print "parameterGroup.comp: defaults not set, returning empty dict."
            return shared
        for k, v in self.params.items():
            if k in self.defaults.keys():
                if self.defaults[k]!=v:
                    shared[k] = v
        return shared

    def readChanges(self, df):
        # turn df to a simple dictionary
        if 'comment' in df.columns:
            pars, values = list(df.index), list(df['value'])
            newpdict = dict(zip(pars, values))
        else:
            newpdict = df.T.to_dict("records")[0]
        # parse values to avoid 'int'
        parsedv = [bu.fortParse(str(x)) for x in newpdict.values()]
        self.params.update(dict(zip(newpdict.keys(), parsedv)))
        # refresh docked values and delete banged(!) value parameters.
        self.addDocs()
        self.deleteSkipped()

    def deleteSkipped(self):
        for k, v in self.params.items():
            if _delchar in v:
                print "[deleteSkipped]: Removing {}".format(k)
                del self.params[k]
                del self.docked[k]


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
               geometry='cylindrical', maxbl=400, debug=False, clobber=True, portable=True):
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
    elif clobber:
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


def readSetupParams(filename):
    pardict = cl.OrderedDict()
    setp = cl.OrderedDict()
    comment = []
    for line in reversed(open(filename).readlines()):
        if line.startswith('        '):
            comment.append(line.strip())
        elif not line.startswith((' ', '\n')):
            fam = line.strip('\n ')
            for k in pardict.keys():
                pardict[k]['family'] = fam
            #setp[fam] = pardict
            setp.update(pardict)
            pardict = cl.OrderedDict()
        elif line.startswith('    '):
            par = line.strip().split()[0]
            pardict[par] = {}
            pardict[par]['value'] = ""
            pardict[par]['default'] = line.strip().split()[-1]
            pardict[par]['comment'] = " ".join(reversed(comment))
            comment = []
    return setp


def getListedDefs(supradict):
    pars, vals, defs, docs, fams = [], [], [], [], []
    for par in supradict.keys():
        pars.append(par)
        docs.append(supradict[par]['comment'])
        vals.append(supradict[par]['value'])
        defs.append(supradict[par]['default'])
        fams.append(supradict[par]['family'])
    return pars, vals, defs, docs, fams


def makeParDict(parfile):
    """ get a dictionary of parameters from a flash.par file"""
    pars = []
    with open(parfile, 'r') as par:
        for line in par:
            l = line.strip('\n')
            if l.startswith('#'):
                continue
            elif '=' in l:
                lsplit = l.split('=')
                # get name of parameter (left of '=')
                p = lsplit[0].strip()
                # remove comment from value (right of '=')
                v = lsplit[-1].split('#')[0].strip()
                pars.append((p,v))
    return dict(pars)


def writeDictionary(indict, outfile, meta=False, terse=False):
    if meta:
        ps, vs, ds, dcs, fms = getListedDefs(indict)
        dat = zip(ps,vs,ds,dcs,fms)
        dat.sort(key=lambda x: x[4])
        groups = sorted(set(fms))
        with open(outfile, 'w') as o:
            for g in groups:
                o.write("\n##### {} #####\n\n".format(g))
                vals = [x for x in dat if x[4]==g]
                maxlen = max([len(x[0]) for x in vals])
                for (p, v, d, dc, fm) in sorted(vals):
                    if terse:
                        o.write("{:{length}} = {} # {} \n".format(p, v, d, length=maxlen))
                    else:
                        o.write("{:{length}} = {} # {} {}\n".format(p, v, d, dc, length=maxlen))
    else:
        maxlen = max([len(x) for x in indict.keys()])
        with open(outfile, 'w') as o:
            for key, val in sorted(indict.items()):
                o.write("{:{pal}} = {:}\n".format(key, bu.fortParse(str(val)), pal=maxlen))


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
