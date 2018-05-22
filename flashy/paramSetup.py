import pandas as pd
from .IOutils import cl, np, fortParse, os, _cdxfolder, getTITANtime, writePBSscript
_FLASHdefaults = 'setup_params'


class parGroup(object):
    def __setattr__(self, att, val):
        if hasattr(self, att):
            comp = getattr(self, att)
            if isinstance(comp, dict):  # here were at a defaults obj (values are dicts)
                if isinstance(val, dict):  # init/overwriting a value
                    super().__setattr__(att, val)
                else:  # setting the value field
                    comp['value'] = val
                    super().__setattr__(att, comp)
            else:  # this is a param obj (only keys and values)
                super().__setattr__(att, val)
        else:  # new att, initialize
            super().__setattr__(att, val)
    
    def __init__(self, dictionary):
        for k, v in dictionary.items():
            setattr(self, k, v)
    
    def update(self, dictionary):
        for k, v in dictionary.items():
            setattr(self, k, v)
            
    def items(self):
        keys = self.__dict__
        values = [getattr(self, att) for att in keys]
        return zip(keys, values)


class parameterGroup(object):
    def __init__(self, parfile):
        """fillcode is a workaround to avoid creating empty parGroups."""
        self.meta = getMeta(parfile)
        if not self.meta['code']:
            dc = makeParDict(parfile)
            self.params = parGroup(dc)
        elif self.meta['code']==1:
            dc = readSetupParams(parfile)
            self.defaults = parGroup(dc)
        else:
            dc = makeParDict(parfile)
            self.params = parGroup(dc)
            dc = readSetupParams(os.path.join(self.meta['cdxpath'], _FLASHdefaults))
            self.defaults = parGroup(dc)
            self.mergeValues()

    def setPars(self, parfile):
        """sets the parameters from a file in the object."""
        newmeta = getMeta(parfile)
        if newmeta['default']:  # add or update defaults
            if not self.meta['code']:
                self.defaults = readSetupParams(parfile)
                self.meta = newmeta
            else:
                self.defaults.update(readSetupParams(parfile))
        else:  # overwrite params (safer)
            self.params = parGroup(makeParDict(parfile))
            if self.meta['code']:
                self.mergeValues()
    
    def mergeValues(self):
        """adds parameter values to the defaults dictionary."""
        if not self.defaults or not self.params:
            print('flashy.parameterGroup: Params-Defaults pair not found. Returning.')
            return 1
        else:
            for k, v in self.params.items():
                setattr(self.defaults, k, v)
            self.meta['code'] = 2
    
    def tabulate(self):
        if not self.meta['code']:  # return params
            A = pd.DataFrame(list(self.params.items()), columns=['Parameter', 'Value'])
            return A.set_index('Parameter')
        elif self.meta['code']==1:  # return defaults
            A = pd.DataFrame(dict(self.defaults.items()))
        else:  # return 'docked' params
            docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
            A = pd.DataFrame(dict(docked))
        A = A.transpose()
        A.index.name = 'Parameter'
        return A

    def writeParfile(self, outfile='', terse=False):
        if outfile:
            outpath = os.path.abspath(outfile)
            print(outpath)
            try:
                docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
                writeDictionary(dict(docked), outpath, meta=True, terse=terse)
            except Exception as e:
                print('Failed to write documentation, writing params only.')
                writeDictionary(dict(self.params.items()), outpath, meta=False)
            print('Wrote: {}'.format(outpath))
        else:  # default is to assume 'docked' params and write to runfolder/otp
            self.vuvuzela()
            docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
            cdx = self.meta['cdxpath']
            cpname = os.path.join(cdx, 'flash.par')
            writeDictionary(dict(docked), os.path.abspath(cpname), meta=True, terse=terse)
            print('Wrote: {}'.format(cpname))
            otpf = self.defaults.output_directory['value']
            outpath = os.path.join(cdx, otpf)
            if not os.path.exists(outpath):
                os.makedirs(outpath)
            cpname = os.path.join(outpath, 'flash.par')
            writeDictionary(dict(docked), os.path.abspath(cpname), meta=True, terse=terse)
            print('Wrote: {}'.format(cpname))

    def readChanges(self, df):
        # turn df to a simple dictionary
        if 'comment' in df.columns:
            pars, values = list(df.index), list(df['value'])
            newpdict = dict(zip(pars, values))
        else:
            newpdict = df.T.to_dict("records")[0]
        # parse values to avoid 'int'
        parsedv = [fortParse(str(x), dec=False) for x in newpdict.values()]
        self.params.update(dict(zip(newpdict.keys(), parsedv)))
        # refresh docked values and remove empty value parameters for when retabulating.
        self.mergeValues()
        
    def vuvuzela(self):
        """Sound the horn of ERROR."""
        dkeys = [z[0] for z in self.defaults.items() if len(str(z[1]['value']))>0]
        geom = self.defaults.geometry
        if self.meta['geometry']!=geom['value']:
            print("BZZZZZZZZZZZZ: GEOMETRY DOESN'T MATCH: "\
                  "setup:{} parfile:{}".format(self.meta['geometry'], geom['value']))
        for k in getEssential(self.meta['dimension']):
            if k not in dkeys:
                print("BZZZZZZZZZZZZ: {} NOT SET!".format(k))
    
    def readMeta(self):
        """returns dimension, cells per block, and maxblocks from a 'docked' parfile"""
        dim = int(self.meta['dimension'])
        cells = self.meta['cells'].split('x')
        cells[-1] = cells[-1].replace('cells', '')
        cells = list(map(int, cells))
        maxblocks = float(self.meta['maxblocks'].replace('maxb',''))
        return dim, cells, maxblocks

    def readEssential(self):
        """returns nblocks, minima, and maxima from a 'docked' parfile"""
        dim = int(self.meta['dimension'])
        keys = getEssential(dim)
        step = int(len(keys)/dim)
        nblocks = [getattr(self.defaults, k)['value'] for k in keys[0::step]]
        mins = [getattr(self.defaults, k)['value'] for k in keys[1::step]]
        maxs = [getattr(self.defaults, k)['value'] for k in keys[2::step]]
        return [float(n) for n in nblocks], [float(m) for m in mins], [float(m) for m in maxs]
    
    def probeSimulation(self, frac=0.6):
        dim, cells, maxbl = self.readMeta()
        nblocks, mins, maxs = self.readEssential()
        tblcks, tcells = 1.0, 1.0
        rmax = float(self.defaults.lrefine_max['value'])
        sotp = []
        for i in range(dim):
            dname = {0:'x', 1:'y', 2:'z'}[i]
            span = maxs[i]-mins[i]
            limblcks = np.power(2, (rmax-1))*float(nblocks[i])
            limcells = limblcks*cells[i]
            minspan = span/limcells
            sotp.append('{} span: {:E}'.format(dname, span))
            tblcks*=limblcks
            tcells*=limcells
            # print(limblcks, limcells)
        maxPEs = tblcks/maxbl
        sotp.append('Max Refinement: {:0.0f}'.format(rmax))
        sotp.append('Resolution: {:E}'.format(minspan))
        sotp.append('Maximum cells: {:E}'.format(tcells))
        sotp.append('Maximum Blocks: {:E}'.format(tblcks))
        sotp.append('Max Blocks per PE: {:0.0f}'.format(maxbl))
        sotp.append('Maximum PEs: {:0.0f}'.format(maxPEs))
        sotp.append('Optimistic alloc (60%): {:0.0f}'.format(maxPEs*frac))
        print('\n'.join(sotp))
        return int(maxPEs*frac+1)
    
    def writeSubmit(self, recommended=True, frac=0.6, j1=False,
                    time='02:00:00', nodes=16, ompth=16):
        qsubfold, qsubname = os.path.split(submitpath)
        runf = os.path.abspath(self.meta['cdxpath'])
        code = []
        code.append('export QSUBFOLD={}'.format(os.path.abspath(qsubfold)))
        code.append('export QSUBNAME={}'.format(qsubname))
        code.append('cd {}'.format(runf))
        code.append('bash iterator {} flash.par'.format(self.defaults.output_directory['value']))
        code.append('wait')
        if recommended:
            nodes = self.probeSimulation()
            time = getTITANtime(nodes)
        if j1:
            nodes*=2
            ompth = 8
            code.append('aprun -n{} -d{} -j1 ./flash4 &'.format(nodes, ompth))
        else:
            code.append('aprun -n{} -d{} ./flash4 &'.format(nodes, ompth))
        code.append('wait')
        code.append('cd $QSUBFOLD')
        code.append('qsub $QSUBNAME')
        writePBSscript(submitpath, code, time=time, nodes=nodes, ompth=ompth,
                       proj='csc198', mail='rivas.aguilera@gmail.com',abe='a')
                
def getMeta(filepath):
    """Infer required properties of the run from runfolder name 
    created by flashy.setupFLASH.
    """
    path, file = os.path.split(filepath)
    # path is either ../runcode/cdx or ../runcode/otp_***
    runpath, fldr = os.path.split(path)
    _, runname = os.path.split(runpath)
    meta = {}
    meta['cdxpath'] = os.path.join(runpath, _cdxfolder)
    deffile = os.path.join(meta['cdxpath'], _FLASHdefaults)
    if file==_FLASHdefaults:
        meta['default'] = 1
        meta['code'] = 1  # read only defaults
    elif os.path.isfile(deffile):
        meta['default'] = 0
        meta['code'] = 2  # read both params and defaults
    else:
        meta['default'] = 0
        meta['code'] = 0  # read only params
    try:
        net, geom, cells, maxblocks = runname.split('.')
    except Exception as e:
        return meta
    dimension = len(cells.split('x'))
    keys = ['network', 'geometry', 'cells', 'maxblocks', 'dimension']
    newkeys = dict(zip(keys, [net, geom, cells, maxblocks, dimension]))
    meta.update(newkeys)
    return meta


def readSetupParams(filename):
    pardict = {}
    setp = {}
    comment = []
    try:
        for line in reversed(open(filename).readlines()):
            if line.startswith('        '):
                comment.append(line.strip())
            elif not line.startswith((' ', '\n')):
                fam = line.strip('\n ')
                for k in pardict.keys():
                    pardict[k]['family'] = fam
                #setp[fam] = pardict
                setp.update(pardict)
                pardict = {}
            elif line.startswith('    '):
                par = line.strip().split()[0]
                pardict[par] = {}
                pardict[par]['value'] = ""
                pardict[par]['default'] = line.strip().split()[-1]
                pardict[par]['comment'] = " ".join(reversed(comment))
                comment = []
    except FileNotFoundError:
        print('paramSetup.readSetupParams: Defaults file not found, returning empty dict.')
        pass  # return an empty dictionary
    return setp


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
                v = lsplit[1].split('#')[0].strip()
                pars.append((p,v))
    return dict(pars)


def getEssential(dim):
    """Returns parsed names for essential parameters in a simulation.
    Bleeding edge of inference here, careful with changing order of 
    variables...
    
    """
    dnames = {1:['x'], 2:['x', 'y'], 3:['x', 'y', 'z']}[dim]
    keys = []
    for dn in dnames:        
        line = 'nblock{0},{0}min,{0}max,'\
               '{0}l_boundary_type,{0}r_boundary_type'.format(dn)
        keys += line.split(',')
    return keys


def getListedDefs(supradict):
    pars, vals, defs, docs, fams = [], [], [], [], []
    for par in supradict.keys():
        pars.append(par)
        docs.append(supradict[par]['comment'])
        vals.append(supradict[par]['value'])
        defs.append(supradict[par]['default'])
        fams.append(supradict[par]['family'])
    return pars, vals, defs, docs, fams


def writeDictionary(indict, outfile, meta=False, terse=False):
    if meta:
        ps, vs, ds, dcs, fms = getListedDefs(indict)
        dat = list(zip(ps,vs,ds,dcs,fms))
        sorted(dat, key=lambda x: x[4])
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
                o.write("{:{pal}} = {:}\n".format(key, fortParse(str(val)), pal=maxlen))


# def 
