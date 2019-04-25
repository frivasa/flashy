import pandas as pd
from .IOutils import cl, np, fortParse, os, getWalltime, getMachineSubJob, writeSchedulerScript
from .IOutils import _cdxfolder, _otpfolder, _logfile, _statsfile
from .utils import cart2sph
_FLASHdefaults = 'setup_params'
_strmax = 49
pd.set_option('display.max_colwidth', 0)

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
        # fillcode is a workaround to avoid creating empty parGroups.
        # 0: pars only, 1: +defaults, 2: both and merged
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
    
    def purgeGroup(self):
        """clears default value fields and removes self.params"""
        if self.meta['code']==0:
            print('Erasing params.')
            self.params = parGroup({})
        elif self.meta['code']==1:
            print('Clearing defaults.')
            for k, v in self.defaults.items():
                setattr(self.defaults, k, '')
        else:
            print('Clearing parameterGroup')
            self.params = parGroup({})
            for k, v in self.params.items():
                setattr(self.defaults, k, '')
    
    def mergeValues(self):
        """adds parameter values to the defaults dictionary."""
        if not self.defaults or not self.params:
            print('flashy.parameterGroup: Params-Defaults pair not found. Returning.')
            return 1
        else:
            for k, v in self.params.items():
                setattr(self.defaults, k, fortParse(v, dec=False))
            self.meta['code'] = 2
    
    def tabulate(self, allpars=False):
        if not self.meta['code']:  # return params
            A = pd.DataFrame(list(self.params.items()), columns=['Parameter', 'Value'])
            return A.set_index('Parameter')
        elif self.meta['code']==1 or allpars:  # return defaults
            A = pd.DataFrame(dict(self.defaults.items()))
        else:  # return 'docked' params
            try:
                docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
            except TypeError:
                print('parGroup.tabulate.error: string index error in docked.'\
                      ' flash.par includes unknown params to simulation.')
                return pd.DataFrame()
            A = pd.DataFrame(dict(docked))
        A = A.transpose()
        A.index.name = 'Parameter'
        if 'comment' in A.columns:
            return A[['value', 'default', 'family', 'comment']]
        else:
            return A

    def getStyledTable(self, allpars=False, stylerprops={}, tableprops={}):
        #if not stylerprops:
        #    stylerprops = {'background-color':'#111111', 'color': '#dbe1ea'}
        #if not tableprops:
        #    tableprops = {'selector':"tr:hover", 'props':[("background-color", '#6f5757')]}
        if self.meta['code']==0:
            print('No defaults found, returning params only')
            A = self.tabulate()
            return A.style
        elif self.meta['code']==1:
            print('No set values, returning defaults')
            A = self.tabulate()
            return A.style
        else:
            A = self.tabulate(allpars=allpars)
            S = A.style
            #S.set_properties(**stylerprops)
            #S.set_table_styles([tableprops])
            redind = A.index[A['value']!=A['default']].tolist()
            S.applymap(stylerTest, subset=pd.IndexSlice[redind, ['value']])
            return S

    def writeParfile(self, outfile='', terse=False):
        if outfile:
            outpath = os.path.abspath(outfile)
            print(outpath)
            try:
                writeDictionary(self.getDocked(), outpath, meta=True, terse=terse)
            except Exception as e:
                print('Failed to write documentation, writing params only.')
                writeDictionary(dict(self.params.items()), outpath, meta=False)
            print('Wrote: {}'.format(outpath))
        else:  # default is to assume 'docked' params and write to runfolder/otp
            self.vuvuzela()
            cdx = self.meta['cdxpath']
            cpname = os.path.join(cdx, 'flash.par')
            writeDictionary(self.getDocked(), os.path.abspath(cpname), meta=True, terse=terse)
            print('Wrote: {}'.format(cpname))
            # create checkpoint folder
            otpf = self.defaults.output_directory['value']
            outpath = os.path.join(cdx, otpf)
            os.makedirs(outpath, exist_ok=True)
            # write flash.par to runfolder
            otpf = "../{}/".format(self.meta['runname'])
            outpath = os.path.join(cdx, otpf)
            cpname = os.path.join(outpath, 'flash.par')
            writeDictionary(self.getDocked(), os.path.abspath(cpname), meta=True, terse=terse)
            print('Wrote: {}'.format(cpname))

    def vuvuzela(self):
        """Sound the horn of ERROR."""
        try:
            dkeys = [z[0] for z in self.defaults.items() if len(str(z[1]['value']))>0]
        except TypeError:
            raise Exception('Parsing Error: this is due to set parameters '\
                  'which do not exist in params_setup (e.g.: Xnet vs ap13)')
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
        # print(mins, maxs, nblocks)
        return [float(n) for n in nblocks], [float(m) for m in mins], [float(m) for m in maxs]
    
    def probeSimulation(self, frac=0.4, forcePEs=0, verbose=True):
        dim, cells, maxbl = self.readMeta()
        nblocks, mins, maxs = self.readEssential()
        area, tblcks, tcells = 1.0, 1.0, 1.0
        rmax = float(self.defaults.lrefine_max['value'])
        sotp = []
        for i in range(dim):
            dname = {0:'x', 1:'y', 2:'z'}[i]
            span = maxs[i]-mins[i]
            limblcks = np.power(2, (rmax-1))*float(nblocks[i])
            limcells = limblcks*cells[i]
            #minspan = span/limcells
            sotp.append('{} span: {:2.4E}, resolution: {:2.4E}, cpb: {}'.format(dname, span, span/limcells, cells[i]))
            #area*=span
            tblcks*=limblcks
            tcells*=limcells
        # mult(spans)/mult(nblocks)/mult(cells)/2^(ref-1)/2^(ref-1) = area of cell
        # ref 1 is nblocks, therefore ref-1
        # cells per block might change if blocks are rectangular (not cubes)
        sotp.append('Max Blocks per PE: {:>12.0f}'.format(maxbl))
        sotp.append('Max Refinement:{:>16.0f}'.format(rmax))
        sotp.append('Max blocks/cells:{:>14.4E}/{:.4E}'.format(tblcks, tcells))
        #sotp.append('Resolution: {:E}'.format(np.sqrt(area/tcells)))  # this is not correct
        
        if forcePEs:
            sotp.append('forced PEs: {:0.0f}'.format(forcePEs))
            nodes = forcePEs
        else:
            maxPEs = tblcks/maxbl
            sotp.append('Maximum PEs: {:>18.0f}'.format(maxPEs))
            sotp.append('Optimistic alloc ({:>4.0%}): {:>6.2f}'.format(frac, maxPEs*frac))
            nodes = int(maxPEs*frac)+1
        if verbose:
            print('\n'.join(sotp))
        return nodes, sotp
    
    def writeRunFiles(self, frac=0.4, forcePEs=0, terse=True, ddt=False,
                      multisub=True, prefix='', suffix='', IOwindow=120, forceWallT='',
                      proj='csc198', machine='titan', **kwargs):
        """Probes the parameters, sets up required resources, and writes 
        necessary files based on a stringent structure.
        
        Args:
            frac(float): reduce allocation by frac.
            terse(bool): add descriptions to parameters in the par file.
            multisub(bool): activate iterator (see flashy.IOutils).
            ddt(bool): enable arm-forge debugging connection.
            prefix(str): pass a prefix to runname generator.
            suffix(str): pass a suffix to runname generator.
            IOwindow(int): seconds to extract from walltime to write Checkpoints.
            proj(str): allocation project code.
            machine(str): machine being used for batch submit.
            forcePEs(int): force a number of nodes.
            forceWallT(str): force a walltime (hh:mm:ss).
        
        """
        # sets otp_directory and runname key in meta
        self.generateRunName(suffix=suffix, prefix=prefix, checkpath=self.meta['cdxpath'])
        # estimate allocation
        nodes, _ = self.probeSimulation(frac=frac,  forcePEs=forcePEs)
        if forceWallT:
            time = forceWallT
        else:
            time = getWalltime(nodes, machine=machine)
        inputs = np.array([float(x) for x in time.split(':')])
        factrs = np.array([3600.0, 60.0, 1.0])
        seconds = sum(inputs*factrs) - IOwindow  # time for last checkpoint
        self.defaults.wall_clock_time_limit = int(seconds)
        # write parfiles in both cdx and otp folders 
        self.writeParfile(terse=terse)
        # write submit script at otp folder
        auxpath = '../{}/'.format(self.meta['runname'])
        outpath = os.path.join(self.meta['cdxpath'], auxpath)
        subpath = os.path.join(outpath, '{}'.format(self.meta['runname']))
        scheduler, ext = self.writeSubmit(subpath, proj, machine=machine, nodes=nodes, 
                                          time=time, multisub=multisub, ddt=ddt, **kwargs)
        return '{} {}{}'.format(scheduler, subpath, ext)
    
    def writeSubmit(self, submitpath, proj, machine='titan', time='02:00:00', nodes=16, ompth=16, 
                    multisub=True, ddt=False, **kwargs):
        qsubfold, qsubname = os.path.split(submitpath)
        runf = os.path.abspath(self.meta['cdxpath'])
        otpf = self.defaults.output_directory['value']
        auxf = '../{}'.format(self.meta['runname'])
        code = []
        # move where the action is and get the corresponding flash.par
        code.append('cd {}'.format(runf))
        code.append('cp {} .'.format(os.path.join(auxf, 'flash.par')))
        # not sure if this still goes for alpine...
        if self.meta['dimension'] > 1:
            if nodes>512:  # hear the warnings, set limit to 512
                code.append('lfs setstripe -c 512 {}'.format(os.path.join(otpf)))
            else:
                code.append('lfs setstripe -c {} {}'.format(nodes, os.path.join(otpf)))
        launcher, scheduler, ext, schedulercode = getMachineSubJob(machine, proj, time, nodes, ompth, 
                                                                   ddt, os.path.join(runf, auxf), **kwargs)
        # print(os.path.join(runf, auxf))
        if ddt:
            code.append('module load forge/18.3')  # specify version to latest
            launcher = 'ddt --connect {}'.format(launcher)
        otpname = submitpath + ext
        code.append(launcher)
        if multisub:
            nendestimate = self.defaults.tmax['value']/self.defaults.checkpointFileIntervalTime['value']
            if nendestimate<1.0:
                nendestimate = 10
            # export chaining arguments and apply iterator to set the file number
            code.insert(len(code)-1, 'export QSUBFOLD={}'.format(os.path.abspath(qsubfold)))
            code.insert(len(code)-1, 'export QSUBNAME={}{}'.format(qsubname, ext))
            code.insert(len(code)-1, 'bash iterator {} flash.par {}'.format(otpf, int(nendestimate)))
            code.insert(len(code)-1, 'wait')
            code.append('wait')
            code.append('cd $QSUBFOLD')
            code.append('{} $QSUBNAME'.format(scheduler))
        code.append('wait')
        writeSchedulerScript(otpname, code, schedulercode)
        print('Wrote: {}'.format(otpname))
        return scheduler, ext
        
    def generateRunName(self, prefix='', suffix='', checkpath=""):
        basename = getProfilePrefix(self.defaults.initialWDFile['value'])
        match = (self.defaults.x_match['value'], 
                 self.defaults.y_match['value'], 
                 self.defaults.z_match['value'])
        direc = cart2sph(*match)
        rout, rin = self.defaults.r_match_outer['value'], self.defaults.r_match_inner['value']
        size = rout-rin
        dim = self.meta['dimension']
        if sum(direc)==0.0:
            basename += '_r{}'.format(int(self.defaults.r_match_inner['value']/1e5))
        elif dim==3:
            basename += '_p{}_{}_{}'.format(int(direc[0]/1e5), int(direc[1]), int(direc[2]))
        elif dim==2:  # 2D is x-y, not x-z as in canonical spherical, hence direc[2]
            basename += '_p{}_{}'.format(int(direc[0]/1e5), int(direc[2]))
        else:
            basename += '_p{}'.format(int(direc[0]/1e5))
        temp = self.defaults.t_ignite_outer['value']
        # some params might be str from setPars >> float(temp)
        runname = '{}_ms{}_t{:.1f}'.format(basename, int(size/1e5), float(temp)/1e9)
        if prefix:
            runname = "_".join([prefix, runname])
        if suffix:
            runname = "_".join([runname, suffix])    
        print ('Run name generated: {}'.format(runname))
        self.defaults.run_comment = runname
        # check length of runname for FLASH output str limit and add a number to the run
        if len(runname)> _strmax:
            print('paramSetup.generateRunName.WARNING: '\
                  'Run name too long for FLASH. clipping to {} chars.'.format(_strmax))
            runname = runname[:_strmax]
        num = len(os.listdir(os.path.split(checkpath)[0]))-1
        runname = "{:02}{}".format(num, runname)
        print ('Run name used: {}'.format(runname))
        
        self.defaults.geometry = self.meta['geometry']
        # separate auxiliary files from checkpoints to stripe otp folder.
        self.meta['runname'] = runname
        self.defaults.output_directory = "../{}/{}/".format(runname, _otpfolder)
        self.defaults.log_file = "../{}/{}".format(runname, _logfile)
        self.defaults.stats_file = "../{}/{}".format(runname, _statsfile)
        
    def getDocked(self, onlyvals=False):
        """returns default values that have been changed through 
        reading pars or manually changed
        """
        group = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
        if onlyvals:
            return dict([(a, b['value']) for (a, b) in group])
        else:
            return dict(group)


def comPars(par1, par2):
    """compare two flash .par files
    
    Args:
        par1(str): reference .par filename.
        par1(str): comparison .par filename.
        
    Returns:
        (pandas.Styler): highlighted pandas dataframe.
    
    """
    d1 = makeParDict(par1)
    d2 = makeParDict(par2)
    return comDicts(d1, d2)


def comDicts(d1, d2):
    """compare two parameter dictionaries.
    
    Args:
        d1(str): reference dict.
        d2(str): comparison dict.
        
    Returns:
        (pandas.Styler): highlighted pandas dataframe.
    
    """
    A = pd.DataFrame(list(d1.items()), columns=['Parameter', 'Value'])
    B = pd.DataFrame(list(d2.items()), columns=['Parameter', 'Value'])
    joined = pd.merge(A, B, on='Parameter')
    redind = joined.index[joined['Value_x']!=joined['Value_y']].tolist()
    S = joined.style
    S.applymap(stylerTest, subset=pd.IndexSlice[redind, ['Value_x']])
    return S


def getProfilePrefix(string):
    """returns a prefix for a runfolder from a given profile filename
    profile format: wdSource_shellSource_mass_shellmass.dat
    
    """
    profilename = os.path.basename(string[:-4])
    try:
        wd, shell, wdm, shm = profilename.split('_')
    except ValueError:
        return 'custom'
    code = wd[0]+shell[0]
    mass = getFloat(wdm)
    shma = getFloat(shm)
    return '{}_m{:.0f}_sh{:.0f}'.format(code, mass*1e3, shma*1e3)


def getFloat(string):
    """returns rightmost float in a string."""
    for i in range(len(string)):
        try:
            val = float(string[i:])
            return val
        except ValueError:
            continue
    return 0.0


def stylerTest(value):
    """stub to change colors in selected cells."""
    #return 'color: #ec971f'  # orange
    return 'color: #FF0000'  # red


def getMeta(filepath):
    """Infer required properties of the run from runfolder name 
    created by the method flashy.setupFLASH.
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
                setp.update(pardict)
                pardict = {}
            elif line.startswith('    '):
                par = line.strip().split()[0]
                pardict[par] = {}
                pardict[par]['value'] = ""
                defa = line.strip().split('[')[-1].strip(' ]')
                pardict[par]['default'] = fortParse(defa)
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
    """split a defaults dictionary into arrays."""
    pars, vals, defs, docs, fams = [], [], [], [], []
    for par in supradict.keys():
        pars.append(par)
        docs.append(supradict[par]['comment'])
        vals.append(supradict[par]['value'])
        defs.append(supradict[par]['default'])
        fams.append(supradict[par]['family'])
    return pars, vals, defs, docs, fams


def writeDictionary(indict, outfile, meta=False, terse=False):
    """write a formatted parameter list to a file."""
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
                    if '"' in str(fortParse(str(v))):  # yeah this is ugly
                        strlen = len(p)
                    else:
                        strlen = maxlen
                    if terse:
                        o.write("{:{length}} = {} # {} \n".format(p, fortParse(str(v)), d, length=strlen))
                    else:
                        o.write("{:{length}} = {} # {} {}\n".format(p, fortParse(str(v)), d, dc, length=strlen))
    else:
        maxlen = max([len(x) for x in indict.keys()])
        with open(outfile, 'w') as o:
            for key, val in sorted(indict.items()):
                if '"' in str(fortParse(str(val))):  # yeah this is ugly
                    strlen = len(key)
                else:
                    strlen = maxlen
                o.write("{:{pal}} = {:}\n".format(key, fortParse(str(val)), pal=strlen))

