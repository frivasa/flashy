import pandas as pd
from .IOutils import cl, np, fortParse, os
_delchar = u'!'


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
        if "setup_params" in parfile:
            dc = readSetupParams(parfile)
            self.defaults = parGroup(dc)
            self.fillcode = 1
            self.meta = getMeta(parfile)
        else:
            dc = makeParDict(parfile)
            self.params = parGroup(dc)
            self.fillcode = 0
            self.meta = {}

    def setPars(self, parfile, defaults=False):
        """sets the parameters from a file in the
        object
        """
        if defaults:
            if self.fillcode==0:
                self.defaults = readSetupParams(parfile)
                self.fillcode = 2
                self.meta = getMeta(parfile)
            else:
                self.defaults.update(readSetupParams(parfile))
                self.fillcode = 1
                self.meta = getMeta(parfile)
        else:
            if self.fillcode==1:
                self.params = parGroup(makeParDict(parfile))
                self.mergeValues()
            else:
                self.params.update(makeParDict(parfile))
                self.fillcode = 0
    
    def mergeValues(self):
        """adds parameter values to the defaults dictionary."""
        if not self.defaults or not self.params:
            print('flashy.parameterGroup: Params-Defaults pair not found. Returning.')
            return 1
        else:
            for k, v in self.params.items():
                setattr(self.defaults, k, v)
            self.fillcode = 2
    
    def tabulate(self):
        if self.fillcode==0:  # return params
            A = pd.DataFrame(list(self.params.items()), columns=['Parameter', 'Value'])
            return A.set_index('Parameter')
        elif self.fillcode==1:  # return defaults
            A = pd.DataFrame(dict(self.defaults.items()))
        else:  # return 'docked' params
            docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
            A = pd.DataFrame(dict(docked))
        A = A.transpose()
        A.index.name = 'Parameter'
        return A

    def writeParfile(self, outfile, terse=False):
        outpath = os.path.abspath(outfile)
        if self.meta:
            self.vuvuzela()
        if self.fillcode>1:
            docked = [z for z in self.defaults.items() if len(str(z[1]['value']))>0]
            writeDictionary(dict(docked), outfile, meta=True, terse=terse)
        else:
            writeDictionary(dict(self.params.items()), outfile, meta=False)

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
        geom = getattr(self.defaults, 'geometry')['value']
        if self.meta['geometry']!=geom:
            print("BZZZZZZZZZZZZ: GEOMETRY DOESN'T MATCH: "\
                  "setup:{} parfile:{}".format(self.meta['geometry'], geom))
        for k in getEssential(self.meta['dimension']):
            if k not in dkeys:
                print("BZZZZZZZZZZZZ: {} NOT SET!".format(k))

                
def getMeta(filepath):
    """Infer required properties of the run from runfolder name 
    created by flashy.setupFLASH.
    """
    path, file = os.path.split(filepath)
    runname = path.split('/')[-2]
    meta = {}
    try:
        net, geom, cells, maxblocks = runname.split('.')
    except Exception as e:
        return meta
    dimension = len(cells.split('x'))
    keys = ['network', 'geometry', 'cells', 'maxblocks', 'dimension']
    return dict(zip(keys, [net, geom, cells, maxblocks, dimension]))


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
                v = lsplit[-1].split('#')[0].strip()
                pars.append((p,v))
    return dict(pars)


def yell(k, sv, pv):
    print("#"*100)
    print("BZZZZZZZZZZZZ: {} DOESN'T MATCH:".format(k))
    print("setup:{} parfile:{}".format(sv, pv))


def getEssential(dim):
    dnames = {1:['x'], 2:['x', 'y'], 3:['x', 'y', 'z']}[dim]
    keys = []
    for dn in dnames:        
        line = 'nblock{0},{0}max,{0}min,'\
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
