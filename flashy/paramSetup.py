import pandas as pd
from IOutils import cl, np, fortParse, os
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
        parsedv = [fortParse(str(x)) for x in newpdict.values()]
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
                o.write("{:{pal}} = {:}\n".format(key, fortParse(str(val)), pal=maxlen))
