from .IOutils import _cdxfolder, _juptree, os
from .utils import chunks, scientify, pd
from .simulation import simulation
from .datahaul.plainText import dataMatrix
import flashy.profile_workshop as pw
_errortags = [
    'reached max wall clock time',
    'DRIVER_ABORT: [XNet] Evolution failed to converge',
    'DRIVER_ABORT: [Eos] Error: '
    'too many Newton-Raphson iterations in eos_helmholtz',
    'DRIVER_ABORT: [eos_helm] ERROR: abar is negative.'
]


def checkCodeFolder(folder, names=['r_match_outer', 't_ignite_outer'],
                    probesim=False, succinct=True):
    """extracts metadata from all runs found within a code folder."""
    allvalues = []
    for i, f in enumerate(sorted(os.listdir(folder))):
        if f == _cdxfolder:
            continue
        runpath = os.path.join(folder, f)
        print('\t', os.path.basename(runpath), '\n')
        tgs, values = getRunMeta(runpath, names=names,
                                 probesim=probesim, succinct=succinct)
        print('\t', values[-1], '\n')  # this is the end condition
        values.append('{}{}'.format(_juptree, runpath[3:]))
        values.append(os.path.basename(runpath))
        allvalues += values
    tgs.append('url')
    tgs.append('name')
    table = pd.DataFrame()
    for vals in chunks(allvalues, len(values)):
        subT = pd.DataFrame([vals], columns=tgs)
        table = table.append(subT)
    return table.set_index('name')


def getRunMeta(runpath, names=['t_ignite_outer', 'r_match_outer'],
               probesim=True, succinct=True):
    """return named parameters plus simulation status."""
    tags = names.copy()
    values = []
    sim = simulation(runpath)
    if probesim:
        nodes, info = sim.pargroup.probeSimulation(verbose=False)
        print("\n".join(info[:-3]))
        print(sim.quickLook(refsimtime=0.05, refstep=80))
    # get properties of profile
    if not succinct:
        pprof = ['TotM', 'CoreM', 'EnvM', 'Rho_c', 'Rho_he', 'WDR', 'matchHeight']
        pobj = dataMatrix(sim.profile)
        intfpos = pw.getInterfacePosition(pobj)
        interface = pobj.radius[intfpos]
        corems = pw.getSummedMasses(pobj, range=(None, intfpos))
        shellm = pw.getSummedMasses(pobj, range=(intfpos, None))
        values = [pobj.masses[-1],
                  sum(list(corems.values())),
                  sum(list(shellm.values())),
                  pobj.density[0], pobj.density[intfpos], interface]
        # get the match position and height
        x = getattr(sim.pargroup.defaults, 'x_match')['value']
        y = getattr(sim.pargroup.defaults, 'y_match')['value']
        z = getattr(sim.pargroup.defaults, 'z_match')['value']
        r_match = (x*x+y*y+z*z)**0.5
        values.append(r_match - interface)
        tags = pprof+tags
    # fill with required parameters first since flash.par is always there
    for n in names:
        fl = getattr(sim.pargroup.defaults, n)['value']
        values.append(fl)
    # add additional metadata
    tags.append('simtime')
    tags.append('Ending Condition')
    maxt = float(sim.pargroup.defaults.tmax['value'])
    if not sim.steps:  # no steps read, so no .log
        simt = 0.0
        endc = 'No .log file (sim has not run)'
    else:
        lastst = sim.steps[-1]
        endt = lastst.dt+lastst.t
        # check if finished
        if (maxt-endt) > 1e-6:  # log file precision
            errtag, errlog = showerror(sim.otpfiles[-1])
            if '[Eos]' in errtag:
                comm = 'EoS: {:e} {:e} {:e}'.format(*eosLocation(errlog))
            else:
                comm = errtag
            simt = endt
            endc = comm
        else:
            simt = endt
            endc = 'Completed sim.'
    values.append(simt)
    values.append(endc)
    return tags, scientify(values)


def showerror(filename, searchlines=20):
    """artisanal tail to avoid using subprocess"""
    with open(filename, 'r') as f:
        data = f.read()
    lines = data.split('\n')
    # reduce to searchlines
    errsearch = lines[-searchlines:]
    textblock = '\n'.join(errsearch)
    for r in _errortags:
        if r in textblock:
            return r, textblock
    return 'Error not tabulated.', textblock


def eosLocation(errblock):
    """return location for the specific [Eos] type error."""
    for l in errblock.split('\n'):
        if '(x, y, z)' in l:
            x, y, z = [float(x) for x in l.split()[-3:]]
    return x, y, z
