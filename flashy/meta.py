from .IOutils import _cdxfolder, _juptree, os
from .utils import chunks, scientify
from .simulation import simulation
from .datahaul.plainText import dataMatrix
import flashy.profile_workshop as pw

_errortags = [
'reached max wall clock time',
'DRIVER_ABORT: [XNet] Evolution failed to converge',
'DRIVER_ABORT: [Eos] Error: too many Newton-Raphson iterations in eos_helmholtz',
'DRIVER_ABORT: [eos_helm] ERROR: abar is negative.'
]

def checkCodeFolder(folder, names=['r_match_outer', 't_ignite_outer'], probesim=False):
    """extracts metadata from all runs found within a code folder."""
    if folder[-1]!='/':
        filename = '{}.csv'.format(os.path.basename(folder))
    else:
        filename = '{}.csv'.format(os.path.basename(os.path.dirname(folder)))
    allvalues = []
    for i, f in enumerate(sorted(os.listdir(folder))):
        if f==_cdxfolder:
            continue
        runpath = os.path.join(folder, f)
        print('Checking:', os.path.basename(runpath))
        tgs, values = getRunMeta(runpath, names=names, probesim=probesim)
        print('\t', values[-2])  # this is the end condition
        # add relative path, this should break if run from elsewhere than 10.yt_FLASH
        values.append('=HYPERLINK("{}{}")'.format(_juptree, runpath[3:]))  # remove '../'
        allvalues += values
    tgs.append('url')
    with open(filename, 'w') as otp:
        otp.write(','.join(tgs))
        otp.write('\n')
        for vals in chunks(allvalues, len(values)):
            otp.write(','.join([str(x) for x in scientify(vals)]))
            otp.write('\n')
    print('Wrote:', filename)


def getRunMeta(runpath, names=['t_ignite_outer', 'r_match_outer'], probesim=True):
    """return named parameters plus simulation status."""
    tags = names.copy()
    values = []
    sim = simulation(runpath)
    if probesim:
        nodes, info = sim.pargroup.probeSimulation(verbose=False)
        print("\n".join(info[:-3]))
    # get properties of profile
    pprof = ['Tmass', 'HeMass', 'Rho_c', 'Rho_he', 'WDR', 'matchHeight']
    pobj = dataMatrix(sim.profile)
    intfpos = pw.getInterfacePosition(pobj)
    masses = pw.getSummedMasses(pobj)
    interface = pobj.radius[intfpos]
    values = [sum(list(masses.values())), masses['he4'], 
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
    tags.append('Comment')
    maxt = float(sim.pargroup.defaults.tmax['value'])
    if not sim.steps:  # no steps read, so no .log
        simt = 0.0
        endc = 'No .log file'
        comm = 'Simulation has not been run.'
    else:
        lastst = sim.steps[-1]
        endt = lastst.dt+lastst.t
        # check if finished
        if (maxt-endt)>1e-6:  # log file precision
            # print('Simtime/Max simtime: {:.2f}/{:.2f} s. Error below:'.format(endt, maxt))
            errtag, errlog = showerror(sim.otpfiles[-1])
            if '[Eos]' in errtag:
                comm = '{:e} {:e} {:e}'.format(*eosLocation(errlog))
            else:
                comm = ''
            simt = endt
            endc = errtag
        else:
            simt = endt
            endc = 'Completed sim.'
            comm = ''
    values.append(simt)
    values.append(endc)
    values.append(comm)
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
    