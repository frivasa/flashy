from ..simulation import readTiming
from ..IOutils import pairGen
from .nucplot import plotNuclideGrid, plotReacNet
from .globals import *
_foe = 1.0e51


def plotNetwork(sim, dpi=100):
    """build a figure with enabled rates in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    f, ax = plt.subplots(figsize=(5,5), dpi=dpi)
    matsh = os.path.join(sim.netpath, 'matr_shape')
    sunet = os.path.join(sim.netpath, 'sunet')
    plotReacNet(ax, sunet, matsh, forcedZ=1e4, step=7)
    ax.tick_params(axis='both', which='both', 
                   bottom=False, right=False,
                   direction='out', length=2, width=0.5, pad=0.05)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    for side in ['bottom','right','top','left']:
        ax.spines[side].set_visible(False)
    f.tight_layout()
    return f


def plotSpeciesGrid(sim, **kwargs):
    """build a figure with nuclides included in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    f, ax = plt.subplots()
    sunet = os.path.join(sim.netpath, 'sunet')
    with open(sunet) as fl:
        species = fl.readlines()
    species = [s.strip() for s in species]
    plotNuclideGrid(ax, species, **kwargs)
    f.tight_layout()
    return f


def plotSlowCoords(sim):
    """build a figure with the evolution of the slowest cell
    through time.

    """
    xm = np.array(sim.getSlowCoord(direction='x'))
    ym = np.array(sim.getSlowCoord(direction='y'))
    zm = np.array(sim.getSlowCoord(direction='z'))
    rm = np.sqrt(xm**2+ym**2+zm**2)
    f, ax = plt.subplots()
    ax.semilogy(sim.time[:len(rm)], rm, label='r')
    ax.set_ylabel('Radius (cm)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotRayleighs(sim):
    if not sim.CJ:
        print("plot.simplot: No CJ output in object.")
        return None
    times, xins, raylins, xouts, raylouts = np.genfromtxt(sim.CJ[0], unpack=True)
    dt = np.diff(times)
    dri = np.diff(xins)
    dro = np.diff(xouts)

    f, ax = plt.subplots(figsize=(7,5))
    ind = 1  # remove zeroth checkpoint
    ax.semilogy(times[ind:], abs(dri/dt), label='Inward Shock', c='#0087ff')
    ax.semilogy(times[ind:], raylins[ind:], 
                label='CJ', c='#0087ff', lw=1, marker='o', ms=4, alpha=0.6)
    ax.semilogy(times[ind:], abs(dro/dt), label='Outward Shock', c='#d75f5f')
    ax.semilogy(times[ind:], raylouts[ind:], 
                label='CJ', c='#d75f5f',lw=1, marker='o', ms=4, alpha=0.6)
    ax.set_title('Speed')
    ax.set_ylabel('cm/s')
    ax.set_xlabel('s')
    ax.legend()
    return f


def plotTsteps(sim, burndtfactor=0.0, range=[0,0]):
    """build a figure with timestep size vs evolution time."""
    if sum(range)!=0:
        cut = slice(*range)
    else:
        cut = slice(0, None, None)
    f, ax = plt.subplots()
    ax.loglog(sim.time[cut], sim.getStepProp('dt', range=range), 
              color='k', label='picked', marker='x', ls=':')
    step = sim.steps[0]
    try:
        ax.loglog(sim.time[cut], sim.getStepProp('dt_hydro', range=range), 
                  color='b', ls=':', alpha=0.8, label='hydro')
        if not burndtfactor:
            burndtfactor = float(sim.pargroup.defaults.enucDtFactor['value'])
        ax.loglog(sim.time[cut], np.array(sim.getStepProp('dt_Burn', range=range))*burndtfactor, 
                  color='r', alpha=0.8, ls=':', label='burn*{:.2e}'.format(burndtfactor))
    except Exception as e:
        print(e)
        pass
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('Time (s)')
    ax.legend()
    f.tight_layout()
    return f


def plotBlockUse(sim, range=[0,0]):
    """build a figure with requested blocks vs evolution time.
    
    Args:
        sim(simulation): target simulation object.
        cutoff(int): plots data beyond cutoff (data[cutoff:]).
        
    Returns:
        (mpl.figure)
    
    """
    times = sim.getStepProp('t', range=range, src='refs')
    blocks = sim.getStepProp('totblocks', range=range, src='refs')
    minblk = sim.getStepProp('minblocks', range=range, src='refs')
    maxblk = sim.getStepProp('maxblocks', range=range, src='refs')
    f = plt.figure()
    layout = (3,1)
    totax = plt.subplot2grid(layout, (0, 0), rowspan=2)
    totax.loglog(times, blocks, c='k', label='total')
    totax.set_ylabel('Blocks')
    for tick in totax.get_xticklabels():
        tick.set_visible(False)
    totax.legend()
    othax = plt.subplot2grid(layout, (2, 0), sharex=totax)
    othax.loglog(times, maxblk, c='r', label='max')
    othax.loglog(times, minblk, c='b', label='min')
    othax.set_ylabel('Blocks')
    othax.set_xlabel('Time (s)')
    othax.legend()
    return f


def plotStats(sim):
    """build a figure with energy components vs sim time."""
    f, ax = plt.subplots()
    ax.loglog(sim.time, sim.getStepProp('E_total')/_foe, label='Total', color='k')
    ax.loglog(sim.time, sim.getStepProp('E_kinetic')/_foe, label='Kin', color='r', alpha=0.8, ls=':')
    ax.loglog(sim.time, sim.getStepProp('E_internal')/_foe, label='Pot', color='b', alpha=0.8, ls=':')
    ax.legend()
    ax.set_ylabel('Energy (foe)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotTiming(sim, which=0, fsize=5, ringsize=0.2, startingr=0.5):
    """Plot a nested pie chart with timings nested outwards
    
    Args:
        sim(flashy.sim): simulation object to extract data from. 
        which(int): which timing to plot.
        fsize(float): size of figure.
        ringsize(float): thickness of each ring.
        startingr(float): donut hole size.
    
    Returns:
        (mpl.fig)
    
    """
    if not sim.timings:
        print ("No timing information.")
        return None
    time, steps, data = readTiming(sim.timings[which])

    # split block and get some shape info from it
    units, deps, uts, cs, tavs, tpcs = zip(*data)
    maxdep = max(deps)+1
    delt = np.diff(deps)
    delt = np.append(0, delt)
    mm = np.zeros((len(data), maxdep), dtype="S30")
    mmval = np.zeros((len(data), maxdep))
    # build a matrix with leveled names and values
    for i in range(maxdep):
        layer = np.array(deps)-i
        pick = np.where(layer==0)[0]
        for j in range(len(units)):
            if j in pick:
                mov = sum(delt[:j+1])
                if mov>0:
                    mm[j-mov][i] = units[j].strip()
                    mmval[j-mov][i] = uts[j]
                else:
                    mm[j][i] = units[j].strip()
                    mmval[j][i] = uts[j]
            else:
                continue
    # remove emptied lines and mark roots
    pmm, pmmv = [], []
    for i, row in enumerate(mm):
        if any([bool(x) for x in row]):
            pmm.append(row)
            pmmv.append(mmval[i,:])
        else:
            continue
    pmm = np.array(pmm)
    pmmv = np.array(pmmv)
    # fill empty row values to the right and fill first non-zero to the left
    for i, row in enumerate(pmm):
        valrange = [k for k, x in enumerate(row) if bool(x)]
        for j in range(valrange[-1], maxdep):
            pmmv[i][j] = pmmv[i][valrange[-1]]
        for j in range(1, valrange[0]):
            pmmv[i][j] = pmmv[i][valrange[0]]
    rootloc = [i for i, row in enumerate(pmmv) if bool(row[0])]
    rootloc.append(len(pmmv))
    
    slices = []
    for a,b in pairGen(rootloc):
        slices.append(slice(a,b))
    names = [x.decode('utf-8') for x in pmm[:,0]]
    vals = pmmv[:,0]
    maincols = []
    for i, slc in enumerate(slices):
        grad = linear_gradient(colors[i], n=len(range(*slc.indices(1000))))
        maincols = maincols + grad['hex']

    ringnames, ringvals = [], []
    slc = slice(len(maincols))
    for i in range(1, maxdep):
        rnam = [x.decode('utf-8') for x in pmm[slc,i]] 
        if any([bool(x)  for x in rnam]):
            ringnames.append(rnam)
            ringvals.append(pmmv[slc,i])
        else:
            break
    fsize = 5
    ringsize = 0.2
    startingr = 0.5
    layout = (1, 3)
    fig = plt.figure(figsize=(2.6*fsize, fsize))
    ax = plt.subplot2grid(layout, (0, 0), colspan=1)
    lax = plt.subplot2grid(layout, (0, 2))
    lax.axis('off')
    ax.axis('equal')
    pies = []
    # always normalize so that mpl.pie plots up to the values encountered 
    # instead of filling the circle with proportional fractions.
    # https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.axes.Axes.pie.html
    pies, _ = ax.pie(vals/np.max(vals), radius=startingr, 
                     colors=maincols,
                     startangle=30,
                     wedgeprops=dict(width=ringsize, edgecolor='w', lw=0.2))
    for i in range(len(ringnames)):
        outerp, _ = ax.pie(ringvals[i]/np.max(ringvals[i]), radius=startingr+ringsize*(i+1),
                           colors=maincols, rotatelabels=True, 
                           startangle=30, counterclock=False,
                           wedgeprops=dict(width=ringsize, edgecolor='w', lw=0.2))
        pies = pies + outerp
        names = np.append(names, ringnames[i])

    a = ax.annotate("Steps: {}. Runtime: {:.5f} s".format(steps, time),
                xy=(0.0, 0.0), xytext=(0.4, 0.10), size=10,
                textcoords='figure fraction', xycoords='figure fraction', 
                bbox=dict(boxstyle='round', alpha=0.8, fc='w', ec='k'))

    # handle handles ha ha
    codes = [i for i, x in enumerate(names) if bool(x)]
    pickedhandles = [pies[i] for i in codes]
    stub = ax.legend(handles=pickedhandles, labels=[x for x in names if bool(x)], 
                     loc='center left', bbox_to_anchor=(1.2, 0.6), 
                     ncol=3, prop={'size': 10}, facecolor='w', framealpha=1.0)
    return fig