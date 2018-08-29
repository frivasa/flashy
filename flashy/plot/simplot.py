# from ..simulation import simulation
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


def plotSpeciesGrid(sim):
    """build a figure with nuclides included in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    f, ax = plt.subplots()
    sunet = os.path.join(sim.netpath, 'sunet')
    with open(sunet) as fl:
        species = fl.readlines()
    species = [s.strip() for s in species]
    plotNuclideGrid(ax, species)
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
#     ax.semilogy(sim.time[1:], xm, label='x')
#     ax.semilogy(sim.time[1:], ym, label='y')
    ax.semilogy(sim.time[:len(rm)], rm, label='r')
#     ax.semilogy(sim.time[1:], rm, color='k')
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


def plotTsteps(sim, burndtfactor=0.0):
    """build a figure with timestep size vs evolution time."""
    f, ax = plt.subplots()
    ax.loglog(sim.time, sim.getStepProp('dt'), color='k', label='simulation')
    try:
        ax.loglog(sim.time, sim.getStepProp('dthydro'), 
                  color='b', ls=':', alpha=0.8, label='hydro')
        if not burndtfactor:
            burndtfactor = float(sim.pargroup.defaults.enucDtFactor['default'])
            print(burndtfactor)
        ax.loglog(sim.time, np.array(sim.getStepProp('dtburn'))/burndtfactor, 
                  color='r', alpha=0.8, ls=':', label='burn/{:.2e}'.format(burndtfactor))
    except Exception as e:
        print(e)
        pass
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('Time (s)')
    ax.legend()
    f.tight_layout()
    return f


def plotBlockUse(sim):
    """build a figure with requested blocks vs evolution time."""
    f, ax = plt.subplots()
    ax.semilogx(sim.time, sim.getStepProp('blocks'), color='k')
    ax.set_ylabel('Blocks Requested')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotStats(sim):
    """build a figure with energy components through evolution time."""
    f, ax = plt.subplots()
    ax.loglog(sim.time, sim.getStepProp('E_total')/_foe, label='Total', color='k')
    ax.loglog(sim.time, sim.getStepProp('E_kinetic')/_foe, label='Kin', color='r', alpha=0.8, ls=':')
    ax.loglog(sim.time, sim.getStepProp('E_internal')/_foe, label='Pot', color='b', alpha=0.8, ls=':')
    ax.legend()
    ax.set_ylabel('Energy (foe)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotTiming(sim, depth=2, cut=0.01):
    if not sim.timings:
        print ("No timing information.")
        return None
    steps, time, delta, md = sim.mergeTimings()
    fig = plt.figure(figsize=(9, 9))
    layout = (1, 3)
    # Species
    ax1 = plt.subplot2grid(layout, (0, 0), colspan=2)
    a = ax1.annotate("Steps: {}\nRuntime: {:.5f} s\nAvg.Step:{:.5e} s".format(steps, time, delta),
                    xy=(0.0, 0.0), xytext=(0.6, 0.10), size=12,
                    textcoords='figure fraction', xycoords='figure fraction', 
                    bbox=dict(boxstyle='round', fc='w', ec='k'))
    for i, k in enumerate(md.keys()):
        probe = probeAtDepth(md, k, depth)
        for (lab, val) in probe:
            if val < cut:
                continue
            ax1.bar(i, val, label=lab)
    ax1.set_ylabel('Percentage(%)')
    ax1.set_xticks(range(len(md.keys())))
    ax1.set_xticklabels(md.keys(), rotation=-30, ha='left', rotation_mode='anchor')
    ax1.tick_params(pad=6)
    ax1.set_yscale('log')
    lgd = ax1.legend(ncol=1, loc='upper left', bbox_to_anchor=(1.00, 1.05), 
                     columnspacing=0.0,
                     numpoints=3)
    plt.tight_layout()
    return fig


def probeAtDepth(tree, stem, depth, value='percent'):
    props = ['secs', 'calls', 'percent']
    vals = []
    for k in tree[stem]:
        if k in props:
            continue
        if tree[stem][k]['depth']==depth:
            vals.append((k, tree[stem][k][value]))
    return vals

