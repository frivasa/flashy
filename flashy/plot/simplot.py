# from ..simulation import simulation
from .nucplot import plotNuclideGrid, plotReacNet
from .oneDim import plainTprofile
from .globals import *
_foe = 1.0e51

def plotInitialProfile(sim):
    return plainTprofile(sim.profile)


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


def plotTsteps(sim):
    """build a figure with timestep size vs evolution time."""
    f, ax = plt.subplots()
    ax.semilogx(sim.time, sim.getTsteps(), color='k')
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotBlockUse(sim):
    """build a figure with requested blocks vs evolution time."""
    f, ax = plt.subplots()
    ax.semilogx(sim.time, sim.getBlocks(), color='k')
    ax.set_ylabel('Blocks Requested')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotStats(sim):
    """build a figure with energy components through evolution time."""
    f, ax = plt.subplots()
    ax.loglog(sim.time, sim.E_total/_foe, label='Total')
    ax.loglog(sim.time, sim.E_kinetic/_foe, label='Kin')
    ax.loglog(sim.time, sim.E_internal/_foe, label='Pot')
    ax.legend()
    ax.set_ylabel('Energy (foe/bethe)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f
