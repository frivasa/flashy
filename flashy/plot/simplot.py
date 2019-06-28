from ..simulation import readTiming
from ..IOutils import pairGen
from .nucplot import plotNuclideGrid, plotReacNet
from .globals import *
import plotly.graph_objs as gob
from plotly.offline import (download_plotlyjs, init_notebook_mode, iplot)
_foe = 1.0e51


def plotNetwork(sim, dpi=100):
    """build a figure with enabled rates in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    f, ax = plt.subplots(figsize=(5, 5), dpi=dpi)
    matsh = os.path.join(sim.netpath, 'matr_shape')
    sunet = os.path.join(sim.netpath, 'sunet')
    plotReacNet(ax, sunet, matsh, forcedZ=1e4, step=7)
    ax.tick_params(axis='both', which='both',
                   bottom=False, right=False,
                   direction='out', length=2, width=0.5, pad=0.05)
    ax.xaxis.set_label_position('top')
    ax.xaxis.tick_top()
    for side in ['bottom', 'right', 'top', 'left']:
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
    times, xins, raylins, xouts, raylouts = np.genfromtxt(sim.CJ[0],
                                                          unpack=True)
    dt = np.diff(times)
    dri = np.diff(xins)
    dro = np.diff(xouts)

    f, ax = plt.subplots(figsize=(7, 5))
    ind = 1  # remove zeroth checkpoint
    ax.semilogy(times[ind:], abs(dri/dt), label='Inward Shock', c='#0087ff')
    ax.semilogy(times[ind:], raylins[ind:],
                label='CJ', c='#0087ff', lw=1, marker='o', ms=4, alpha=0.6)
    ax.semilogy(times[ind:], abs(dro/dt), label='Outward Shock', c='#d75f5f')
    ax.semilogy(times[ind:], raylouts[ind:],
                label='CJ', c='#d75f5f', lw=1, marker='o', ms=4, alpha=0.6)
    ax.set_title('Speed')
    ax.set_ylabel('cm/s')
    ax.set_xlabel('s')
    ax.legend()
    return f


def plotTsteps(sim, burndtfactor=0.0, range=[0, 0]):
    """build a figure with timestep size vs evolution time."""
    if sum(range) != 0:
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
        burnarray = np.array(sim.getStepProp('dt_Burn', range=range))
        ax.loglog(sim.time[cut], burnarray*burndtfactor, color='r',
                  alpha=0.8, ls=':', label='burn*{:.2e}'.format(burndtfactor))
    except Exception as e:
        print(e)
        pass
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('Time (s)')
    ax.legend()
    f.tight_layout()
    return f


def plotBlockUse(sim, range=[0, 0]):
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
    layout = (3, 1)
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
    ax.loglog(sim.time, sim.getStepProp('E_total')/_foe,
              label='Total', color='k')
    ax.loglog(sim.time, sim.getStepProp('E_kinetic')/_foe,
              label='Kin', color='r', alpha=0.8, ls=':')
    ax.loglog(sim.time, sim.getStepProp('E_internal')/_foe,
              label='Pot', color='b', alpha=0.8, ls=':')
    ax.legend()
    ax.set_ylabel('Energy (foe)')
    ax.set_xlabel('Time (s)')
    f.tight_layout()
    return f


def plotTiming(sim, which=0, avgProc=True, column=0, clrshift=1):
    """
    Plot a nested pie chart with timings nested outwards
    # setup plotly from notebook:
    import plotly
    plotly.tools.set_credentials_file(username='NAME', api_key='APKEY')
    avgProc False cols:
        time sec  num calls   secs avg  time pct
    avgProc True cols:
        max/proc (s)  min/proc (s) avg/proc (s)   num calls

    Args:
        sim(flashy.sim): simulation object to extract data from.
        which(int): which timing to plot.
        avgProc(bool): averaged procs (True) or only proc 0 (False).
        column(int): column from timing to plot.

    Returns:
        (mpl.fig)

    """
    if not sim.timings:
        print("No timing information.")
        return None
    try:
        init_notebook_mode(connected=True)
    except:
        print('Plotly not configured, see docstring.')
        return None
    time, steps, colnames, data = readTiming(sim.timings[which],
                                             avgProc=avgProc)
    # split block and get some shape info from it
    units, deps, fulldata = zip(*data)
    colname = colnames[column+1]
    uts = [v[column] for v in fulldata]
    maxdep = max(deps)+1
    delt = np.diff(deps)
    delt = np.append(0, delt)

    roots = [i for i, d in enumerate(deps) if d == 0]
    rootblockGen = pairGen(roots+[len(deps)])

    slices = []
    for a, b in rootblockGen:
        slices.append(slice(a, b))

    traces, mcolors = [], []
    shift = 0
    for i, slc in enumerate(slices):
        names = units[slc]
        values = uts[slc]
        sdeps = deps[slc]
        delts = delt[slc]
        drange = range(1, np.max(sdeps)+1)
        if delts[0] < 0:
            delts[0] = 0
        meta = []
        parents = [names[0].strip()]
        for j, n in enumerate(delts):
            if n > 0:
                parents.append(names[j-n].strip())
            elif n < 0:
                for k in range(abs(n)):
                    del parents[-1]
            meta.append(parents[-1].strip())
        mcolors.append(colors[i+clrshift])
        if not len(drange):
            shift += 1
            continue
        for m, d in enumerate(drange):
            subvals = [values[k] for k, j in enumerate(sdeps) if j == d]
            subnames = [names[k].strip()
                        for k, j in enumerate(sdeps) if j == d]
            submeta = [meta[k] for k, j in enumerate(sdeps) if j == d]
            grad = linear_gradient(colors[i+clrshift], n=len(meta))
            trace = gob.Pie(values=subvals, labels=subnames,
                            sort=True, direction='clockwise', opacity=0.9,
                            text=submeta, hole=None,
                            hovertemplate="<b>%{label}</b><br>" +
                                          "%{text}<br>" +
                                          "%{percent}<br><extra></extra>",
                            title=dict(text=parents[0],
                                       font=dict(color='black', size=20),
                                       position='bottom center'),
                            domain={'row': i-shift, 'column': m},
                            textinfo='value', rotation=20,
                            textposition='inside',
                            marker={'colors': grad['hex'],
                                    'line': {'color': 'black',
                                             'width': 1.0}})
            traces.append(trace)

    # finally overplot depth 0
    hnames = [units[i] for i, d in enumerate(deps) if d == 0]
    hvals = [uts[i] for i, d in enumerate(deps) if d == 0]

    htrace = gob.Pie(values=hvals, labels=hnames,
                     sort=True, direction='clockwise', opacity=1.0,
                     hole=0.2, hoverinfo='label+value',
                     title=dict(text='Root Calls',
                                font=dict(color='black', size=20),
                                position='top center'),
                     domain=dict(row=0, column=maxdep-2),
                     textinfo='percent', rotation=20,
                     textposition='inside',
                     marker={'colors': mcolors,
                             'line': {'color': 'black', 'width': 1.0}})
    traces.append(htrace)
    cstr = "Variable: {}<br>".format(colname)
    tstr = " Total Time:{:.4f}<br>".format(time)
    sstr = "Steps:{:d}".format(steps)
    layout = gob.Layout(showlegend=False,
                        grid=dict(rows=len(slices)-shift, columns=5),
                        title=dict(text=cstr + tstr + sstr,
                                   font=dict(color='black', size=20)),
                        autosize=False,
                        height=250*(i+1-shift), width=200*(maxdep-1))
    fig = gob.FigureWidget(data=traces, layout=layout)
    iplot(fig, filename='test')
