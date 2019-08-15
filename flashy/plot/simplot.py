from ..simulation import readTiming
from ..IOutils import pairGen
from .nucplot import plotNuclideGrid, plotReacNet
from .globals import *
from ..paramSetup import pd
from astropy.time import Time
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
    """plot cj vs rayleigh for calculated 1D (Lineout) speeds."""
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
    """build a figure with timestep size vs evolution time.

    Args:
        sim(simulation): target simulation object.
        burdtfactor(float): scale burn dt manually.
        range(int list): index cutoffs for data.

    Returns:
        (mpl.figure)

    """
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


def plotIRLsteps(sim, sigma=15.0, checkpoints=False):
    """plot walltime steps vs size of step.

    Args:
        sim(simulation): target simulation object.

    Returns:
        (mpl.figure)

    """
    timestamps = sim.getStepProp('timestamp')
    deltas = sim.getStepProp('irldelta')
    nums = sim.getStepProp('n')
    astTtype = [Time(t, format='datetime') for t in timestamps]
    astroTs =  [a.mjd for a in astTtype]
    tsteps = [d.seconds for d in deltas]
    f, ax = plt.subplots()
    ax.scatter(astroTs, tsteps, color='peru', marker='.')
    # overplot step numbers for outliers
    locs = np.where(np.array(tsteps) > sigma*sim.irlstep.seconds)[0]
    for l in locs:
        ax.text(astroTs[l], tsteps[l]+1.0,  # offset
                "{}".format(int(nums[l])), fontsize=10,
                color='royalblue')
    # overplot checkpoints
    if checkpoints:
        timestamps = sim.getStepProp('timestamp', src='chkp')
        types = sim.getStepProp('filetype', src='chkp')
        numbers = sim.getStepProp('number', src='chkp')
        astTtype = [Time(t, format='datetime') for t in timestamps]
        astroTs =  [a.mjd for a in astTtype]
        for p, t, n in zip(astroTs, types, numbers):
            if t == 'checkpoint':
                c = 'r'
                l = 'checkpoints'
                fac = 2.6
            else:
                c = 'g'
                l = 'plotfiles'
                fac = 2.0
            ax.axvline(p, c=c, alpha=0.2, label=l)
            ax.text(p, fac*sim.irlstep.seconds,
                    "{}".format(n), fontsize=10,
                    alpha=0.2, color=c)
        handles, labels = ax.get_legend_handles_labels()
        handle_list, label_list = [], []
        for handle, label in zip(handles, labels):
            if label not in label_list:
                handle_list.append(handle)
                label_list.append(label)
        ax.legend(handle_list, label_list)
#     ax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('MJD')
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
    if np.max(blocks)/np.min(blocks) < 10.0:
        totax.semilogx(times, blocks, c='k', label='total')
    else:
        totax.loglog(times, blocks, c='k', label='total')
    totax.set_ylabel('Blocks')
    stub = [tick.set_visible(False) for tick in totax.get_xticklabels()]
    stub = [tick.set_visible(False)
            for tick in totax.get_xticklabels(minor=True)]
    totax.legend()
    othax = plt.subplot2grid(layout, (2, 0), sharex=totax)
    othax.semilogx(times, maxblk, c='r', label='max')
    othax.semilogx(times, minblk, c='b', label='min')
    othax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    othax.set_ylabel('Blocks')
    othax.set_xlabel('Time (s)')
    # othax.legend()
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


def plotSunburstTiming(sim, which=0, avgProc=True, column=0):
    """Plot a nested pie chart with timings nested outwards
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
    time, steps, colnames, data = readTiming(sim.timings[which],
                                             avgProc=avgProc)
    # split block and get some shape info from it
    units, deps, fulldata = zip(*data)
    units = [u.strip() for u in units]
    colname = colnames[column+1]
    allvals = [v[column] for v in fulldata]
    delt = np.diff(deps)
    delt = np.append(0, delt)
    # find the families
    roots = [i for i, d in enumerate(deps) if d == 0]
    rootblockGen = pairGen(roots+[len(deps)])
    slices = []
    for a, b in rootblockGen:
        slices.append(slice(a, b))
    # sort the data
    labels, parents, ids = [], [], []
    values, metav = [], []
    for i, slc in enumerate(slices):
        names = units[slc]
        itvals = allvals[slc]
        sdeps = deps[slc]
        delts = delt[slc]
        # traverse delts extracting parent trees
        currp = [""]
        for j, n in enumerate(delts):
            if n > 0:
                currp.append(names[j-n].strip())
            elif n < 0:
                for k in range(abs(n)):
                    del currp[-1]
            par = '-'.join(currp).strip('-')
            id = '-'.join(currp[1:]+[names[j]]).strip('-')
            parents.append(par)
            ids.append(id)
        maxval = np.max(itvals)
        normalization = np.array(itvals)/maxval*100
        normalization[0] = maxval/time*100
        values = values + itvals
        metav = metav + list(normalization)
        labels = labels + names
    df = pd.DataFrame(data={'labels': labels, 'parents': parents,
                            'vals': values, 'meta': metav,
                            'ids': ids})
    trace = gob.Sunburst(ids=df.ids, labels=df.labels, values=df.vals,
                         parents=df.parents, meta=df.meta,
                         maxdepth=5, branchvalues='remainder',
                         hovertemplate="<b>%{label}</b><br>" +
                                       "%{value:.0f}<br>" +
                                       "%{meta:.4f} <extra></extra>")
    tsize = 16  # making space for title
    forstr = "Var:{}<br>Total:{:.2f} Steps:{:d}"
    tag = forstr.format(colname, time, steps)
    layout = gob.Layout(margin=gob.layout.Margin(t=tsize*5, l=0, r=0, b=0),
                        sunburstcolorway=["#636efa", "#EF553B", "#00cc96",
                                          "#ab63fa", "#19d3f3", "#e763fa",
                                          "#FECB52", "#FFA15A", "#FF6692"],
                        extendsunburstcolors=True,
                        title=tag)
    fig = gob.Figure(data=[trace], layout=layout)
    return fig


def plotPizzaTiming(sim, which=0, avgProc=True, column=0, clrshift=1):
    """Plot a nested pie chart with timings nested outwards
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
        clrshift(int): change color list start (change colorscheme).

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
    return fig
