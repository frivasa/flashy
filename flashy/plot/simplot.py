from ..simulation import readTiming
from ..IOutils import pairGen
from ..paramSetup import pd
from .nucplot import plotNuclideGrid, plotReacNet
from .globals import (linear_gradient, log,
                      np, os, plt, colors, rollingAverage,
                      StrMethodFormatter)
from astropy.time import Time
#import plotly.graph_objs as gob
#from plotly.offline import (download_plotlyjs, init_notebook_mode, iplot)
_foe = 1.0e51


def plotNetwork(sim, netpath='', dpi=100, step=4, aspect=1.0, cmap='Blues'):
    """build a figure with enabled rates in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    elif netpath:
        path = netpath
    else:
        path = sim.netpath
    f, ax = plt.subplots(figsize=(9, 9), dpi=dpi)
    matsh = os.path.join(path, 'matr_shape')
    sunet = os.path.join(path, 'sunet')
    plotReacNet(ax, sunet, matsh, step=step, aspect=aspect)
    # for side in ['bottom', 'right', 'top', 'left']:
    #     ax.spines[side].set_visible(False)
    return f


def plotSpeciesGrid(sim, netpath='', **kwargs):
    """build a figure with nuclides included in a network."""
    if not sim.netpath:
        print("plot.simplot: Run doesn't use XNet.")
        return None
    elif netpath:
        path = netpath
    else:
        path = sim.netpath
    f, ax = plt.subplots()
    sunet = os.path.join(path, 'sunet')
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
    astroTs = [a.mjd for a in astTtype]
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
        astroTs = [a.mjd for a in astTtype]
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
    ax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    ax.set_ylabel('Step size (s)')
    ax.set_xlabel('MJD')
    return f


def plotImprint(sim, pePerNode=42, rollAvgStep=80, dpi=80,
                blkprank=512, bot=4e-3, showcheckpoints=False,
                irltime=False):
    """plot all block data alongside ranks, timesteps and submit number
    for all the simtime achieved by the simulation.

    Args:
        sim(flashy.simulation): simulation object to probe.
        pePerNode(int): ranks per node on machine ("-r x -a").
        rollAvgStep(int): width of cell averaging for IRL timesteps.
        dpi(int): dpi of figure returned.
        blkprank(int): max blocks per rank (sets ylim of lower plot).
        bot(float): top plot bottom limit.
        showcheckpoints(bool): plot checkpoints/plotfiles lines.
        irltime(bool): use timestamps for x-axis.

    Returns:
        (mpl.figure)

    """
    legdict = {'ncol': 1, 'loc': 'upper left', 'columnspacing': 0.0,
               'labelspacing': 0.1, 'numpoints': 2, 'handletextpad': 0.2,
               'bbox_to_anchor': (1.0, 1.02)}
    totblk = sim.getStepProp('totblocks')
    totblk = totblk.astype(float)
    minblk = sim.getStepProp('minblocks')
    maxblk = sim.getStepProp('maxblocks')
    # nanoseconds is the default for np.timedelta
    tsteps = [s.item()*1e-9 for s in sim.getStepProp('irldelta')]
    subs = sim.getStepProp('submit_number')
    PE = sim.getStepProp('mpi_ranks')  # this is what FLASH sees.
    if irltime:
        time = sim.getStepProp('timestamp')
        xlab = 'Timestamps (date/time)'
#     time = [t.to_julian_date() for t in sim.steps['timestamp']]
    else:
        time = sim.getStepProp('t')
        xlab = 'Simulation Time (s)'
    wallsteps = rollingAverage(tsteps, step=rollAvgStep)
    totalwall = np.sum(wallsteps)/3600
    dt = 1e3*sim.getStepProp('dt')  # milliseconds
    # FLASH cannot know the layout of the machine
    # so pePerNode has to be provided to calculate the nodes used.
    nodes = PE/pePerNode
    cost = 0
    # sum walltimes per submit
    for subnum in np.unique(subs):
        mask = (subs == subnum)
        walltime = np.sum(wallsteps[mask])
        nodes = np.unique(PE[mask])[0]/pePerNode
        cost += walltime*nodes/3600
    # find max number for limits of plot
    lims = [np.max(dt), np.max(PE/1000), np.max(wallsteps), np.max(subs)]
    top = 1.2*np.max(lims)
    f = plt.figure(dpi=dpi)
    layout = (3, 1)
    # main properties
    alpha = 0.9
    mainx = plt.subplot2grid(layout, (0, 0), rowspan=2)
    mainx.semilogy(time, dt, label=u'tstep(ms)', alpha=alpha)
    mainx.semilogy(time, wallsteps, label=u'IRL(s)', alpha=alpha)
    mainx.semilogy(time, PE/1000, label=u'ranks($10^3$)', alpha=alpha)
    # mainx.semilogy(time, sim.getStepProp('n')/1e6, label=u'Step($10^6$)', alpha=0.3)
    # mainx.semilogy(time, subs, label='Submit', alpha=0.3)
    mainx.set_ylabel(u'(ms, s, unitless)')
    mainx.set_ylim([bot, top])
    plt.setp(mainx.get_xticklabels(), visible=False)
    mainx.legend(**legdict)
    pc = mainx._get_lines.prop_cycler
    # blocks
    blkax = plt.subplot2grid(layout, (2, 0), sharex=mainx, rowspan=2)
    fac = 2 if np.max(totblk) < 1e3 else 3
    flab = u'Tot($10^{:d}$)'.format(fac)
    if np.min(minblk) != 0:
        blkax.semilogy(time, minblk, label='min', color=next(pc)['color'], alpha=alpha)
        blkax.semilogy(time, maxblk, label='max', color=next(pc)['color'], alpha=alpha)
        totblk /= 10**fac
        blkax.semilogy(time, totblk, label=flab, color=next(pc)['color'], alpha=alpha)
        # -5: space for tot blocks when scaled
        spacer = np.min(minblk)-5
        low = 1 if spacer <= 0 else spacer
        blkax.set_ylim([low, blkprank])
    else:
        blkax.plot(time, minblk, label='min', color=next(pc)['color'], alpha=alpha)
        blkax.plot(time, maxblk, label='max', color=next(pc)['color'], alpha=alpha)
        blkax.plot(time, totblk/10**fac, label=flab, color=next(pc)['color'], alpha=alpha)
        blkax.set_ylim([-20, blkprank])
    blkax.yaxis.set_major_formatter(StrMethodFormatter('{x:.0f}'))
    blkax.set_ylabel('Blocks')
    if not irltime:
        blkax.xaxis.set_major_formatter(StrMethodFormatter('{x:.2f}'))
    blkax.set_xlabel(xlab)
    if irltime:
        plt.subplots_adjust(bottom=0.2)
        plt.xticks(rotation=25)
    blkax.legend(**legdict)
    if showcheckpoints:
        for p in sim.chkp:
            if irltime:
                x = p.timestamp
            else:
                try:
                    loc = np.where(sim.steps['timestamp'] < p.timestamp)[0][-1]
                    x = sim.steps['t'][loc]
                except IndexError:  # chk/plt 0 are behind step 1
                    continue
            if p.filetype.split()[0]=='checkpoint':
                mainx.axvline(x, color='green', alpha=0.2, ls='--')
                blkax.text(x, 1.6*blkprank, p.number, color='green', size=8)
            else:
                mainx.axvline(x, color='black', alpha=0.2)
                blkax.text(x, 1.1*blkprank, p.number, color='black', size=8)
    tfmt = 'IRL(h):{:.2f},Subs:{:d},Node h:{:.2f}'
    stub = mainx.set_title(tfmt.format(totalwall, subs[-1], cost))
    f.tight_layout()
    f.subplots_adjust(hspace=0.0)
    return f


def plotStats(sim, xlims=[1e-8, 6.0], ylims=[1e-16, 3], logx=True):
    """build a figure with energy components vs sim time."""
    f, ax = plt.subplots()
    if logx:
        ax.loglog(sim.time, sim.getStepProp('E_total')/_foe,
                  label='Total', color='k')
        ax.loglog(sim.time, sim.getStepProp('E_kinetic')/_foe,
                  label='Kin', color='r', alpha=0.8, ls=':')
        ax.loglog(sim.time, sim.getStepProp('E_internal')/_foe,
                  label='Pot', color='b', alpha=0.8, ls=':')
    else:
        ax.semilogy(sim.time, sim.getStepProp('E_total')/_foe,
                    label='Total', color='k')
        ax.semilogy(sim.time, sim.getStepProp('E_kinetic')/_foe,
                    label='Kin', color='r', alpha=0.8, ls=':')
        ax.semilogy(sim.time, sim.getStepProp('E_internal')/_foe,
                    label='Pot', color='b', alpha=0.8, ls=':')
    ax.legend()
    ax.set_ylabel('Energy (foe)')
    ax.set_xlabel('Time (s)')
    ax.set_xlim(xlims)
    ax.set_ylim(ylims)
    f.tight_layout()
    return f


# def plotSunburstTiming(sim, which=0, avgProc=True, column=0):
#     """Plot a nested pie chart with timings nested outwards
#     # setup plotly from notebook:
#     import plotly
#     plotly.tools.set_credentials_file(username='NAME', api_key='APKEY')
#     avgProc False cols:
#         time sec  num calls   secs avg  time pct
#     avgProc True cols:
#         max/proc (s)  min/proc (s) avg/proc (s)   num calls

#     Args:
#         sim(flashy.sim): simulation object to extract data from.
#         which(int): which timing to plot.
#         avgProc(bool): averaged procs (True) or only proc 0 (False).
#         column(int): column from timing to plot.

#     Returns:
#         (mpl.fig)

#     """
#     time, steps, colnames, data = readTiming(sim.timings[which],
#                                              avgProc=avgProc)
#     # split block and get some shape info from it
#     units, deps, fulldata = zip(*data)
#     units = [u.strip() for u in units]
#     colname = colnames[column+1]
#     allvals = [v[column] for v in fulldata]
#     delt = np.diff(deps)
#     delt = np.append(0, delt)
#     # find the families
#     roots = [i for i, d in enumerate(deps) if d == 0]
#     rootblockGen = pairGen(roots+[len(deps)])
#     slices = []
#     for a, b in rootblockGen:
#         slices.append(slice(a, b))
#     # sort the data
#     labels, parents, ids = [], [], []
#     values, metav = [], []
#     for i, slc in enumerate(slices):
#         names = units[slc]
#         itvals = allvals[slc]
#         sdeps = deps[slc]
#         delts = delt[slc]
#         # traverse delts extracting parent trees
#         currp = [""]
#         for j, n in enumerate(delts):
#             if n > 0:
#                 currp.append(names[j-n].strip())
#             elif n < 0:
#                 for k in range(abs(n)):
#                     del currp[-1]
#             par = '-'.join(currp).strip('-')
#             id = '-'.join(currp[1:]+[names[j]]).strip('-')
#             parents.append(par)
#             ids.append(id)
#         maxval = np.max(itvals)
#         normalization = np.array(itvals)/maxval*100
#         normalization[0] = maxval/time*100
#         values = values + itvals
#         metav = metav + list(normalization)
#         labels = labels + names
#     df = pd.DataFrame(data={'labels': labels, 'parents': parents,
#                             'vals': values, 'meta': metav,
#                             'ids': ids})
#     trace = gob.Sunburst(ids=df.ids, labels=df.labels, values=df.vals,
#                          parents=df.parents, meta=df.meta,
#                          maxdepth=5, branchvalues='remainder',
#                          hovertemplate="<b>%{label}</b><br>" +
#                                        "%{value:.0f}<br>" +
#                                        "%{meta:.4f} <extra></extra>")
#     tsize = 16  # making space for title
#     forstr = "Var:{}<br>Total:{:.2f} Steps:{:d}"
#     tag = forstr.format(colname, time, steps)
#     layout = gob.Layout(margin=gob.layout.Margin(t=tsize*5, l=0, r=0, b=0),
#                         sunburstcolorway=["#636efa", "#EF553B", "#00cc96",
#                                           "#ab63fa", "#19d3f3", "#e763fa",
#                                           "#FECB52", "#FFA15A", "#FF6692"],
#                         extendsunburstcolors=True,
#                         title=tag)
#     fig = gob.Figure(data=[trace], layout=layout)
#     return fig


# def plotPizzaTiming(sim, which=0, avgProc=True, column=0, clrshift=1):
#     """Plot a nested pie chart with timings nested outwards
#     # setup plotly from notebook:
#     import plotly
#     plotly.tools.set_credentials_file(username='NAME', api_key='APKEY')
#     avgProc False cols:
#         time sec  num calls   secs avg  time pct
#     avgProc True cols:
#         max/proc (s)  min/proc (s) avg/proc (s)   num calls

#     Args:
#         sim(flashy.sim): simulation object to extract data from.
#         which(int): which timing to plot.
#         avgProc(bool): averaged procs (True) or only proc 0 (False).
#         column(int): column from timing to plot.
#         clrshift(int): change color list start (change colorscheme).

#     Returns:
#         (mpl.fig)

#     """
#     if not sim.timings:
#         print("No timing information.")
#         return None
#     try:
#         init_notebook_mode(connected=True)
#     except:
#         print('Plotly not configured, see docstring.')
#         return None
#     time, steps, colnames, data = readTiming(sim.timings[which],
#                                              avgProc=avgProc)
#     # split block and get some shape info from it
#     units, deps, fulldata = zip(*data)
#     colname = colnames[column+1]
#     uts = [v[column] for v in fulldata]
#     maxdep = max(deps)+1
#     delt = np.diff(deps)
#     delt = np.append(0, delt)
#     roots = [i for i, d in enumerate(deps) if d == 0]
#     rootblockGen = pairGen(roots+[len(deps)])
#     slices = []
#     for a, b in rootblockGen:
#         slices.append(slice(a, b))
#     traces, mcolors = [], []
#     shift = 0
#     for i, slc in enumerate(slices):
#         names = units[slc]
#         values = uts[slc]
#         sdeps = deps[slc]
#         delts = delt[slc]
#         drange = range(1, np.max(sdeps)+1)
#         if delts[0] < 0:
#             delts[0] = 0
#         meta = []
#         parents = [names[0].strip()]
#         for j, n in enumerate(delts):
#             if n > 0:
#                 parents.append(names[j-n].strip())
#             elif n < 0:
#                 for k in range(abs(n)):
#                     del parents[-1]
#             meta.append(parents[-1].strip())
#         mcolors.append(colors[i+clrshift])
#         if not len(drange):
#             shift += 1
#             continue
#         for m, d in enumerate(drange):
#             subvals = [values[k] for k, j in enumerate(sdeps) if j == d]
#             subnames = [names[k].strip()
#                         for k, j in enumerate(sdeps) if j == d]
#             submeta = [meta[k] for k, j in enumerate(sdeps) if j == d]
#             grad = linear_gradient(colors[i+clrshift], n=len(meta))
#             trace = gob.Pie(values=subvals, labels=subnames,
#                             sort=True, direction='clockwise', opacity=0.9,
#                             text=submeta, hole=None,
#                             hovertemplate="<b>%{label}</b><br>" +
#                                           "%{text}<br>" +
#                                           "%{percent}<br><extra></extra>",
#                             title=dict(text=parents[0],
#                                        font=dict(color='black', size=20),
#                                        position='bottom center'),
#                             domain={'row': i-shift, 'column': m},
#                             textinfo='value', rotation=20,
#                             textposition='inside',
#                             marker={'colors': grad['hex'],
#                                     'line': {'color': 'black',
#                                              'width': 1.0}})
#             traces.append(trace)
#     # finally overplot depth 0
#     hnames = [units[i] for i, d in enumerate(deps) if d == 0]
#     hvals = [uts[i] for i, d in enumerate(deps) if d == 0]
#     htrace = gob.Pie(values=hvals, labels=hnames,
#                      sort=True, direction='clockwise', opacity=1.0,
#                      hole=0.2, hoverinfo='label+value',
#                      title=dict(text='Root Calls',
#                                 font=dict(color='black', size=20),
#                                 position='top center'),
#                      domain=dict(row=0, column=maxdep-2),
#                      textinfo='percent', rotation=20,
#                      textposition='inside',
#                      marker={'colors': mcolors,
#                              'line': {'color': 'black', 'width': 1.0}})
#     traces.append(htrace)
#     cstr = "Variable: {}<br>".format(colname)
#     tstr = " Total Time:{:.4f}<br>".format(time)
#     sstr = "Steps:{:d}".format(steps)
#     layout = gob.Layout(showlegend=False,
#                         grid=dict(rows=len(slices)-shift, columns=5),
#                         title=dict(text=cstr + tstr + sstr,
#                                    font=dict(color='black', size=20)),
#                         autosize=False,
#                         height=250*(i+1-shift), width=200*(maxdep-1))
#     fig = gob.FigureWidget(data=traces, layout=layout)
#     return fig


def homologateRefsToBulk(sim, property):
    baseSteps = sim.getStepProp('n')
    refSteps = sim.getStepProp('n', src='refs')
    extended = np.zeros(len(sim.steps))
    refProp = sim.getStepProp(property, src='refs')
    for i, n in enumerate(refSteps):
        stencil = np.where(baseSteps-n > 0)
        if not i:
            extended[:] = refProp[0]
        else:
            extended[stencil] = refProp[i]
    return extended
