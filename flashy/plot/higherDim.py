from .globals import (np, os, _ldpi, AxesGrid, plt, log, reformatTag, 
                      SymLogNorm, ScalarFormatter, StrMethodFormatter,
                      writeFig, resizeText, ImageGrid)
from flashy.datahaul.hdf5yt import getFields, yt
from flashy.IOutils import pairGen, getFileList
from flashy.datahaul.plainText import dataMatrix
from yt.funcs import mylog  # avoid yt warnings
import flashy.datahaul.ytfields as ytf
import flashy.datahaul.ytDatasetTools as ytt
import flashy.datahaul.tmat as tmat
from flashy.plot.oneDim import plotSpecies
from mpl_toolkits.axes_grid1.axes_grid import CbarAxes
mylog.setLevel(50)
log.name = __name__


def slice_cube(fname, grids=False, batch=False,
               frame=1e9, center=[0.0, 0.0, 0.0], sliceAxis='z',
               fields=['density', 'pressure', 'temperature'],
               linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
               maxs=[6e7, 3e+25, 8e9], mark=[], cmap='', linthresh=1e14,
               particles=False, returnArgDict=False, pwidth=1e6, comm='', 
               filetag='ytprops', dpi=80, lineouts=[],
               metaprops={'radius': 4e8}):
    """YT 2D plots of a specified list of fields through a defined slice.
    WARN: does not account for geometries other than 3D cartesian.

    Args:
        fname(str): filename to plot.
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float or tuple list): total width of frame square or 
            list of tuples specifying frame.
        center(float list): center for slice.
        sliceAxis(str): normal axis for slice.
        fields(list of str): list of named fields to plot.
        linear(bool list): set linear or log scale(false) for each plot.
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        mark(float list): mark a (x,y) coordinate in the plot.
        cmap(str): matplotlib colormap for the plot.
        linthresh(float): symlog linear area around 0.
        particles(bool): overplot particles on the leftmost plot.
        returnArgDict(bool): print locals() turned into a batch submit.
        pwidth(float): particle sliver thickness.
        comm(str): comment for figure (appears on top)
        filetag(str): override default filetag (ytprops).
        dpi(int): force dpi of figure
        lineouts(list of tuples): point-pair list to draw lineouts.

    Returns:
        (mpl.figure or None)

    """
    if returnArgDict:
        kwargs = locals()
        kwargs['batch'] = True
        del(kwargs['fname'])
        del(kwargs['returnArgDict'])
        print(kwargs)
        return None
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(7*size, 7), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # make a box because reg does not inherit derived fields
    if isinstance(frame, float):
        left_edge = [center[0]-frame, center[1] - frame,
                    center[2] - frame]
        right_edge = [center[0] + frame, center[1] + frame,
                      center[2] + frame]
        xs = slice(center[0] - frame, center[0] + frame)
        ys = slice(center[1] - frame, center[1] + frame)
        zs = slice(center[2] - frame, center[2] + frame)
    else:
        left_edge = [frame[0][0], frame[1][0], frame[2][0]]
        right_edge = [frame[0][1], frame[1][1], frame[2][1]]
        xs = slice(frame[0][0], frame[0][1])
        ys = slice(frame[1][0], frame[1][1])
        zs = slice(frame[2][0], frame[2][1])
    try:
        reg = ds.box(left_edge, right_edge)
    except Exception as e:
        reg = ds.r[xs, ys]
# get metadata
    metatxt = []
    # in situ calculations and particle data
    if batch:
        # make a box because reg does not inherit derived fields
        ad = reg
        meta = getattr(ytf, '_speed')
        ds.add_field(("flash", 'speed'), function=getattr(ytf, 'speed'), **meta)
        # ad = ds.all_data()
        mvpos = ad['speed'].argmax()
        maxsp = ad['speed'].value[mvpos]
        mvlocx, mvlocy = ad['x'].value[mvpos], ad['y'].value[mvpos]
        delT = ds.parameters['dt']
        mline = '{:.10E} {:.10E}'.format(float(ds.current_time), delT)
        metatxt.append(mline)
        mline = '{:.10E} {:.10E} {:.10E}'.format(maxsp, mvlocx, mvlocy)
        metatxt.append(mline)
        # get energy release
        log.warning('Assuming cartesian 3D slice')
        dx = ad[("flash", "path_element_x")].value
        dy = ad[("flash", "path_element_y")].value
        dz = ad[("flash", "path_element_z")].value
        vol = dx*dy*dz
        cell_masses = vol*ad['density'].value
        energyRelease = ad['enuc'].value*cell_masses  # erg/s
        # lump sum of positive energy
        energymask = energyRelease > 0.0
        lump = np.sum(energyRelease[energymask]*delT)
        mline = '{:.10E} {:.10E}'.format(lump, 0.0)
        metatxt.append(mline)
        # calculate energy maxima for both core and shell
        mline = '{:.10E} {:.10E}'
        mline = mline.format(metaprops['radius'], metaprops['radius'])
        metatxt.append(mline)
        rads = np.sqrt(ad['x'].value**2 + ad['y'].value**2 + ad['z'].value**2)
        # rads = ad['radius'].value  # these values are slightly different (?)
        # SHELL
        mask = rads >= metaprops['radius']
        mepos = ad['enuc'][mask].argmax()
        medens, metemp = ad['dens'][mask].v[mepos], ad['temp'][mask].v[mepos]
        melocxS, melocyS = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
        meloczS = ad['z'][mask].v[mepos]
        mevelx, mevely = ad['velx'][mask].v[mepos], ad['vely'][mask].v[mepos]
        mevelz = ad['velz'][mask].v[mepos]
        emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
        toten = np.sum(energyRelease[mask & energymask]*delT)
        cvol = vol[mask][mepos]
        megamc = ad['gamc'][mask].v[mepos]
        mepres = ad['pres'][mask].v[mepos]
        mecs = np.sqrt(megamc*mepres/medens)
        # emax in erg/cc/s
        mline = 'shell ' + '{:.10E} '*12
        mline = mline.format(toten, emax, melocxS, melocyS, meloczS,
                             medens, metemp, mevelx, mevely, mevelz,
                             cvol, mecs)
        metatxt.append(mline)
        # CORE
        mask = rads < metaprops['radius']
        try:
            mepos = ad['enuc'][mask].argmax()
            medens = ad['dens'][mask].v[mepos]
            metemp = ad['temp'][mask].v[mepos]
            melocxC, melocyC = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
            meloczC = ad['z'][mask].v[mepos]
            mevelx = ad['velx'][mask].v[mepos]
            mevely = ad['vely'][mask].v[mepos]
            mevelz = ad['velz'][mask].v[mepos]
            emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
            toten = np.sum(energyRelease[mask & energymask]*delT)
            cvol = vol[mask][mepos]
            megamc = ad['gamc'][mask].v[mepos]
            mepres = ad['pres'][mask].v[mepos]
            mecs = np.sqrt(megamc*mepres/medens)
        except ValueError:
            # eventually all the domain is below thresh so set this to "shell"
            melocxC, melocyC, meloczC = melocxS, melocyS, meloczS
        mline = 'core ' + '{:.10E} '*12
        mline = mline.format(toten, emax, melocxC, melocyC, meloczC,
                             medens, metemp, mevelx, mevely, mevelz,
                             cvol, mecs)
        metatxt.append(mline)
        for f in fields:
            fstr = '{} range: {:.10E} {:.10E}'
            exts = fstr.format(f, *ad.quantities.extrema(f).value)
            log.debug(exts)
            metatxt.append(exts)
# build yt plot
    p = yt.SlicePlot(ds, sliceAxis, list(fields),
                     origin='native', center=center, data_source=reg)
    p.set_axes_unit('km')
    dnum = {'x': 0, 'y': 1, 'z': 2}
    header = '{:2.8f}s, {}={:.2E}'
    axes = {
        'x': ('y', 'z'),
        'y': ('z', 'x'),
        'z': ('x', 'y')
    }.get(sliceAxis, ('x', 'y'))
    if isinstance(frame, float):
        p.set_width((frame, frame))
    else:
        rectangle = []
        for a in axes:
            pick = frame[dnum[a]]
            delt = pick[1] - pick[0]
            rectangle.append(delt)
        p.set_width(tuple(rectangle))
    if mark:
        xm, ym = np.array(mark)/1e5
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color': 'green', 's': 30,
                                     'facecolors': "None"})
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs, linear)
    for i, (f, mi, mx, lin) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(float(ds.current_time),
                                           sliceAxis.capitalize(),
                                           center[dnum[sliceAxis]]))
        if cmap:
            p.set_cmap(f, cmap)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            if mi < 0 or mx < 0:
                log.debug('{} with limit <0: {:E} {:E}'.format(f, mi, mx))
                if not np.nanmax(reg[f].v):
                    log.debug('{} only zero/nan, setting linear.'.format(f))
                    # there's only 0 or nans so fall back to linear scale
                    p.set_zlim(f, 0, 1)
                    p.set_log(f, False)
                else:
                    p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    if lineouts:
        pargs = {'color': 'black', 'ls': '--', 'linewidth': 1}
        for pos in lineouts:
            p0, p1 = np.array(pos[0])/1e5, np.array(pos[1])/1e5
            p.annotate_line(p0, p1, coord_system='plot', plot_args=pargs)
    if batch:
        # only tag the shell value
        grp = [melocxS, melocyS, meloczS]
        x, y = grp[dnum[axes[0]]]/1e5, grp[dnum[axes[1]]]/1e5
        # x, y = melocxS/1e5, melocyS/1e5
        log.warning('shell max enuc (red): {:E} {:E}'.format(x, y))
        # p.annotate_marker((x, y), coord_system='plot', marker='o',
        #                   plot_args={'color': 'yellow', 's': 250,
        #                              'linewidth': 2, 'facecolors': "None"})
        # p.annotate_sphere([0.0, 0.0, 0.0], 
        #                   radius=(metaprops['radius']/1e5, 'km'),
        #                   circle_args={'color': 'white',
        #                                'ls': '--', 'linewidth': 1})
    p._setup_plots()
    
    
    for i, ax in enumerate(fig.axes):
        if not i:
            fig.suptitle(comm, va='bottom')
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            # particles are their own thing. so only plot for checkpoints
            if particles:
                if 'plt' in fname:
                    log.warning("Filename not a checkpoint, skipping particles")
                else:
                    ad = ds.all_data(fields=[('all', 'particle_position')])
                    allp = ad[('all', 'particle_position')].v
                    ax2num = {'x': 0, 'y': 1, 'z': 2}
                    n = ax2num[sliceAxis]
                    refh = center[n]
                    subs = [p for p in allp if p[n] < (refh + 0.5*pwidth)]
                    subs = [p for p in subs if p[n] > (refh - 0.5*pwidth)]
                    lm = 'Total/Drawn Particles: {} {}'.format(len(allp), len(subs))
                    log.debug(lm)
                    prx = [p[ax2num[axes[0]]] for p in subs]
                    pry = [p[ax2num[axes[1]]] for p in subs]
                    ax.scatter(prx, pry, marker='o', s=10, c='green')
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def bulk_energy(fname, grids=False,
                frame=1e9, center=[0.0, 0.0, 0.0], sliceAxis='z',
                fields=['temperature'], linear=[], mins=[], maxs=[], mark=[], cmap='',
                linthresh=1e14, particles=False, returnArgDict=False,
                pwidth=1e6, comm='', filetag='energy_meta', dpi=80,
                metaprops={'radius': 4e8}):
    """Copy of energy_meta, only calculates bulk energy.
    avoids any large arrays
    
    Args: **see slice_cube above**
    
    Returns:
        (None)

    """
    batch=True
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(7*size, 7), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # make a box because reg does not inherit derived fields
    if isinstance(frame, float):
        left_edge = [center[0] - frame, center[1] - frame,
                    center[2] - frame]
        right_edge = [center[0] + frame, center[1] + frame,
                      center[2] + frame]
        xs = slice(center[0] - frame, center[0] + frame)
        ys = slice(center[1] - frame, center[1] + frame)
        zs = slice(center[2] - frame, center[2] + frame)
    else:
        left_edge = [frame[0][0], frame[1][0], frame[2][0]]
        right_edge = [frame[0][1], frame[1][1], frame[2][1]]
        xs = slice(frame[0][0], frame[0][1])
        ys = slice(frame[1][0], frame[1][1])
        zs = slice(frame[2][0], frame[2][1])
    try:
        reg = ds.box(left_edge, right_edge)
    except Exception as e:
        reg = ds.r[xs, ys]

# get metadata
    metatxt = []
    # in situ calculations and particle data
    if batch:
        # make a box because reg does not inherit derived fields
        ad = reg
        # ad = ds.all_data()
        delT = ds.parameters['dt']
        mline = '{:.10E} {:.10E}'.format(float(ds.current_time), delT)
        metatxt.append(mline)
        mline = '{:.10E} {:.10E} {:.10E}'.format(0.0, 0.0, 0.0)
        metatxt.append(mline)
        # get energy release
        log.warning('Assuming cartesian 3D slice')
        # cartesian calculates volumes/masses correctly
        cell_masses = ad[('gas', 'cell_mass')].value
        energyRelease = ad[('flash', 'enuc')].value*cell_masses  # erg/s
        # lump sum of positive energy
        energymask = energyRelease > 0.0
        lump = np.sum(energyRelease[energymask]*delT)
        mline = '{:.10E} {:.10E}'.format(lump, 0.0)
        metatxt.append(mline)
        # calculate energy maxima for both core and shell
        mline = '{:.10E} {:.10E}'
        mline = mline.format(metaprops['radius'], metaprops['radius'])
        metatxt.append(mline)
        rads = ad[('index', 'radius')]
        # rads = np.sqrt(ad['x'].value**2 + ad['y'].value**2 + ad['z'].value**2)
        # rads = ad['radius'].value  # these values are slightly different (?)
        # SHELL
        mline = 'shell ' + '{:.10E} '*12
        mvals = [0.0]*12
        mline = mline.format(*mvals)
        metatxt.append(mline)
        # CORE
        mline = 'core ' + '{:.10E} '*12
        mvals = [0.0]*12
        mline = mline.format(*mvals)
        metatxt.append(mline)
        for f in fields:
            fstr = '{} range: {:.10E} {:.10E}'
            exts = fstr.format(f, *[0.0, 0.0])
            log.debug(exts)
            metatxt.append(exts)
    # skip building the plot, only do meta
    writeFig(fig, os.path.join(ds.fullpath, ds.basename),
             filetag, meta='\n'.join(metatxt))


def energy_meta(fname, grids=False,
                frame=1e9, center=[0.0, 0.0, 0.0], sliceAxis='z',
                fields=['temperature'], linear=[], mins=[], maxs=[], mark=[], cmap='',
                linthresh=1e14, particles=False, returnArgDict=False,
                pwidth=1e6, comm='', filetag='energy_meta', dpi=80,
                metaprops={'radius': 4e8}):
    """Copy of slice_cube, using a dummy image to calculate whole domain
    meta properties. Only works on batch mode.
    
    Args: **see slice_cube above**
    
    Returns:
        (None)

    """
    batch=True
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(7*size, 7), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # make a box because reg does not inherit derived fields
    if isinstance(frame, float):
        left_edge = [center[0] - frame, center[1] - frame,
                    center[2] - frame]
        right_edge = [center[0] + frame, center[1] + frame,
                      center[2] + frame]
        xs = slice(center[0] - frame, center[0] + frame)
        ys = slice(center[1] - frame, center[1] + frame)
        zs = slice(center[2] - frame, center[2] + frame)
    else:
        left_edge = [frame[0][0], frame[1][0], frame[2][0]]
        right_edge = [frame[0][1], frame[1][1], frame[2][1]]
        xs = slice(frame[0][0], frame[0][1])
        ys = slice(frame[1][0], frame[1][1])
        zs = slice(frame[2][0], frame[2][1])
    try:
        reg = ds.box(left_edge, right_edge)
    except Exception as e:
        reg = ds.r[xs, ys]

# get metadata
    metatxt = []
    # in situ calculations and particle data
    if batch:
        # make a box because reg does not inherit derived fields
        ad = reg
        meta = getattr(ytf, '_speed')
        ds.add_field(("flash", 'speed'), function=getattr(ytf, 'speed'), **meta)
        # ad = ds.all_data()
        mvpos = ad['speed'].argmax()
        maxsp = ad['speed'].value[mvpos]
        mvlocx, mvlocy = ad['x'].value[mvpos], ad['y'].value[mvpos]
        delT = ds.parameters['dt']
        mline = '{:.10E} {:.10E}'.format(float(ds.current_time), delT)
        metatxt.append(mline)
        mline = '{:.10E} {:.10E} {:.10E}'.format(maxsp, mvlocx, mvlocy)
        metatxt.append(mline)
        # get energy release
        log.warning('Assuming cartesian 3D slice')
        dx = ad[("flash", "path_element_x")].value
        dy = ad[("flash", "path_element_y")].value
        dz = ad[("flash", "path_element_z")].value
        vol = dx*dy*dz
        cell_masses = vol*ad['density'].value
        energyRelease = ad['enuc'].value*cell_masses  # erg/s
        # lump sum of positive energy
        energymask = energyRelease > 0.0
        lump = np.sum(energyRelease[energymask]*delT)
        mline = '{:.10E} {:.10E}'.format(lump, 0.0)
        metatxt.append(mline)
        # calculate energy maxima for both core and shell
        mline = '{:.10E} {:.10E}'
        mline = mline.format(metaprops['radius'], metaprops['radius'])
        metatxt.append(mline)
        rads = np.sqrt(ad['x'].value**2 + ad['y'].value**2 + ad['z'].value**2)
        # rads = ad['radius'].value  # these values are slightly different (?)
        # SHELL
        mask = rads >= metaprops['radius']
        mepos = ad['enuc'][mask].argmax()
        medens, metemp = ad['dens'][mask].v[mepos], ad['temp'][mask].v[mepos]
        melocxS, melocyS = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
        meloczS = ad['z'][mask].v[mepos]
        mevelx, mevely = ad['velx'][mask].v[mepos], ad['vely'][mask].v[mepos]
        mevelz = ad['velz'][mask].v[mepos]
        emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
        toten = np.sum(energyRelease[mask & energymask]*delT)
        cvol = vol[mask][mepos]
        megamc = ad['gamc'][mask].v[mepos]
        mepres = ad['pres'][mask].v[mepos]
        mecs = np.sqrt(megamc*mepres/medens)
        # emax in erg/cc/s
        mline = 'shell ' + '{:.10E} '*12
        mline = mline.format(toten, emax, melocxS, melocyS, meloczS,
                             medens, metemp, mevelx, mevely, mevelz,
                             cvol, mecs)
        metatxt.append(mline)
        # CORE
        mask = rads < metaprops['radius']
        try:
            mepos = ad['enuc'][mask].argmax()
            medens = ad['dens'][mask].v[mepos]
            metemp = ad['temp'][mask].v[mepos]
            melocxC, melocyC = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
            meloczC = ad['z'][mask].v[mepos]
            mevelx = ad['velx'][mask].v[mepos]
            mevely = ad['vely'][mask].v[mepos]
            mevelz = ad['velz'][mask].v[mepos]
            emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
            toten = np.sum(energyRelease[mask & energymask]*delT)
            cvol = vol[mask][mepos]
            megamc = ad['gamc'][mask].v[mepos]
            mepres = ad['pres'][mask].v[mepos]
            mecs = np.sqrt(megamc*mepres/medens)
        except ValueError:
            # eventually all the domain is below thresh so set this to "shell"
            melocxC, melocyC, meloczC = melocxS, melocyS, meloczS
        mline = 'core ' + '{:.10E} '*12
        mline = mline.format(toten, emax, melocxC, melocyC, meloczC,
                             medens, metemp, mevelx, mevely, mevelz,
                             cvol, mecs)
        metatxt.append(mline)
        for f in fields:
            fstr = '{} range: {:.10E} {:.10E}'
            exts = fstr.format(f, *ad.quantities.extrema(f).value)
            log.debug(exts)
            metatxt.append(exts)
    # skip building the plot, only do meta
    writeFig(fig, os.path.join(ds.fullpath, ds.basename),
             filetag, meta='\n'.join(metatxt))


def mprops_split(fname, mhead=False, grids=False, batch=False, frame=1e9,
                 center=(0.0, 0.0), returnArgDict=False,
                 fields=['density', 'pressure', 'temperature'],
                 linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                 maxs=[6e7, 3e+25, 8e9], marks=[], lineouts=[],
                 cmaps=['RdBu_r']*3, linthresh=1e15, dpi=90, comm='', 
                 filetag='ytm_radius', display=False, particles=False,
                 metaprops={'radius': 4e8, 'alphasminfactor': 1e-4}):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis.
    Batch builds a meta plaintext file with a core/shell split
    based only on geometric radius.

    time timedelta
    maxspeed xlocation ylocation
    # integrated enuc over a threshold
    (sum of positive energy release) threshold
    # peak energies for shell and core, discerning by limit radius
    radius radius
    'shell' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely cell_vol soundspeed
    'core' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely cell_vol soundspeed
    (listed plot variable ranges)

    Args:
        fname(str): filename to plot.
        mhead(bool): mark the position of the matchhead.
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): physical extent of plot in x (twice for y) .
        center(float tuple): physical center of plot.
        fields(str list): list of named fields to plot.
        linear(bool list): set linear or log scale(false).
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        marks(float list): mark a list of (x,y) coordinate.
        lineouts(float list): list of lineouts by point pairs.
        cmaps(str): matplotlib colormap for the plot.
        linthresh(float): symlog linear area around 0.
        dpi(float): dpi of figure returned/saved.
        filetag(str): base filename for ouput (if written).
        comm(str): add custom mesage to plot.
        display(bool): plot without meta markers, skipping calculation.
        returnArgDict(bool): print locals() turned into a batch submit.
        metaprops(dict): brittle meta file properties.
            'radius': radius for shell/core separation.
            'alphasminfactor': fraction of max to set zlim bottom.
                negative values override this to set limits by hand.

    Returns:
        (mpl.figure or None)

    """
    if returnArgDict:
        kwargs = locals()
        kwargs['batch'] = True
        kwargs['display'] = False
        del(kwargs['fname'])
        del(kwargs['returnArgDict'])
        print(kwargs)
        return None
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # in situ calculations and particle data
    if not display or particles:
        ad = ds.all_data()
    if particles:
        allp = ad[('all', 'particle_position')].v
    if not display:
        mvpos = ad['speed'].argmax()
        maxsp = ad['speed'].value[mvpos]
        mvlocx, mvlocy = ad['x'].value[mvpos], ad['y'].value[mvpos]
        delT = ds.parameters['dt']
        mline = '{:.10E} {:.10E}'.format(float(ds.current_time), delT)
        metatxt.append(mline)
        mline = '{:.10E} {:.10E} {:.10E}'.format(maxsp, mvlocx, mvlocy)
        metatxt.append(mline)
        # get energy release
        log.warning('Assuming cylindrical slice')
        dx = ad['path_element_x'].value
        dy = ad['path_element_y'].value
        r = ad['x'].value
        cylvol = 2.0*np.pi*dy*dx*r
        cell_masses = cylvol*ad['density'].value
        energyRelease = ad['enuc'].value*cell_masses  # erg/s
        # lump sum of positive energy
        energymask = energyRelease > 0.0
        lump = np.sum(energyRelease[energymask]*delT)
        mline = '{:.10E} {:.10E}'.format(lump, 0.0)
        metatxt.append(mline)
        # calculate energy maxima for both core and shell
        mline = '{:.10E} {:.10E}'
        mline = mline.format(metaprops['radius'], metaprops['radius'])
        metatxt.append(mline)
        # SHELL
        rads = np.sqrt(ad['x'].value**2 + ad['y'].value**2)
        mask = rads >= metaprops['radius']
        mepos = ad['enuc'][mask].argmax()
        medens, metemp = ad['dens'][mask].v[mepos], ad['temp'][mask].v[mepos]
        melocxS, melocyS = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
        mevelx, mevely = ad['velx'][mask].v[mepos], ad['vely'][mask].v[mepos]
        emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
        toten = np.sum(energyRelease[mask & energymask]*delT)
        cvol = cylvol[mask][mepos]
        megamc = ad['gamc'][mask].v[mepos]
        mepres = ad['pres'][mask].v[mepos]
        mecs = np.sqrt(megamc*mepres/medens)
        # emax in erg/cc/s
        mline = 'shell ' + '{:.10E} '*10
        mline = mline.format(toten, emax, melocxS, melocyS,
                             medens, metemp, mevelx, mevely,
                             cvol, mecs)
        metatxt.append(mline)
        # CORE
        mask = rads < metaprops['radius']
        try:
            mepos = ad['enuc'][mask].argmax()
            medens = ad['dens'][mask].v[mepos]
            metemp = ad['temp'][mask].v[mepos]
            melocxC, melocyC = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
            mevelx = ad['velx'][mask].v[mepos]
            mevely = ad['vely'][mask].v[mepos]
            emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
            toten = np.sum(energyRelease[mask & energymask]*delT)
            cvol = cylvol[mask][mepos]
            megamc = ad['gamc'][mask].v[mepos]
            mepres = ad['pres'][mask].v[mepos]
            mecs = np.sqrt(megamc*mepres/medens)
        except ValueError:
            # eventually all the domain is below thresh so set this to "shell"
            melocxC, melocyC = melocxS, melocyS
        mline = 'core ' + '{:.10E} '*10
        mline = mline.format(toten, emax, melocxC, melocyC,
                             medens, metemp, mevelx, mevely,
                             cvol, mecs)
        metatxt.append(mline)
        for f in fields:
            fstr = '{} range: {:.10E} {:.10E}'
            exts = fstr.format(f, *ad.quantities.extrema(f).value)
            log.debug(exts)
            metatxt.append(exts)
    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    if mhead:
        x_match = ds.parameters['x_match']/1e5
        y_match = ds.parameters['y_match']/1e5
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color': 'green', 's': 100})
    if marks:
        pargs = {'color': 'green', 's': 250, 
                 'linewidth': 2, 'facecolors': "None"}
        for mark in marks:
            xm, ym = np.array(mark)/1e5
            log.warning('mark (green): {:E} {:E}'.format(xm, ym))
            p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                              plot_args=pargs)
    if not display:
        x, y = mvlocx/1e5, mvlocy/1e5
        log.warning('max speed (white): {:E} {:E}'.format(x, y))
        # p.annotate_marker((x, y), coord_system='plot', marker='o',
        #                   plot_args={'color': 'white', 's': 250,
        #                              'linewidth': 2, 'facecolors': "None"})
        x, y = melocxS/1e5, melocyS/1e5
        log.warning('shell max enuc (red): {:E} {:E}'.format(x, y))
        # p.annotate_marker((x, y), coord_system='plot', marker='o',
        #                   plot_args={'color': 'yellow', 's': 250,
        #                              'linewidth': 2, 'facecolors': "None"})
        x, y = melocxC/1e5, melocyC/1e5
        log.warning('core max enuc (red): {:E} {:E}'.format(x, y))
        # p.annotate_marker((x, y), coord_system='plot', marker='o',
        #                   plot_args={'color': 'red', 's': 250,
        #                              'linewidth': 2, 'facecolors': "None"})
    if 'spheres' in metaprops:
        for sph in metaprops['spheres']:
            cent = sph[0]
            rad = sph[1]
            pargs = {'color': 'red', 'ls': '--', 'linewidth': 1}
            p.annotate_sphere(cent, radius=(rad, 'cm'),
                              circle_args=pargs)
    p.set_font({'family': 'monospace'})
    pargs = {'color': 'white', 'ls': '--', 'linewidth': 1}
    # p.annotate_sphere([0.0, 0.0, 0.0], radius=(metaprops['radius'], 'cm'),
    #                   circle_args=pargs)
    if lineouts:
        pargs = {'color': 'black', 'ls': '--', 'linewidth': 1}
        for pos in lineouts:
            p0, p1 = np.array(pos[0])/1e5, np.array(pos[1])/1e5
            p.annotate_line(p0, p1, coord_system='plot', plot_args=pargs)
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs, linear, cmaps)
    spstags = []
    for i, (f, mi, mx, lin, cm) in enumerate(pvars):
        if cm:
            p.set_cmap(f, cm)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        if f in ytf._alphas:
            region = ds.r[0:frame, -frame:frame]
            peak = region.max(f)
            if metaprops['alphasminfactor'] < 0.0:
                p.set_zlim(f, mi, mx)
            else:
                nmi = peak*metaprops['alphasminfactor']
                p.set_zlim(f, nmi, peak)
            ft = "Max:{:.3E}".format(peak.v)
            # yt annotates over all plots so this is a 
            # deep pass to matplotlib to set the peak tag
            spstags.append((i, ft))
        else:
            p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            if mi < 0 or mx < 0:
                log.debug('{} with limit <0: {:E} {:E}'.format(f, mi, mx))
                # this is prone to error so check the subset
                cent = np.array((center[0], center[1], 0.0))
                delt = np.array((frame*0.5, frame, ds.domain_dimensions[2]))
                le = cent-delt
                re = cent+delt
                viewbox = ds.region(cent, le, re)
                if not np.nanmax(viewbox[f].v):
                    # there's only 0 or nans so fall back to linear scale
                    p.set_zlim(f, 0, 1)
                    p.set_log(f, False)
                else:
                    p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    
    if display:
        datb, datb2 = ytt.get_plotTitles(ds)
        ll = len(datb.split('\n')[-1])
        n, v, t = datb2.split()
        mid = str(ll - len(n) - 2)
        restr = '{}{:>' + mid + '} {}'
        restr = restr.format(n, v, t)
        fig.axes[0].set_title('\n'.join([datb, restr]))
    else:
        extras = {'maxspeed': maxsp}
        datb, datb2 = ytt.get_plotTitles(ds, comment=comm, extras=extras)
        fig.axes[0].set_title(datb)
        fig.axes[1].set_title(datb2)

    spstagnums = [s[0] for s in spstags]
    for i, ax in enumerate(fig.axes):
        if i in spstagnums:
            ax.annotate(spstags[spstagnums.index(i)][-1],
                        xy=(0.55, 0.95), xycoords='axes fraction')
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            if particles:
                uniq, unpos = np.unique(allp, axis=0, return_index=True)
                subset = [x for i, x in enumerate(uniq) if i in unpos]
                prx, pry, _ =list(zip(*subset))
                ax.scatter(prx, pry, marker='.', s=5, c='green')
                subset = [x for i, x in enumerate(uniq) if i not in unpos]
                if subset:
                    prx, pry, _ =list(zip(*subset))
                    ax.scatter(prx, pry, marker='.', s=7, c='yellow')
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                if metaprops['alphasminfactor'] < 0.0:
                    ax.set_ylabel(tag + ' Mass Fraction')
                else:
                    ax.set_ylabel(tag + ' (% from max)')
    if not batch:
        return fig
    else:
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def contourProps(fname, grids=False, batch=False, frame=1e9,
                 center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
                 linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                 maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
                 linthresh=1e15, dpi=90,
                 contp={'contour': 'enuc', 'maxbased':True,
                        'values': [0.006, 0.01]}):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis. utilizes contours to get data of cells
    upwards of the smallest contour.

    Args:
        fname(str): filename to plot.
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): physical extent of plot in x (twice for y) .
        center(float tuple): physical center of plot.
        fields(str list): list of named fields to plot.
        linear(bool list): set linear or log scale(false).
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        cmaps(str): matplotlib colormap for the plot.
        linthresh(float): symlog linear area around 0.
        dpi(float): dpi of figure returned/saved.
        contp(dict): contour values and field.
        maxbased(bool): contours based on max (values become fractions)

    Returns:
        (mpl.figure or None)

    """
    ccolw = ['green', 'yellow', 'red']
    wh = len(ccolw)
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    if contp['maxbased']:
        ad = ds.all_data()
        pkp = ad[contp['contour']].argmax()
        peaken = ad[contp['contour']].value[pkp]
        mline = '# ' + '{:.10E} '*2
        delT = ds.parameters['dt']
        mline = mline.format(float(ds.current_time), delT)
        metatxt.append(mline)
        cut = peaken*np.max(contp['values'])
        log.warning('Peak, cutoff: {:.3E} {:.3E}'.format(peaken, cut))
        msk = ad['enuc'].value > cut
        r = ad['x'][msk].value
        dns = ad['dens'][msk].value
        ts = ad['temp'][msk].value
        ens = ad['enuc'][msk].value
        xs, ys = ad['x'][msk].value, ad['y'][msk].value
        vxs, vys = ad['velx'][msk].value, ad['vely'][msk].value
        dxs = ad['path_element_x'][msk].value
        dys = ad['path_element_y'][msk].value
        log.warning('# cells picked: {}'.format(len(xs)))
        log.warning('Assuming cylindrical slice')
        for dn, t, en, x, y, vx, vy, dx, dy in zip(dns, ts, ens, xs, ys,
                                                   vxs, vys, dxs, dys):
            vl = 2.0*np.pi*dy*dx*x
            ms = vl*dn
            sp = np.sqrt(vx**2+vy**2)
            rd = np.sqrt(x**2+y**2)
            values = (en, x, y, rd, dn, t, vl, ms, vx, vy, sp)
            strnums = ["{:.10E}".format(k) for k in values]
            mline = ' '.join(strnums)
            metatxt.append(mline)

    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs, linear, cmaps)
    for i, (f, mi, mx, lin, cm) in enumerate(pvars):
        if cm:
            p.set_cmap(f, cm)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            if mi < 0 or mx < 0:
                log.debug('{} with limit <0: {:E} {:E}'.format(f, mi, mx))
                # this is prone to error so check the subset
                cent = np.array((center[0], center[1], 0.0))
                delt = np.array((frame*0.5, frame, ds.domain_dimensions[2]))
                le = cent-delt
                re = cent+delt
                viewbox = ds.region(cent, le, re)
                if not np.nanmax(viewbox[f].v):
                    # there's only 0 or nans so fall back to linear scale
                    p.set_zlim(f, 0, 1)
                    p.set_log(f, False)
                else:
                    p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]

    for i, v in enumerate(contp['values']):
        if contp['maxbased']:
            val = peaken*v
            print(val, peaken, ccolw[i%wh])
        p.annotate_contour(contp['contour'], ncont=1,
                           clim=(val, val), label=False,
                           text_args={'colors': 'k'},
                           plot_args={'colors': ccolw[i%wh], 
                          'linewidths': 1, 'linestyles':'-'})
    p._setup_plots()

    t = 'Contours: {}'.format(contp['contour'])
    vals = ['{:.3E}'.format(x) for x in contp['values']]
    t +='\nVals: {}'.format(','.join(vals))
    datb, datb2 = ytt.get_plotTitles(ds, comment=t)
    fig.axes[0].set_title(datb)
    fig.axes[1].set_title(datb2)

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    for i, ax in enumerate(fig.axes):
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        filetag = 'contour'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def sphereProps(fname, sphR=1e5, sphC=(0.0, 0.0, 0.0),
                grids=False, batch=False, frame=1e9,
                center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
                linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
                linthresh=1e15, dpi=90):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis.
    cuts a sphere in the domain, extracting all cell data to a metafile.
    
    Args:
        sphR(float): radius of sphere.
        sphC(float tuple): center of sphere.
        *see contourProps*
    
    Returns:
        (mpl.figure or None)

    """
    ccolw = ['green', 'yellow', 'red']
    wh = len(ccolw)
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)

    ad = ds.sphere(sphC, (sphR, 'cm'))
    metatxt.append('# time timedelta')
    mline = '# ' + '{:.10E} '*2
    delT = ds.parameters['dt']
    mline = mline.format(float(ds.current_time), delT)
    metatxt.append(mline)
    # make a mask for all cells in the sphere
    msk = ad['temp'].value > -1.0

    r = ad['x'][msk].value
    dns = ad['dens'][msk].value
    ts = ad['temp'][msk].value
    ens = ad['enuc'][msk].value
    prs = ad['pres'][msk].value
    xs, ys = ad['x'][msk].value, ad['y'][msk].value
    vxs, vys = ad['velx'][msk].value, ad['vely'][msk].value
    dxs = ad['path_element_x'][msk].value
    dys = ad['path_element_y'][msk].value
    log.warning('# cells picked: {}'.format(len(xs)))
    log.warning('Assuming cylindrical slice')
    vls = 2.0*np.pi*dys*dxs*xs
    mss = vls*dns
    sps = np.sqrt(vxs**2+vys**2)
    rds = np.sqrt(xs**2+ys**2)
    els = []
    for el in ytf._alphas:
        els.append(ad[el][msk].value)
    chungus = zip(ens, xs, ys, rds, dns, ts, vls, prs, 
                  mss, vxs, vys, sps, *els)
    mline = '# enuc x y r dens temp vol pres mass velx vely speed '
    mline += ' '.join(ytf._alphas) 
    metatxt.append(mline)
    for line in chungus:
        strnums = ["{:.10E}".format(k) for k in line]
        mline = ' '.join(strnums)
        metatxt.append(mline)

    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs, linear, cmaps)
    for i, (f, mi, mx, lin, cm) in enumerate(pvars):
        if cm:
            p.set_cmap(f, cm)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            if mi < 0 or mx < 0:
                log.debug('{} with limit <0: {:E} {:E}'.format(f, mi, mx))
                # this is prone to error so check the subset
                cent = np.array((center[0], center[1], 0.0))
                delt = np.array((frame*0.5, frame, ds.domain_dimensions[2]))
                le = cent-delt
                re = cent+delt
                viewbox = ds.region(cent, le, re)
                if not np.nanmax(viewbox[f].v):
                    # there's only 0 or nans so fall back to linear scale
                    p.set_zlim(f, 0, 1)
                    p.set_log(f, False)
                else:
                    p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p.annotate_sphere(sphC, radius=(sphR, 'cm'),
                      circle_args={'color': 'black', 'ls': '--', 'linewidth': 1})
    p._setup_plots()

    t = 'Contour R: {:.2E}'.format(sphR)
    t +='\nCenter: {:.1E},{:.1E},{:.1E}'.format(*sphC)
    datb, datb2 = ytt.get_plotTitles(ds, comment=t)
    fig.axes[0].set_title(datb)
    fig.axes[1].set_title(datb2)

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    for i, ax in enumerate(fig.axes):
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        filetag = 'sphereContour'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def zoneProps(fname, ledge=(0.0, 0.0, 0.0), redge=(40e5, 50e5, 0.0175),
              grids=False, batch=False, frame=1e9,
              center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
              linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
              maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
              linthresh=1e15, dpi=90):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis.
    cuts a box in the domain, extracting all cell data to a metafile.
    
    Args:
        ledge(float tuple): left edge of zone to focus on.
        redge(float tuple): right edge of zone.
        *see contourProps*
    
    Returns:
        (mpl.figure or None)

    """
    ccolw = ['green', 'yellow', 'red']
    wh = len(ccolw)
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)

    ad = ds.box(ledge, redge)
    metatxt.append('# time timedelta')
    mline = '# ' + '{:.10E} '*2
    delT = ds.parameters['dt']
    mline = mline.format(float(ds.current_time), delT)
    metatxt.append(mline)
    # make mask all cells in the sphere
    msk = ad['temp'].value > -1.0

    r = ad['x'][msk].value
    dns = ad['dens'][msk].value
    ts = ad['temp'][msk].value
    ens = ad['enuc'][msk].value
    prs = ad['pres'][msk].value
    xs, ys = ad['x'][msk].value, ad['y'][msk].value
    vxs, vys = ad['velx'][msk].value, ad['vely'][msk].value
    dxs = ad['path_element_x'][msk].value
    dys = ad['path_element_y'][msk].value
    log.warning('# cells picked: {}'.format(len(xs)))
    log.warning('Assuming cylindrical slice')
    vls = 2.0*np.pi*dys*dxs*xs
    mss = vls*dns
    sps = np.sqrt(vxs**2+vys**2)
    rds = np.sqrt(xs**2+ys**2)
    els = []
    for el in ytf._alphas:
        els.append(ad[el][msk].value)
    chungus = zip(ens, xs, ys, rds, dns, ts, vls, prs, 
                  mss, vxs, vys, sps, *els)
    mline = '# enuc x y r dens temp vol pres mass velx vely speed '
    mline += ' '.join(ytf._alphas) 
    metatxt.append(mline)
    for line in chungus:
        strnums = ["{:.10E}".format(k) for k in line]
        mline = ' '.join(strnums)
        metatxt.append(mline)

    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs, linear, cmaps)
    for i, (f, mi, mx, lin, cm) in enumerate(pvars):
        if cm:
            p.set_cmap(f, cm)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            if mi < 0 or mx < 0:
                log.debug('{} with limit <0: {:E} {:E}'.format(f, mi, mx))
                # this is prone to error so check the subset
                cent = np.array((center[0], center[1], 0.0))
                delt = np.array((frame*0.5, frame, ds.domain_dimensions[2]))
                le = cent-delt
                re = cent+delt
                viewbox = ds.region(cent, le, re)
                if not np.nanmax(viewbox[f].v):
                    # there's only 0 or nans so fall back to linear scale
                    p.set_zlim(f, 0, 1)
                    p.set_log(f, False)
                else:
                    p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    # draw rectangle by a loop of lines
    #        ---- redge
    #       |   |
    # Ledge ---- 
    ledge = ledge/1e5
    redge = redge/1e5
    p.annotate_line((ledge[0], ledge[1]), (redge[0], ledge[1]), coord_system='plot',
                    plot_args={'color': 'black', 'ls': '--', 'linewidth': 1})
    p.annotate_line((redge[0], ledge[1]), (redge[0], redge[1]), coord_system='plot',
                    plot_args={'color': 'black', 'ls': '--', 'linewidth': 1})
    p.annotate_line((redge[0], redge[1]), (ledge[0], redge[1]), coord_system='plot',
                    plot_args={'color': 'black', 'ls': '--', 'linewidth': 1})
    p.annotate_line((ledge[0], redge[1]), (ledge[0], ledge[1]), coord_system='plot',
                    plot_args={'color': 'black', 'ls': '--', 'linewidth': 1})
    p._setup_plots()

    t = 'Left edge: {:.1E},{:.1E},{:.1E}'.format(*ledge)
    t +='\nRight edge: {:.1E},{:.1E},{:.1E}'.format(*redge)
    datb, datb2 = ytt.get_plotTitles(ds, comment=t)
    fig.axes[0].set_title(datb)
    fig.axes[1].set_title(datb2)

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    for i, ax in enumerate(fig.axes):
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        filetag = 'squareContour'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def delt2D(bview, fname, fields=['speed', 'velx', 'vely'],
           subset=[0.0, 1e9, -1e9, 1e9], cmap='RdBu_r',
           vmin=-1e5, vmax=1e5, batch=False,
           nticks=4, wedges=8):
    """plots deltas with contiguous checkpoint/plotfile
    for a list of fields.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file path.
        fields(str list): fields to process.
        subset(float list): domain window for plotting.
        cmap(str): colormap name.
        vmin(float): minimum value.
        vmax(float): maximum value.
        batch(bool): print to file toggle.
        nticks(int): ticks for xaxis.
        wedges(int): slices for data extraction.

    Returns:
        (mpl.figure or None)

    """
    dat = []
    for f in fields:
        t0, dt, ext, dattmat = tmat.get2dPane(bview, fname, f, wedges=wedges)
        dat.append(dattmat)
    shift = '{}{:04d}'.format(fname[:-4], int(fname[-4:])+1)
    shiftDat = []
    for f in fields:
        t1, dt, ext, dattmat = tmat.get2dPane(bview, shift, f)
        shiftDat.append(dattmat)
    dt = t1 - t0
    delts = [y-x for (x, y) in zip(dat, shiftDat)]
    acc = [x/dt for x in delts]
    fig, axs = plt.subplots(nrows=1, ncols=len(fields), dpi=140,
                            sharey=True, sharex=True, constrained_layout=True)
    fig.suptitle('Acceleration ({:.3f})'.format(t1))
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(acc[i], extent=ext, cmap=cmap,
                               norm=SymLogNorm(linthresh=1e2,
                                               linscale=1.0,
                                               vmin=vmin, vmax=vmax))
        ax.set_title(fields[i])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_ticks(np.arange(0.0, np.max(subset),
                                     np.max(subset)/nticks))
    cbar = fig.colorbar(mshow, ax=[axs[:]], location='bottom', extend='both')

    if not batch:
        return fig
    else:
        filetag = '_'.join(['acc'] + fields)
        writeFig(fig, fname, filetag)


def par_delt_runs_2D(bview, fname, fields=['speed', 'velx', 'vely'],
                     compdir='', cmap='RdBu_r',
                     subset=[0.0, 1e9, -1e9, 1e9], linthresh=1e7,
                     vmin=-1e9, vmax=1e9, batch=False,
                     fac=8, wedges=10):
    """Parallel. Plots deltas in checkpoint/plotfile
    for a list of fields and two runs.
    
    """
    dat = []
    fields = ['x', 'y'] + fields
    for f in fields:
        t0, dt, ext, dattmat = tmat.get2dPane(bview, fname, f, wedges=wedges)
        dat.append(dattmat)
    log.debug('Finished first file')
    shift = os.path.join(compdir, os.path.basename(fname))
    shiftDat = []
    for f in fields:
        t1, dt, ext, dattmat = tmat.SINGLE_get2dPane(shift, f)
        shiftDat.append(dattmat)
    log.debug('Finished second file')
    if dat[0].shape == shiftDat[0].shape:
        log.debug('Shapes are equal, directly calculating deltas.')
        delts = [y-x for (x, y) in zip(dat[2:], shiftDat[2:])]
    else:
        # find out which shape is largest and build outputs
        if dat[0].shape > shiftDat[0].shape:
            refx, refy = dat[0], dat[1]
            x, y = shiftDat[0], shiftDat[1]
        else:
            refx, refy = shiftDat[0], shiftDat[1]
            x, y = dat[0], dat[1]
        delts = [np.zeros(refx.shape) for x in fields[2:]]
        for i, yval in enumerate(refy[:, 0]):
            if yval == tmat._nullval:
                for d in delts:
                    d[i] = tmat._nullval
                continue
            if yval < 0:
                row = np.where(y[:, 0] >= -yval)[0][-1]
            else:
                row = np.where(y[:, 0] >= yval)[0][-1]
            # treat first segment separately (pairGen is object agnostic)
            mask = np.where(refx[i] <= x[row][0])
            for k, d in enumerate(delts, 2):
                d[i, mask] = dat[k][row, 0] - shiftDat[k][row, 0]
            for j, (r0, r1) in enumerate(pairGen(x[row]), 1):
                mask = np.where((refx[i] >= r0) & (refx[i] <= r1))
                for k, d in enumerate(delts, 2):
                    d[i, mask] = dat[k][row, j] - shiftDat[k][row, j]
    log.debug('Deltas done. Plotting')
    cols = len(fields[2:])
    fig, axs = plt.subplots(nrows=1, ncols=cols, dpi=90,
                            sharey=True, sharex=True, constrained_layout=True)
    ext = ext/1e5
    if isinstance(axs, type(np.array)):  # restore this after fidgeting
        for i, ax in enumerate(axs.flat):
            mshow = axs[i].matshow(delts[i], extent=ext, cmap=cmap,
                                   norm=SymLogNorm(linthresh=linthresh,
                                                   linscale=0.8,
                                                   vmin=vmin, vmax=vmax))
            ax.set_title(fields[i+2])
            ax.axis(subset)
            ax.xaxis.set_ticks_position('bottom')
    else:
        mshow = axs.matshow(delts[0], extent=ext, cmap=cmap,
                            norm=SymLogNorm(linthresh=linthresh,
                            linscale=0.8, vmin=vmin, vmax=vmax))
        axs.set_title(fields[2])
        axs.axis(subset)
        axs.xaxis.set_ticks_position('bottom')

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))

    if isinstance(axs, type(np.array)):
        cbar = fig.colorbar(mshow, ax=[axs[:]],
                            location='right', extend='both')
    else:
        cbar = fig.colorbar(mshow, ax=axs, location='right', extend='both')

    log.warning('Assuming filename folder structure: {}'.format(fname))
    finn = '/'.join(fname.split('/')[-4:-2])
    jake = '/'.join(compdir.split('/')[-3:-1])
    suptitle = '{}\n{}\n'.format(finn, jake)
    suptitle += '\n{:.5f}'.format(t0)
    fig.suptitle(suptitle)

    if not batch:
        return fig
    else:
        filetag = '_'.join(['delta'] + [s.strip() for s in fields])
        writeFig(fig, fname, filetag)


def delt_runs_2D(fname, fields=['speed', 'velx', 'vely'],
                 compdir='', cmap='RdBu_r',
                 subset=[0.0, 1e9, -1e9, 1e9], linthresh=1e7,
                 vmin=-1e9, vmax=1e9, batch=False):
    """(Non-parallel). Plots deltas with contiguous checkpoint/plotfile
    for a list of fields.

    Args:
        fname(str): file path.
        fields(str list): fields to process.
        compdir(str): folderpath for comparison run.
        cmap(str): colormap name.
        subset(float list): domain window for plotting ([x0,x1,y0,y1]).
        linthresh(float): linear range around zero.
        vmin(float): minimum value.
        vmax(float): maximum value.
        batch(bool): print to file toggle.

    Returns:
        (mpl.figure or None)

    """
    dat = []
    fields = ['x', 'y'] + fields
    for f in fields:
        try:  # F4 compatibility
            t0, dt, ext, dattmat = tmat.SINGLE_get2dPane(fname, f)
            dat.append(dattmat)
        except:
            t0, dt, ext, dattmat = tmat.SINGLE_get2dPane(fname, f.strip())
            dat.append(dattmat)
    log.debug('Finished first file')
    print(os.path.basename(fname))
    shift = os.path.join(compdir, os.path.basename(fname))
    shiftDat = []
    for f in fields:
        t1, dt, ext, dattmat = tmat.SINGLE_get2dPane(shift, f)
        shiftDat.append(dattmat)
    log.debug('Finished second file')
    if dat[0].shape == shiftDat[0].shape:
        log.debug('Shapes are equal, directly calculating deltas.')
        delts = [y-x for (x, y) in zip(dat[2:], shiftDat[2:])]
    else:
        # find out which shape is largest and build outputs
        if dat[0].shape > shiftDat[0].shape:
            refx, refy = dat[0], dat[1]
            x, y = shiftDat[0], shiftDat[1]
        else:
            refx, refy = shiftDat[0], shiftDat[1]
            x, y = dat[0], dat[1]
        delts = [np.zeros(refx.shape) for x in fields[2:]]
        for i, yval in enumerate(refy[:, 0]):
            if yval == tmat._nullval:
                for d in delts:
                    d[i] = tmat._nullval
                continue
            if yval < 0:
                row = np.where(y[:, 0] >= -yval)[0][-1]
            else:
                row = np.where(y[:, 0] >= yval)[0][-1]
            # treat first segment separately (pairGen is object agnostic)
            mask = np.where(refx[i] <= x[row][0])
            for k, d in enumerate(delts, 2):
                d[i, mask] = dat[k][row, 0] - shiftDat[k][row, 0]
            for j, (r0, r1) in enumerate(pairGen(x[row]), 1):
                mask = np.where((refx[i] >= r0) & (refx[i] <= r1))
                for k, d in enumerate(delts, 2):
                    d[i, mask] = dat[k][row, j] - shiftDat[k][row, j]
    log.debug('Deltas done. Plotting')
    cols = len(fields[2:])
    fig, axs = plt.subplots(nrows=1, ncols=cols, dpi=90,
                            sharey=True, sharex=True, constrained_layout=True)
#     if isinstance(axs, type(np.array)):  # restore this after fidgeting
    for d in delts:
        print(d.shape, np.max(d), np.min(d))
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(delts[i], extent=ext/1e5, cmap=cmap,
                               norm=SymLogNorm(linthresh=linthresh,
                                               linscale=0.8,
                                               vmin=vmin, vmax=vmax))
        ax.set_title(fields[i+2])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
#     else:
#         mshow = axs.matshow(delts[0], extent=ext/1e5, cmap=cmap,
#                             norm=SymLogNorm(linthresh=linthresh,
#                             linscale=0.8, vmin=vmin, vmax=vmax))
#         axs.set_title(fields[2])
#         axs.axis(subset)
#         axs.xaxis.set_ticks_position('bottom')

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        ax.set_xlabel(u'x (km)')
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
            ax.set_ylabel(u'y (km)')

    if isinstance(axs, type(np.array)):
        cbar = fig.colorbar(mshow, ax=[axs[:]],
                            location='right', extend='both')
    else:
        cbar = fig.colorbar(mshow, ax=axs, location='right', extend='both')

    log.warning('Assuming filename folder structure: {}'.format(fname))
    finn = '/'.join(fname.split('/')[-4:-2])
    jake = '/'.join(compdir.split('/')[-3:-1])
    suptitle = '{}\n{}\n'.format(finn, jake)
    suptitle += '\n{:.5f}'.format(t0)
    fig.suptitle(suptitle)

    if not batch:
        return fig
    else:
        filetag = '_'.join(['delta'] + [s.strip() for s in fields])
        writeFig(fig, fname, filetag)


def single_custom(file):
    """copy paste from 20.plots.  to run as batch."""
    filetag = 'custom'
    field = 'density'
    vmin = 5e4
    vmax = 0.2e8
    cmap = 'twilight_r'
    frame = 0.6e9
    center = (0.3e9, 0.0e8)
    markR = 4.120000E+08
    fac = 8
    linthresh = 1e10
    prad = 0.0
    linear = False
    grid = False
    fig = plt.figure(figsize=(5, 8), dpi=100)
    grid = AxesGrid(fig, 111, nrows_ncols=(1, 1),
                    axes_pad=0.0, aspect=True, label_mode="all",
                    share_all=False, cbar_location="right",
                    cbar_mode="single", cbar_size="10%", cbar_pad="0%")

    ds = yt.load(file)
    if prad:
        rgn = ds.box(left_edge=ds.arr((0.9e8, -0.5e9, 0.0), "code_length"),
                     right_edge=ds.arr((0.5e9, 0.5e9, 0.2), "code_length"))
        p = yt.SlicePlot(ds, 'z', [field], data_source=rgn)
    else:
        p = yt.SlicePlot(ds, 'z', [field])    
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    p.set_cmap(field, cmap)
    p.set_zlim(field, vmin, vmax)
    if linear:
        p.set_log(field, False)
    else:
        if vmin < 0 or vmax < 0:
            log.debug('{} with limit <0: {:E} {:E}'.format(field, vmin, vmax))
            # this is prone to error so check the subset
            cent = np.array((0.15e9, -0.27e9, 0.0))
            delt = np.array((0.3e9, 0.6e9, ds.domain_dimensions[2]))
            le = cent-delt*0.5
            re = cent+delt*0.5
            viewbox = ds.region(cent, le, re)
            if not np.nanmax(viewbox[f].v):
                # there's only 0 or nans so fall back to linear scale
                p.set_zlim(field, 0, 1)
                p.set_log(field, False)
            else:
                p.set_log(field, True, linthresh=linthresh)
        else:
            p.set_log(field, True)
    p.annotate_sphere([0.0, 0.0, 0.0], radius=(markR, 'cm'),
                      circle_args={'color': 'white', 'ls': '--', 'linewidth': 1})
    p.annotate_contour("temperature", ncont=1, clim=(1e9,1e9), label=False,
                       text_args={'colors': 'k'},
                       plot_args={'colors': 'green', 
                                  'linewidths': 1, 'linestyles':'--'})
    plot = p.plots[field]
    plot.axes = grid[0].axes
    plot.cax = grid.cbar_axes[0]
    p._setup_plots()
    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
        if not i:
            ax.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
            timetag = '{:.3f}s'.format(float(ds.current_time))
            ax.set_title(timetag, loc='right')
        elif not isinstance(ax, CbarAxes):
            ax.set_ylabel('')    
        # rotate element tags to the horizontal
        if i:  # skip the leftmost ylabel
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})

    filepath = os.path.join(ds.fullpath, ds.basename)
    writeFig(fig, filepath, filetag)


def delta_yt(filepair, field='temperature', tags=['A', 'B'],
             resolution=3000, width=6000e5, trim=9e8, vmin=-100,
             vmax=100, cmap='tab10', batch=False, dpi=120,
             filetag='delta_2d'):
    """Uses yt to trim and compare two files over a single field.
    Compared by subtracting the second file from the first and 
    dividing by the absolute maximum of the delta.
    Plots as a percentage from said peak.

    Args:
        filepair(str): filepaths (delta = first - second file).
        field(str): field to process.
        tags(str list): tag names for the file pair.
        resolution(float): pixel length of image (correlate with width).
        width(float): physical extent of image (x-axis, y-axis is 2x).
        trim(float): cutoff for yt data to reduce memory usage.
        vmin/max(floats): range for colorbar.
        cmap/filetag(str): colormap and filetag for output.
        dpi(int): resolution of plotted figure.
        batch(bool): print to file toggle (does not return fig).
        

    Returns:
        (mpl.figure or None)

    """
    ds = yt.load(filepair[0])
    cds = yt.load(filepair[1])
    if field in dir(ytf):
        meta = getattr(ytf, '_' + field)
        ds.add_field(("flash", field),
                     function=getattr(ytf, field), **meta)
        cds.add_field(("flash", field),
                      function=getattr(ytf, field), **meta)
    # trim domain
    reg = ds.r[0:trim, -trim:trim]
    creg = cds.r[0:trim, -trim:trim]
    # subset for plot
    cent = (width/2, 0.0)
    res = resolution
    extent = [0.0, width/1e5, -width/1e5, width/1e5]
    slc = ds.slice('z', 0.0, data_source=reg)
    frb = slc.to_frb(width, (res, 2*res),
                     height=2*width, center=cent)
    cslc = cds.slice('z', 0.0, data_source=creg)
    cfrb = cslc.to_frb(width, (res, 2*res),
                       height=2*width, center=cent)
    
    fig, ax  = plt.subplots(dpi=dpi)
    final = np.flip(frb[field].v - cfrb[field].v, axis=0)
    final = final/np.max(abs(final))*100
    msh = ax.imshow(final, vmin=vmin, vmax=vmax,
                    extent=extent, cmap=cmap)
    ax.set_xlabel('y (km)');ax.set_ylabel('x (km)')
    tag = r'{} (raw $\delta$)'.format(field)
    tag = r'{} (% from max)'.format(field.capitalize())
    cbar = fig.colorbar(msh, ax=[ax], label=tag)  # , extend='both', extendfrac=0.03)
    times = [float(ds.current_time), float(cds.current_time)]
    ttl = "{0}({2:.3f}s) - {1}({3:.3f}s)".format(*tags, *times)
    ax.set_title(ttl, pad=10)
#     for axes in fig.axes:
    fig.axes[0].yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    fig.axes[0].xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    fig.axes[1].yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    if batch:
        stem = filepair[0]
        writeFig(fig, stem, filetag)
    else:
        return fig
