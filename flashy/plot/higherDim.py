from .globals import (np, os, AxesGrid, plt, log, reformatTag, SymLogNorm,
                      ScalarFormatter, customFormatter, writeFig)
from flashy.datahaul.hdf5yt import getFields, yt
from flashy.IOutils import pairGen
from yt.funcs import mylog  # avoid yt warnings
import flashy.datahaul.ytfields as ytf
import flashy.datahaul.ytDatasetTools as ytt
import flashy.datahaul.tmat as tmat
mylog.setLevel(50)
log.name = __name__


def slice_cube(fname, grids=False, batch=False,
               frame=1e9, center=[0.0, 0.0, 0.0], sliceAxis='z',
               fields=['density', 'pressure', 'temperature'],
               linear=False, mins=[1.0, 1e+18, 1e7],
               maxs=[6e7, 3e+25, 8e9], mark=[], cmap='', filetag='ytprops'):
    """YT 2D plots of a specified list of fields through a defined slice.

    Args:
        fname(str): filename to plot.
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): total width of frame square.
        center(float list): center for slice.
        sliceAxis(str): normal axis for slice.
        fields(list of str): list of named fields to plot.
        linear(bool): set linear or log scale(false).
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        mark(float list): mark a (x,y) coordinate in the plot.
        cmap(str): matplotlib colormap for the plot.
        filetag(str): override default filetag (ytprops).

    Returns:
        (mpl.figure or None)

    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(7*size, 7))
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    p = yt.SlicePlot(ds, sliceAxis, list(fields),
                     origin='native', center=center)
    p.set_width((frame, frame))
    p.set_axes_unit('cm')
    dnum = {'x': 0, 'y': 1, 'z': 2}[sliceAxis]
    header = '{:3.5f} s | {}-slice @ {:.2e}'
    if mark:
        xm, ym = mark
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color': 'green', 's': 30,
                                     'facecolors': "None"})
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs)
    for i, (f, mi, mx) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(float(ds.current_time),
                                           sliceAxis.capitalize(),
                                           center[dnum]))
        if cmap:
            p.set_cmap(f, cmap)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if linear:
            p.set_log(f, False)
        plot = p.plots[f]
        # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    if not batch:
        return fig
    else:
        writeFig(fig, os.path.join(ds.fullpath, ds.basename), filetag)


def metaProps(fname, mhead=False, grids=False, batch=False, frame=1e9,
              center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
              linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
              maxs=[6e7, 3e+25, 8e9], mark=[], cmaps=['RdBu_r']*3,
              linthresh=1e15, dpi=90, fac=8, comm='', f4=False,
              metaprops={'shelld':1e6, 'cored':3e6, 'enucthresh':1e19}):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis. Batch build a meta plaintext file with
    the structure:
    
    time timedelta
    maxspeed xlocation ylocation
    # integrated enuc over a threshold
    (energy release in dt over threshold) threshold
    # calculate peak energies for shell and core, discerning by density
    shelldens coredens
    'shell' (summed energy) (peak release) X Y dens T velx vely
    'core' (summed energy) (peak release) X Y dens T velx vely
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
        mark(float list): mark a (x,y) coordinate in the plot.
        cmaps(str): matplotlib colormap for the plot.
        linthresh(float): symlog linea area around 0.
        dpi(float): dpi of figure returned/saved.
        fac(int): custom formatter exponent factor(cm to km = 5).
        comm(str): add custom mesage to plot.
        f4(bool): backwards compatibility with F4.
        metaprops(dict): brittle meta file properties.

    Returns:
        (mpl.figure or None)

    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt. force speed
    # because max speed calculation uses it.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)

    # in situ calculations
    ad = ds.all_data()
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
    energyRelease = ad['enuc'].value*cell_masses
    delT = ds.parameters['dt']
    # find the lump sum of energy over a threshold
    energymask = ad['enuc'].value > metaprops['enucthresh']
    lump = np.sum(energyRelease[energymask]*delT)
    mline = '{:.10E} {:.10E}'.format(lump, metaprops['enucthresh'])
    metatxt.append(mline)
    # calculate energy maxima for both core and shell 
    mline = '{:.10E} {:.10E}'.format(metaprops['shelld'], metaprops['cored'])
    metatxt.append(mline)
    mask = ad['dens'].value < metaprops['shelld']
    mepos = ad['enuc'][mask].argmax()
    medens, metemp = ad['dens'][mask].value[mepos], ad['temp'][mask].value[mepos]
    melocxS, melocyS = ad['x'][mask].value[mepos], ad['y'][mask].value[mepos]
    mevelx, mevely = ad['velx'][mask].value[mepos], ad['vely'][mask].value[mepos]
    emax = ad['enuc'][mask].value[mepos]*cell_masses[mask][mepos]
    toten = np.sum(energyRelease[mask&energymask]*delT)
    
    mline = 'shell {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E}'
    mline = mline.format(toten, emax, melocxS, melocyS, medens, metemp, mevelx, mevely)
    metatxt.append(mline)
    mask = ad['dens'].value > metaprops['cored']   
    try:
        mepos = ad['enuc'][mask].argmax()
        medens, metemp = ad['dens'][mask].value[mepos], ad['temp'][mask].value[mepos]
        melocxC, melocyC = ad['x'][mask].value[mepos], ad['y'][mask].value[mepos]
        mevelx, mevely = ad['velx'][mask].value[mepos], ad['vely'][mask].value[mepos]
        emax = ad['enuc'][mask].value[mepos]*cell_masses[mask][mepos]
        toten = np.sum(energyRelease[mask&energymask]*delT)
    except ValueError:
        # eventually all the domain is below threshold so set this equal to "shell"
        melocxC, melocyC = melocxS, melocyS
    mline = 'core {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E}'
    mline = mline.format(toten, emax, melocxC, melocyC, medens, metemp, mevelx, mevely)
    metatxt.append(mline)

    if f4:
        fields = [f.strip() for f in fields]
    for f in fields:
        fstr = '{} range: {:.10E} {:.10E}'
        exts = fstr.format(f, *ad.quantities.extrema(f).value)
        log.debug(exts)
        metatxt.append(exts)
    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family':'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color': 'black', 's': 30})
    if mark:
        xm, ym = mark
        log.warning('mark (green): {:E} {:E}'.format(xm, ym))
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color': 'green', 's': 250, 'linewidth':2,
                                     'facecolors': "None"})
    log.warning('max speed (white): {:E} {:E}'.format(mvlocx, mvlocy))
    p.annotate_marker((mvlocx, mvlocy), coord_system='plot', marker='o',
                      plot_args={'color': 'white', 's': 250, 'linewidth':2,
                                 'facecolors': "None"})
    log.warning('shell max energy output (red): {:E} {:E}'.format(melocxS, melocyS))
    p.annotate_marker((melocxS, melocyS), coord_system='plot', marker='o',
                      plot_args={'color': 'yellow', 's': 250, 'linewidth':2,
                                 'facecolors': "None"})
    log.warning('core max energy output (red): {:E} {:E}'.format(melocxC, melocyC))
    p.annotate_marker((melocxC, melocyC), coord_system='plot', marker='o',
                      plot_args={'color': 'red', 's': 250, 'linewidth':2,
                                 'facecolors': "None"})
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
                p.set_log(f, True , linthresh=linthresh)
            else:
                p.set_log(f, True)    
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()

    extras = {'maxspeed': maxsp}
    datb, datb2 = ytt.get_plotTitles(ds, comment=comm, extras=extras)
    fig.axes[0].set_title(datb)
    fig.axes[1].set_title(datb2)

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            ax.set_ylabel(u'y ($10^{{{}}}$ km)'.format(fac-5))
    for i, ax in enumerate(fig.axes):
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        filetag = 'ytmeta'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def speeds2D(bview, fname, subset=[0.0, 1e4, -1e4, 1e4], wedges=8,
             vmin=-1e5, vmax=1e5, batch=False, nticks=2):
    """plots all velocity variables from a 2d file.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file path.
        subset(float list): domain window for plotting in km.
        wedges(int): slices for data extraction.
        vmin(float): minimum velocity.
        vmax(float): maximum velocity.
        batch(bool): print to file toggle.
        nticks(int): ticks for xaxis.

    Returns:
        (mpl.figure or None)

    """
    fields = ['velx', 'vely', 'speed']
    dat = []
    for f in fields:
        t, dt, ext, dattmat = tmat.get2dPane(bview, fname, f, wedges=wedges)
        dat.append(dattmat)
    fig, axs = plt.subplots(nrows=1, ncols=3, dpi=150,
                            sharey=True, sharex=True, constrained_layout=True)
    fig.suptitle('Simtime: {:.3f}'.format(t))
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(dat[i]/1e5, cmap='RdBu_r',
                               extent=[x/1e5 for x in ext],
                               norm=SymLogNorm(linthresh=1e2,
                                               linscale=1.0,
                                               vmin=vmin, vmax=vmax))
        ax.set_title(fields[i])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_ticks(np.arange(0.0, np.max(subset),
                                     np.max(subset)/nticks))
        if not i:
            ax.set_ylabel('km')
        ax.set_xlabel('km')
    cbar = fig.colorbar(mshow, ax=[axs[:]],
                        location='bottom', extend='neither')
    cbar.set_label(u'$km/s$')
    if not batch:
        return fig
    else:
        filetag = 'speeds'
        writeFig(fig, fname, filetag)


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


def delt_runs_2D(fname, fields=['speed', 'velx', 'vely'],
                 compdir='', cmap='RdBu_r',
                 subset=[0.0, 1e9, -1e9, 1e9], linthresh=1e7,
                 vmin=-1e9, vmax=1e9, batch=False,
                 fac=8):
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
            if yval==tmat._nullval:
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
            for j, (r0, r1) in enumerate(pairGen(x[row]),1):
                mask = np.where((refx[i] >= r0) & (refx[i] <= r1))
                for k, d in enumerate(delts, 2):
                    d[i, mask] = dat[k][row, j] - shiftDat[k][row, j]
    log.debug('Deltas done. Plotting')
    cols = len(fields[2:])
    fig, axs = plt.subplots(nrows=1, ncols=cols, dpi=90,
                            sharey=True, sharex=True, constrained_layout=True)
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
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            ax.set_ylabel(u'y ($10^{{{}}}$ km)'.format(fac-5))
    
    if isinstance(axs, type(np.array)):
        cbar = fig.colorbar(mshow, ax=[axs[:]], location='right', extend='both')
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


def PARA_delt_runs_2D(bview, fname, fields=['speed', 'velx', 'vely'],
           compdir='', cmap='RdBu_r',
           subset=[0.0, 1e9, -1e9, 1e9], linthresh=1e7,
           vmin=-1e9, vmax=1e9, batch=False,
           nticks=4, wedges=8, fac=8):
    """plots deltas with contiguous checkpoint/plotfile
    for a list of fields.

    Args:
        bview(ipp.LoadBalancedView): ipp setup workhorse.
        fname(str): file path.
        fields(str list): fields to process.
        compdir(str): folderpath for comparison run.
        subset(float list): domain window for plotting.
        cmap(str): colormap name.
        linthresh(float): linear range around zero.
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
        t0, dt, ext, datmat = tmat.get2dPane(bview, fname, f, wedges=wedges)
        dat.append(datmat)
    shift = os.path.join(compdir, os.path.basename(fname))
    shiftDat = []
    for f in fields:
        t1, dt, ext, datmat = tmat.get2dPane(bview, shift, f)
        shiftDat.append(datmat)
    delts = [y-x for (x, y) in zip(dat, shiftDat)]
    fig, axs = plt.subplots(nrows=1, ncols=len(fields), dpi=90,
                            sharey=True, sharex=True, constrained_layout=True)
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(delts[i], extent=ext, cmap=cmap,
                               norm=SymLogNorm(linthresh=linthresh,
                                               linscale=0.8,
                                               vmin=vmin, vmax=vmax))
        ax.set_title(fields[i])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            ax.set_ylabel(u'y ($10^{{{}}}$ km)'.format(fac-5))
    cbar = fig.colorbar(mshow, ax=[axs[:]], location='bottom', extend='both')
    
    log.warning('assuming filename folder structure: {}'.format(fname))
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
