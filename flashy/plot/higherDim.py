from .globals import (np, os, AxesGrid, plt, log, reformatTag, SymLogNorm,
                      ScalarFormatter, customFormatter, writeFig, resizeText)
from flashy.datahaul.hdf5yt import getFields, yt
from flashy.IOutils import pairGen, getFileList
from flashy.datahaul.plainText import dataMatrix
from yt.funcs import mylog  # avoid yt warnings
import flashy.datahaul.ytfields as ytf
import flashy.datahaul.ytDatasetTools as ytt
import flashy.datahaul.tmat as tmat
from flashy.plot.oneDim import plotSpecies
mylog.setLevel(50)
log.name = __name__
_alphas = ['he4 ', 'c12 ', 'o16 ', 'ne20',
           'mg24', 'si28', 's32 ', 'ar36',
           'ca40', 'ti44', 'cr48', 'fe52', 'ni56']


def slice_cube(fname, grids=False, batch=False,
               frame=1e9, center=[0.0, 0.0, 0.0], sliceAxis='z',
               fields=['density', 'pressure', 'temperature'],
               linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
               maxs=[6e7, 3e+25, 8e9], mark=[], fac=8, cmap='',
               particles=False, pwidth=1e6, filetag='ytprops', dpi=80):
    """YT 2D plots of a specified list of fields through a defined slice.

    Args:
        fname(str): filename to plot.
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): total width of frame square.
        center(float list): center for slice.
        sliceAxis(str): normal axis for slice.
        fields(list of str): list of named fields to plot.
        linear(bool list): set linear or log scale(false) for each plot.
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        mark(float list): mark a (x,y) coordinate in the plot.
        fac(int): custom formatter exponent factor(cm to km = 5). 
        cmap(str): matplotlib colormap for the plot.
        particles(bool): overplot particles on the leftmost plot.
        pwidth(float): particle sliver thickness.
        filetag(str): override default filetag (ytprops).
        dpi(int): force dpi of figure

    Returns:
        (mpl.figure or None)

    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(7*size, 7), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
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
    pvars = zip(fields, mins, maxs, linear)
    for i, (f, mi, mx, lin) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(float(ds.current_time),
                                           sliceAxis.capitalize(),
                                           center[dnum]))
        if cmap:
            p.set_cmap(f, cmap)
        else:
            p.set_cmap(f, 'RdBu_r')  # fall back to RdBu
        p.set_zlim(f, mi, mx)
        if lin:
            p.set_log(f, False)
        else:
            p.set_log(f, True)  # some fields start linear
        plot = p.plots[f]
        # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    axes = {
        'x': ('y', 'z'),
        'y': ('z', 'x'),
        'z': ('x', 'y')
    }.get(sliceAxis, ('x', 'y'))
    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'{} (km)'.format(axes[0]))
        else:
            lbl = u'{} ($10^{{{}}}$ km)'.format(axes[0], fac-5)
            ax.set_xlabel(lbl)
        if not i:
            if particles:
                if 'plt' in fname:
                    pre, name = os.path.split(fname)
                    nsp = name.split('_')
                    if 'plt' in nsp:
                        del(nsp[-2])
                        nsp[-2] = 'part'
                    newname = '_'.join(nsp)
                    pfname = os.path.join(pre, newname)
                    dsp = yt.load(pfname)
                    ad = dsp.all_data(fields=[('all', 'particle_position')])
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
                ax.scatter(prx, pry, marker='.', s=7, c='green')
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'{} (km)'.format(axes[1], fac-5))
            else:
                lbl = u'{} ($10^{{{}}}$ km)'.format(axes[1], fac-5)
                ax.set_ylabel(lbl)
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        writeFig(fig, os.path.join(ds.fullpath, ds.basename), filetag)


def metaProps(fname, mhead=False, grids=False, batch=False, frame=1e9,
              center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
              linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
              maxs=[6e7, 3e+25, 8e9], mark=[], cmaps=['RdBu_r']*3,
              linthresh=1e15, dpi=90, fac=8, comm='', f4=False, display=False,
              metaprops={'shelld': 1e6, 'cored': 3e6, 'energyThresh': 1e19}):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis. Batch build a meta plaintext file with
    the structure:

    time timedelta
    maxspeed xlocation ylocation
    # integrated enuc over a threshold
    (sum of positive energy release) threshold
    # calculate peak energies for shell and core, discerning by density
    shelldens coredens
    'shell' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely
    'core' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely
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
        display(bool): plot without meta markers, skipping calculation.
        metaprops(dict): brittle meta file properties.

    Returns:
        (mpl.figure or None)

    """
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
    if f4:
        fields = [f.strip() for f in fields]
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            ds.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # in situ calculations
    if not display:
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
        energyRelease = ad['enuc'].value*cell_masses  # erg/s
        # lump sum of positive energy
        energymask = energyRelease > 0.0
        lump = np.sum(energyRelease[energymask]*delT)
        mline = '{:.10E} {:.10E}'.format(lump, 0.0)
        metatxt.append(mline)
        # calculate energy maxima for both core and shell
        mline = '{:.10E} {:.10E}'
        mline = mline.format(metaprops['shelld'], metaprops['cored'])
        metatxt.append(mline)
        mask = ad['dens'].value <= metaprops['shelld']
        mepos = ad['enuc'][mask].argmax()
        medens, metemp = ad['dens'][mask].v[mepos], ad['temp'][mask].v[mepos]
        melocxS, melocyS = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
        mevelx, mevely = ad['velx'][mask].v[mepos], ad['vely'][mask].v[mepos]
        emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
        toten = np.sum(energyRelease[mask & energymask]*delT)
        # emax in erg/cc/s
        mline = 'shell ' + '{:.10E} '*8
        mline = mline.format(toten, emax, melocxS, melocyS,
                             medens, metemp, mevelx, mevely)
        metatxt.append(mline)
        mask = ad['dens'].value > metaprops['cored']
        try:
            mepos = ad['enuc'][mask].argmax()
            medens = ad['dens'][mask].v[mepos]
            metemp = ad['temp'][mask].v[mepos]
            melocxC, melocyC = ad['x'][mask].v[mepos], ad['y'][mask].v[mepos]
            mevelx = ad['velx'][mask].v[mepos]
            mevely = ad['vely'][mask].v[mepos]
            emax = ad['enuc'][mask].v[mepos]*ad['density'][mask].v[mepos]
            toten = np.sum(energyRelease[mask & energymask]*delT)
        except ValueError:
            # eventually all the domain is below thresh so set this to "shell"
            melocxC, melocyC = melocxS, melocyS
        mline = 'core ' + '{:.10E} '*8
        mline = mline.format(toten, emax, melocxC, melocyC,
                             medens, metemp, mevelx, mevely)
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
    p.set_axes_unit('cm')
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color': 'green', 's': 100})
    if mark:
        xm, ym = mark
        log.warning('mark (green): {:E} {:E}'.format(xm, ym))
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color': 'green', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
    if not display:
        log.warning('max speed (white): {:E} {:E}'.format(mvlocx, mvlocy))
        p.annotate_marker((mvlocx, mvlocy), coord_system='plot', marker='o',
                          plot_args={'color': 'white', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
        log.warning('shell max otp (red): {:E} {:E}'.format(melocxS, melocyS))
        p.annotate_marker((melocxS, melocyS), coord_system='plot', marker='o',
                          plot_args={'color': 'yellow', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
        log.warning('core max otp (red): {:E} {:E}'.format(melocxC, melocyC))
        p.annotate_marker((melocxC, melocyC), coord_system='plot', marker='o',
                          plot_args={'color': 'red', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
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
                if any(le < ds.domain_left_edge.v):
                    le = ds.domain_left_edge.v
                if any(re < ds.domain_right_edge.v):
                    re = ds.domain_right_edge.v
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
        datb, datb2 = ytt.get_plotTitles(ds, comment=comm)
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

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
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


def mprops_split(fname, mhead=False, grids=False, batch=False, frame=1e9,
                 center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
                 linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                 maxs=[6e7, 3e+25, 8e9], mark=[], lineouts=[],
                 cmaps=['RdBu_r']*3,
                 linthresh=1e15, dpi=90, fac=8, comm='',
                 display=False, particles=False,
                 metaprops={'radius': 4e8}):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis. Batch build a meta plaintext file with
    the structure:

    time timedelta
    maxspeed xlocation ylocation
    # integrated enuc over a threshold
    (sum of positive energy release) threshold
    # peak energies for shell and core, discerning by limit radius
    radius radius
    'shell' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely cell_vol
    'core' (sum of >0 in zone) (peak erg/cc) X Y dens T velx vely cell_vol
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
        mark(float list): mark a (x,y) coordinate.
        lineouts(float list): mark lineouts from origin (list of angles).
        cmaps(str): matplotlib colormap for the plot.
        linthresh(float): symlog linea area around 0.
        dpi(float): dpi of figure returned/saved.
        fac(int): custom formatter exponent factor(cm to km = 5).
        comm(str): add custom mesage to plot.
        display(bool): plot without meta markers, skipping calculation.
        metaprops(dict): brittle meta file properties.

    Returns:
        (mpl.figure or None)

    """
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
        # emax in erg/cc/s
        mline = 'shell ' + '{:.10E} '*9
        mline = mline.format(toten, emax, melocxS, melocyS,
                             medens, metemp, mevelx, mevely, cvol)
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
        except ValueError:
            # eventually all the domain is below thresh so set this to "shell"
            melocxC, melocyC = melocxS, melocyS
        mline = 'core ' + '{:.10E} '*9
        mline = mline.format(toten, emax, melocxC, melocyC,
                             medens, metemp, mevelx, mevely, cvol)
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
    p.set_axes_unit('cm')
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color': 'green', 's': 100})
    if mark:
        xm, ym = mark
        log.warning('mark (green): {:E} {:E}'.format(xm, ym))
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color': 'green', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
    if not display:
        log.warning('max speed (white): {:E} {:E}'.format(mvlocx, mvlocy))
        p.annotate_marker((mvlocx, mvlocy), coord_system='plot', marker='o',
                          plot_args={'color': 'white', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
        log.warning('shell max otp (red): {:E} {:E}'.format(melocxS, melocyS))
        p.annotate_marker((melocxS, melocyS), coord_system='plot', marker='o',
                          plot_args={'color': 'yellow', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
        log.warning('core max otp (red): {:E} {:E}'.format(melocxC, melocyC))
        p.annotate_marker((melocxC, melocyC), coord_system='plot', marker='o',
                          plot_args={'color': 'red', 's': 250,
                                     'linewidth': 2, 'facecolors': "None"})
    p.annotate_sphere([0.0, 0.0, 0.0], radius=(metaprops['radius'], 'cm'),
                      circle_args={'color': 'white', 'ls': '--', 'linewidth': 1})
    if lineouts:
        for angle in lineouts:
            height = np.sin(np.radians(angle))*2.0*frame*np.sqrt(2)
            length = np.cos(np.radians(angle))*2.0*frame*np.sqrt(2)
            p.annotate_line((0.0, 0.0), (length, height), coord_system='plot',
                            plot_args={'color': 'black', 'ls': '--', 'linewidth': 1})
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
    p._setup_plots()
    
    if display:
        datb, datb2 = ytt.get_plotTitles(ds, comment=comm)
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

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
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
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
                ax.set_ylabel(u'y ($10^{{{}}}$ km)'.format(fac-5))
        # skip the leftmost ylabel
        if i:
            tag, changed = reformatTag(ax.yaxis.get_label_text())
            if changed:
                ax.set_ylabel(tag, {'rotation': 0})
    if not batch:
        return fig
    else:
        filetag = 'ytm_radius'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def contourProps(fname, grids=False, batch=False, frame=1e9,
                 center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
                 linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                 maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
                 linthresh=1e15, dpi=90, fac=8,
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
        fac(int): custom formatter exponent factor(cm to km = 5).
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
    p.set_axes_unit('cm')
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
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
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
        filetag = 'contour'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def sphereProps(fname, sphR=1e5, sphC=(0.0, 0.0, 0.0),
                grids=False, batch=False, frame=1e9,
                center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
                linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
                linthresh=1e15, dpi=90, fac=8):
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
    for el in _alphas:
        els.append(ad[el][msk].value)
    chungus = zip(ens, xs, ys, rds, dns, ts, vls, prs, 
                  mss, vxs, vys, sps, *els)
    mline = '# enuc x y r dens temp vol pres mass velx vely speed '
    mline += ' '.join(_alphas) 
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
    p.set_axes_unit('cm')
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
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
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
        filetag = 'sphereContour'
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))


def zoneProps(fname, ledge=(0.0, 0.0, 0.0), redge=(40e5, 50e5, 0.0175),
              grids=False, batch=False, frame=1e9,
              center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
              linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
              maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
              linthresh=1e15, dpi=90, fac=8):
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
    for el in _alphas:
        els.append(ad[el][msk].value)
    chungus = zip(ens, xs, ys, rds, dns, ts, vls, prs, 
                  mss, vxs, vys, sps, *els)
    mline = '# enuc x y r dens temp vol pres mass velx vely speed '
    mline += ' '.join(_alphas) 
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
    p.set_axes_unit('cm')
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
    # ledge ---- 
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
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
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
        filetag = 'squareContour'
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

        
def adhoc_plot(fname, batch=False, frame=1e9,
               center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
               linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
               maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
               linthresh=1e15, dpi=90, fac=8, comm=''):
    """same as meta props without meta file.

    Args:
        fname(str): filename to plot.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): physical extent of plot in x (twice for y) .
        center(float tuple): physical center of plot.
        fields(str list): list of named fields to plot.
        linear(bool list): set linear or log scale(false).
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        cmaps(str): matplotlib colormap for the plot.
        linthresh(float): symlog linea area around 0.
        dpi(float): dpi of figure returned/saved.
        fac(int): custom formatter exponent factor(cm to km = 5).
        comm(str): add custom mesage to plot.

    Returns:
        (mpl.figure or None)

    """
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

    ledge = (0.0, -frame, 0.0)
    redge = (frame, frame, 0.0175)
    box = ds.box(ledge, redge)
    p = yt.SlicePlot(ds, 'z', list(fields), data_source=box)
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
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
                if any(le < ds.domain_left_edge.v):
                    le = ds.domain_left_edge.v
                if any(re < ds.domain_right_edge.v):
                    re = ds.domain_right_edge.v
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
#     pargs = {'color': 'green', 's': 100}
#     p.annotate_marker((x_match, y_match), coord_system='plot', plot_args=pargs)
#     pargs = {'color': 'green', 's': 250, 'linewidth': 2, 'facecolors': "None"}
#     p.annotate_marker((xm, ym), coord_system='plot', marker='o', plot_args=pargs)
#     p.annotate_grids()
    cargs = {'color': 'white', 'ls': '--', 'linewidth': 1}
    p.annotate_sphere([0.0, 0.0, 0.0], radius=(4120e5, 'cm'), circle_args=cargs)
    pargs = {'colors': 'green', 'linewidths': 1, 'linestyles':'--'}
    p.annotate_contour("temperature", ncont=1, clim=(1e9,1e9), label=False,
                       text_args={'colors': 'k'},
                       plot_args=pargs)
    p._setup_plots()
#     datb, datb2 = ytt.get_plotTitles(ds, comment=comm)
#     fig.axes[0].set_title(datb2)
#     fig.axes[1].set_title(datb2)

    for i, ax in enumerate(fig.axes):
        ax.xaxis.set_major_formatter(customFormatter(fac, prec=0))
        if fac == 5:
            ax.set_xlabel(u'x (km)')
        else:
            ax.set_xlabel(u'x ($10^{{{}}}$ km)'.format(fac-5))
        if not i:
            timetag = '{:.3f}s'.format(float(ds.current_time))
            ax.set_title(timetag, loc='right')
            ax.yaxis.set_major_formatter(customFormatter(fac, prec=0))
            if fac == 5:
                ax.set_ylabel(u'y (km)'.format(fac-5))
            else:
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
        filetag = 'adhoc'
        filepath = os.path.join(ds.fullpath, ds.basename)
        writeFig(fig, filepath, filetag)