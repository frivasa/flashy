from .globals import (np, os, AxesGrid, plt, log, reformatTag, SymLogNorm,
                      ScalarFormatter, StrMethodFormatter,
                      writeFig, resizeText)
from flashy.datahaul.hdf5yt import getFields, yt
from flashy.IOutils import pairGen, getFileList
from flashy.datahaul.plainText import dataMatrix
from yt.funcs import mylog  # avoid yt warnings
import flashy.datahaul.ytfields as ytf
import flashy.datahaul.ytDatasetTools as ytt
from flashy.plot.oneDim import plotSpecies
mylog.setLevel(50)
log.name = __name__


def radialTracing(fname, radius, diskradius,
                  unkpeak='enuc', invertpeak=False,
                  limspec='c12', limfrac=0.02, cylinR=1e9, otptag='radTrac',
                  grids=False, batch=False, frame=1e9, center=(0.0, 0.0),
                  fields=['density', 'pressure', 'temperature'],
                  linear=[False, False, False], mins=[1.0, 1e+18, 1e7],
                  maxs=[6e7, 3e+25, 8e9], cmaps=['RdBu_r']*3,
                  linthresh=1e15, dpi=90, comm=''):
    """plots 2d properties through YT, writing to a meta file the properties of
    all radial cells around the maximum unk variable picked (unkpeak) within a disk radius.
    
    if limspec is 'radius', limfrac becomes the outer radius of an annulus.
    cylinR has a delta at zero to avoid edge effects (deltaR)
    
    time timedelta radius diskradius (minimum c12)
    # (peak unk) (peak cell mass) X Y velx vely dens temp (species)
    note: values are unsorted, sort by y should yield an ordered 'lineout'

    Args:
        fname(str): filename to plot.
        radius(float): exclusion radius (pick points over this).
        diskradius(float): radius around max to print wake.
        unkpeak(str): unk variable to find peak for.
        invertpeak(bool): flip unpeak comparison.
        limspec(str): species to bound outer ring.
        limfrac(float): species fraction threshold (reject points under X_i).
        cylinR(float): limit annulus to an offset from the polar axis.
        otptag(str): filename tag for batch output.
        grids(bool): overplot the grid structure.
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
        comm(str): add custom mesage to plot.

    Returns:
        (mpl.figure or None)

    """
    deltaR = 100e5
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=dpi)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size), axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    metatxt = []
    # check for custom fields and add them to yt. force speed
    # because max speed calculation uses it.
    for f in fields + ['speed']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # start meta file
    delT = ds.parameters['dt']
    mline = '{:.10E} {:.10E} {:.10E} {:.10E} {:.10E}'
    mline = mline.format(float(ds.current_time), delT, radius, diskradius, limfrac)
    metatxt.append(mline)
    ad = ds.all_data()

    # find max enuc outside a radius and inside a minimum picked species fraction
    radii = np.sqrt(ad['x'].value**2 + ad['y'].value**2)
    over = np.where(radii > radius)[0]
    if limspec=='radius':
        rmask = (radii < limfrac) & (ad['x'].value < cylinR) & (ad['x'].value > deltaR)
        under = np.where(rmask)[0]
    else:
        under = np.where(ad[limspec.ljust(4)] > limfrac)
    donut = np.intersect1d(over, under)
    if invertpeak:
        mpos = ad[unkpeak][donut].argmin()
    else:
        mpos = ad[unkpeak][donut].argmax()
    mpos = np.arange(radii.shape[0])[donut][mpos]
    mlocx, mlocy = ad['x'].value[mpos], ad['y'].value[mpos]
    mvelx, mvely = ad['velx'].value[mpos], ad['vely'].value[mpos]

    # get cell masses
    log.warning('Assuming cylindrical slice')
    dx = ad['path_element_x'].value
    dy = ad['path_element_y'].value
    r = ad['x'].value
    cylvol = 2.0*np.pi*dy*dx*r
    cell_masses = cylvol*ad['density'].value
    _, sps = getFields(ds.field_list)
    ml = '# {} '.format(unkpeak) +\
         'cell_mass x y velx vely ' +\
         'dens temp pres gamc eint enuc ' +\
         ' '.join(sps)
    metatxt.append(ml)

    # build a disk around the maximum and
    # pick cells a dx/2 away from the max radius
    prad = np.sqrt(mlocx**2+mlocy**2)
    spot = ds.disk([mlocx, mlocy, 0.0], [0, 0, 1],
                   (diskradius, 'cm'), (1e5, 'cm'))
    spot_r = np.sqrt(spot['x'].v**2 + spot['y'].v**2)
    tol = np.nanmin(spot['dx'].v)*0.5
    dx = spot['path_element_x'].v
    dy = spot['path_element_y'].v
    r = spot['x'].v
    cylvol = 2.0*np.pi*dy*dx*r
    sub_cell_masses = cylvol*spot['density'].v
    ring = np.where((spot_r > prad-tol) & (spot_r < prad+tol))[0]
    for rp in ring:
        en, cm = spot[unkpeak].v[rp], sub_cell_masses[rp]
        x, y = spot['x'].v[rp], spot['y'].v[rp]
        velx, vely = spot['velx'].v[rp], spot['vely'].v[rp]
        dens, temp = spot['temperature'].v[rp], spot['density'].v[rp]
        pres, gamc = spot['pressure'].v[rp], spot['gamc'].v[rp]
        eint, enuc = spot['eint'].v[rp], spot['enuc'].v[rp]
        ml = '{:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} {:.10E} '
        ml += '{:.10E} {:.10E} {:.10E} {:.10E} '
        ml = ml.format(en, cm, x, y, velx, vely, dens, temp, pres,
                       gamc, eint, enuc)
        mls = ['{:.6E}'.format(spot[s.ljust(4)].value[rp]) for s in sps]
        ml += ' '.join(mls)
        metatxt.append(ml)

    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    p.set_axes_unit('km')
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    for rp in ring:
        ppx, ppy = spot['x'].value[rp]/1e5, spot['y'].value[rp]/1e5
        p.annotate_marker((ppx, ppy), coord_system='plot', marker='o',
                          plot_args={'color': 'white', 's': 200,
                                     'linewidth': 1, 'facecolors': "None"})
    p.annotate_marker((mlocx/1e5, mlocy/1e5), coord_system='plot', marker='o',
                      plot_args={'color': 'yellow', 's': 500, 'linewidth': 2,
                                 'facecolors': "None"})
    p.annotate_sphere([0.0, 0.0, 0.0], radius=(radius, 'cm'),
                      circle_args={'color': 'green'})
    p.annotate_sphere([0.0, 0.0, 0.0], radius=(prad, 'cm'),
                      circle_args={'color': 'yellow'})
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
                p.set_log(f, True, linthresh=linthresh)
            else:
                p.set_log(f, True)
        plot = p.plots[f]
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    extras = {'Speed': np.sqrt(mvelx**2 + mvely**2)}
    datb, datb2 = ytt.get_plotTitles(ds, comment=comm, extras=extras)
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
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 otptag, meta='\n'.join(metatxt))


def pointTracing(fname, field='density', low=1e5, high=5e8,
                 frame=400e5, batch=False, t0=1.2, t1=1.8, dpi=180):
    """Plots surrounding area of a point along with bulk properties
    and species.
    Looks for .meta files made by hdf5yt.getPointData for side plots.

    Args:
        fname(str): filename to plot.
        field(str): field to plot alongside composition.
        low(float): field minimum.
        high(float): field maximum.
        frame(float): field plot frame size
        batch(bool): write to file toggle.
        t0(float): initial time for composition plot.
        t1(float): final time for composition plot.
        dpi(int): dpi of figure.

    Returns:
        (mpl.figure or None)

    """
    sep = 0.15
    labelspace = 0.00
    filetag = 'pointTrack'
    metafold = os.path.join(os.path.dirname(os.path.dirname(fname)), filetag)
    legdict = {'ncol': 2, 'loc': 'upper left', 'columnspacing': 0.0,
               'labelspacing': 0.1, 'numpoints': 2, 'handletextpad': 0.2,
               'bbox_to_anchor': (1.05, 1.02), 'prop': {'size': 9}}
    flist = getFileList(metafold, glob='meta', fullpath=True)
    nums = np.array([int(f[-9:-5]) for f in flist])
    refn = int(fname[-4:])
    stop = np.where(refn >= nums)[0][-1]
    data = []
    for f in flist[:stop+1]:
        with open(f, 'r') as ff:
            header = ff.readline()
            dat = ff.readline()
            data.append([float(d) for d in dat.strip().split()])
    vkeys = header.strip('#\n').split()
    dm = dataMatrix([vkeys, np.array(data)])

    ds = yt.load(fname)
    fig = plt.figure(dpi=dpi, figsize=(10, 5))  # resolution and size
    # add axes by hand to maximize customizability and avoid reshaping bugs
    # by yt or subplots_adjust.
    ax1 = fig.add_axes([0.0, 0.0, 0.3, 1.0])
    ax2 = fig.add_axes([0.1, 0.0, 0.3, 1.0])
    ax3 = fig.add_axes([0.5, 0.0, 0.5, 0.5])
    ax4 = fig.add_axes([0.5, 0.5, 0.5, 0.5])

    # domain plot
    p = yt.SlicePlot(ds, 'z', [field])
    p.set_font({'family': 'monospace'})
    p.set_center((dm.posx[0], dm.posy[0]))
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    p.set_cmap(field, 'RdBu_r')  # fall back to RdBu
    p.set_zlim(field, low, high)
    adc = {'color': 'black', 's': 250, 'linewidth': 1, 'facecolors': "None"}
    p.annotate_marker((dm.posx[0]/1e5, dm.posy[0]/1e5),
                      coord_system='plot', marker='o',
                      plot_args=adc)
    plot = p.plots[field]
    plot.axes = ax1
    plot.cax = ax2
    p._setup_plots()

    # adjust cbar size and position
    mainp = ax1.get_position()
    cbp = ax2.get_position()
    cbp.x0 = mainp.x1
    cbp.x1 = mainp.x1 + 0.016
    cbp.y0 = mainp.y0
    cbp.y1 = mainp.y1
    ax2.set_position(cbp)
    p3 = ax3.get_position()
    p3.x0 = cbp.x0 + sep
    p3.x1 = 1.0-labelspace
    p3.y0 = 0.5
    p3.y1 = cbp.y1
    ax3.set_position(p3)
    p4 = ax4.get_position()
    p4.x0 = cbp.x0 + sep
    p4.x1 = 1.0-labelspace
    p4.y0 = cbp.y0
    p4.y1 = 0.5
    ax4.set_position(p4)

    # titles
    datb, datb2 = ytt.get_plotTitles(ds)
    ll = len(datb.split('\n')[-1])
    n, v, t = datb2.split()
    mid = str(ll - len(n) - 2)
    restr = '{}{:>' + mid + '} {}'
    restr = restr.format(n, v, t)
    ax1.set_title('\n'.join([datb, restr]))

    ax1.xaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))
    ax1.yaxis.set_major_formatter(StrMethodFormatter("{x:.0f}"))

    ax3.semilogy(dm.radius, dm.density*1e3, label=u'dens$\cdot 10^{3}$')
    ax3.semilogy(dm.radius, dm.pressure/1e15, label=u'pres$\cdot 10^{-15}$')
    ax3.semilogy(dm.radius, dm.enuc/1e9, label=u'enuc$\cdot 10^{-9}$')
    ax3.semilogy(dm.radius, dm.temperature, label='temp')

    ax3.set_xticklabels([])
    ax3.set_xlim([t0, t1])
    ax3.set_ylim([2e6, 1e10])
    ax3.legend(**legdict)

    # plotSpecies
    skip = plotSpecies(ax4, dm, byMass=False, thresh=1e-3,
                       marker=True, log=False)
    ax4.legend(**legdict)
    ax4.set_yscale('log')
    ax4.set_xlabel('time (s)')
    ax4.set_xlim([t0, t1])

    # adjust text
    for ax in [ax1, ax2, ax3, ax4]:
        resizeText(ax)
    if not batch:
        return fig
    else:
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag)


def debug_plot(fname, maxradius=4e8, grids=False, batch=False, frame=1e9,
               center=(0.0, 0.0), fields=['density', 'temperature'],
               linear=[False, False], mins=[1e3, 1e8], maxs=[6e6, 2e9],
               linthresh=1e10):
    """metaProps clone for debugging, focuses on max enuc, dens, temp
    minimum fields for image: 2
    probed fields = prbf below

    meta file structure:
    # probed fields:
    # probed fields
    # time timedelta
    time timedelta
    # value x y extra_named_fields
    probed fields yield 2 lines each: max and min

    Args:
        fname(str): filename to plot.
        maxradius(float): centered sphere radius to find maxima/minima
        grids(bool): overplot the grid structure.
        batch(bool): if true save figure to file instead of returning it.
        frame(float): physical extent of plot in x (twice for y) .
        center(float tuple): physical center of plot.
        fields(str list): list of named fields to plot.
        linear(bool list): set linear or log scale(false).
        mins(float list): minima of scale for each field.
        maxs(float list): maxima of scale for each field.
        linthresh(float): symlog linea area around 0.

    Returns:
        (mpl.figure or None)

    """
    cmaps = ['tab20c']*len(fields)
    if batch:
        print(fname, fields)
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 8), dpi=90)
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L",
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    # check for custom fields and add them to yt.
    for f in fields + ['speed', 'fermiDeg']:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)
    # build header and tabulate names
    # main edit chunk
    filetag = 'debug_plot'
    # fields to plot in the image and symbols
    prbf = ['density', 'temperature', 'enuc']
    signs = ['o', 'v', 's']
    # sampling values to write
    samp = ['density', 'temperature', 'enuc', 'eint',
            'c12 ', 'he4 ', 'fermiDeg']
    # auto mode
    metatxt = ['# probed fields:']
    mline = '# '+ ' '.join(prbf)
    metatxt.append(mline)
    mline = '# time delta'
    metatxt.append(mline)
    delT = ds.parameters['dt']
    mline = '{:.10E} {:.10E}'.format(float(ds.current_time), delT)
    metatxt.append(mline)
    mline = '# value x y '+ ' '.join(samp)
    metatxt.append(mline)
    ad = ds.sphere([0.0, 0.0, 0.0], maxradius)
    qu = ad.quantities

    mimark, mxmark = [], []
    for f in prbf:
        query = qu.max_location(f)
        mv, x, y, z = [v.value for v in query]
        txt = '{:.10E} {:.10E} {:.10E} '
        mline = txt.format(mv, x, y)
        query = qu.sample_at_max_field_values(f, sample_fields=samp)
        sampvals = [v.value for v in query]
        sampvals = sampvals[1:]  # min/max is repeated
        txt = '{:.10E} '*len(sampvals)
        mline += txt.format(*sampvals)
        metatxt.append(mline)
        mxmark.append((x, y))
        query = qu.min_location(f)
        mv, x, y, z = [v.value for v in query]
        txt = '{:.10E} {:.10E} {:.10E} '
        mline = txt.format(mv, x, y)
        query = qu.sample_at_min_field_values(f, sample_fields=samp)
        sampvals = [v.value for v in query]
        sampvals = sampvals[1:]  # min/max is repeated
        txt = '{:.10E} '*len(sampvals)
        mline += txt.format(*sampvals)
        metatxt.append(mline)
        mimark.append((x, y))

    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_font({'family': 'monospace'})
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('km')
    for (x, y), s, f in zip(mxmark, signs, prbf):
        p.annotate_marker((x/1e5, y/1e5), coord_system='plot', marker=s,
                          plot_args={'color': 'red', 's': 200, 'linewidth': 1,
                                     'facecolors': "None"})
    for (xmax, ymax), (xmin, ymin), s, f in zip(mxmark, mimark, signs, prbf):
        log.warning('{} is marker {}'.format(f, s))
        log.warning('max (red): {:E} {:E}'.format(xmax, ymax))
        log.warning('min (yellow): {:E} {:E}'.format(xmin, ymin))
#     p.annotate_contour('density', ncont=1, label=True, clim=(,))
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
    cmt = u"$\\rho=$O,$T=\\Delta$"
    datb, datb2 = ytt.get_plotTitles(ds, comment=cmt)
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
        writeFig(fig, os.path.join(ds.fullpath, ds.basename),
                 filetag, meta='\n'.join(metatxt))
