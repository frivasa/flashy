# TDL separate datahaul elements from plotting routines
from .globals import *
from flashy.datahaul.hdf5yt import getFields, yt
from yt.funcs import mylog  # avoid yt warnings
mylog.setLevel(50)
import flashy.datahaul.ytfields as ytf
import flashy.paraMan as pman
from flashy.datahaul.tmat import get2dPane

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
        # checkpoint number 'flash_hdf5_chk_0001'
        num = ds.parameter_filename[-5:]
        dest = os.path.join(os.path.dirname(ds.fullpath), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def mainProps(fname, mhead=True, grids=False, batch=False, frame=1e9,
              center=(0.0, 0.0), fields=['density', 'pressure', 'temperature'],
              linear=False, mins=[1.0, 1e+18, 1e7], maxs=[6e7, 3e+25, 8e9],
              mark=[], cmap=''):
    """YT 2D plots of a specified list of fields through a slice
    perpedicular to the z-axis.

    Args:
        fname (str): filename to plot.
        mhead (bool): mark the position of the matchhead.
        grids (bool): overplot the grid structure.
        batch (bool): if true save figure to file instead of returning it.
        fields (list of str): list of named fields to plot.
        linear (bool): set linear or log scale(false).
        mins (float list): minima of scale for each field.
        maxs (float list): maxima of scale for each field.
        mark (float list): mark a (x,y) coordinate in the plot.
        cmap (str): matplotlib colormap for the plot.

    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 10))
    grid = AxesGrid(fig, (0.075, 0.075, 0.85, 0.85),
                    nrows_ncols=(1, size),
                    axes_pad=1.2, label_mode="L", 
                    share_all=True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    # check for custom fields and add them to yt.
    for f in fields:
        if f in dir(ytf):
            meta = getattr(ytf, '_' + f)
            yt.add_field(("flash", f), function=getattr(ytf, f), **meta)

    p = yt.SlicePlot(ds, 'z', list(fields))
    if sum(center) == 0.0:
        p.set_center((frame*0.5, 0.0))
    else:
        p.set_center(center)
    p.set_width((frame, 2*frame))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    header = '{:17.6f} s'
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color': 'black', 's': 30})
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
            p.annotate_title(header.format(float(ds.current_time)))
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
        filetag = 'ytprops'
        # checkpoint number 'flash_hdf5_chk_0001'
        num = ds.parameter_filename[-5:]
        dest = os.path.join(os.path.dirname(ds.fullpath), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def speeds2D(bview, fname, subset=[0.0, 1e9, -1e9, 1e9],
             vmin=1e3, vmax=1e9, batch=False):
    fields = ['velx', 'vely', 'speed']
    nticks = 2
    dat = []
    for f in fields:
        t, dt, ext, tmat = get2dPane(bview, fname, f)
        dat.append(tmat)
    fig, axs = plt.subplots(nrows=1, ncols=3, dpi=140
                            , sharey=True, sharex=True, constrained_layout=True)
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(dat[i], extent=ext, cmap='RdBu_r'
                               , norm=LogNorm(vmin=vmin, vmax=vmax))
        ax.set_title(fields[i])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_ticks(np.arange(0.0, np.max(subset), np.max(subset)/nticks))
    fig.colorbar(mshow, ax=[axs[:]], location='bottom')

    if not batch:
        return fig
    else:
        filetag = 'speeds'
        # checkpoint number 'flash_hdf5_chk_0001'
        num = fname[-5:]
        basedest = os.path.dirname(os.path.dirname(fname))
        dest = os.path.join(basedest, filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def delt2D(bview, fname, fields=['speed', 'velx', 'vely'],
           subset=[0.0, 1e9, -1e9, 1e9],
           vmin=-1e5, vmax=1e5, batch=False):
    nticks = 4
    dat = []
    for f in fields:
        t0, dt, ext, tmat = get2dPane(bview, fname, f)
        dat.append(tmat)
    shift = '{}{:04d}'.format(fname[:-4], int(fname[-4:])+1)
    shiftDat = []
    for f in fields:
        t1, dt, ext, tmat = get2dPane(bview, shift, f)
        shiftDat.append(tmat)
    dt = t1 - t0
    delts = [y-x for (x, y) in zip(dat, shiftDat)]
    acc = [x/dt for x in delts]
    fig, axs = plt.subplots(nrows=1, ncols=len(fields), dpi=140
                            , sharey=True, sharex=True, constrained_layout=True)
    fig.suptitle('Acceleration ({:.3e})'.format(t1))
    for i, ax in enumerate(axs.flat):
        mshow = axs[i].matshow(acc[i], extent=ext, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.set_title(fields[i])
        ax.axis(subset)
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_major_formatter(ScalarFormatter())
        ax.xaxis.set_ticks(np.arange(0.0, np.max(subset), np.max(subset)/nticks))
    fig.colorbar(mshow, ax=[axs[:]], location='bottom', extend='both')

    if not batch:
        return fig
    else:
        filetag = 'accel'
        # checkpoint number 'flash_hdf5_chk_0001'
        num = fname[-5:]
        basedest = os.path.dirname(os.path.dirname(fname))
        dest = os.path.join(basedest, filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))
