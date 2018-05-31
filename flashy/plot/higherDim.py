from .globals import *
import yt  # move this dependency to datahaul
# avoid yt warnings
from yt.funcs import mylog
mylog.setLevel(50)
# COLORMAPS
# https://jiffyclub.github.io/palettable/colorbrewer/
# colorscheme names: http://jiffyclub.github.io/palettable +/colorbrewer /matplotlib /cmocean
import palettable  # cmocean colors are included in palettable
# palettable schemes are not interpolated so make a new one for yt
# yt call: (name, type, num): num<=11 type={diverging, qualitative, sequential}
#_cScheme = ('GnBu', 'sequential', 9)
# import cmocean fails for ipyparallel. so hack it through palettable
cmap = palettable.cmocean.sequential.Tempo_20.mpl_colormap
# import cmocean fails for ipyparallel. 
# so hack it through palettable
#_ytcmap = palettable.cmocean.sequential.Ice_14_r.mpl_colors
#_ytcmap = palettable.cmocean.sequential.Gray_20.mpl_colors
_ytcmap = palettable.cmocean.sequential.Haline_10.mpl_colors
#_ytcmap = palettable.cmocean.diverging.Balance_20.mpl_colors
cols = [tuple([list(x),1]) for x in _ytcmap]
#cols.append(([0.0, 0.0, 0.0], 0))  # black initial color
cols.append(([1.0, 1.0, 1.0], 0))  # white initial color

setcmap = yt.make_colormap(cols, name='custom')

# yt plotting functions
def planeSlice(fname, field, lims, zcut=0.0, linear=False, show=False, width=1.1e9):
    """Makes a slice at zcut"""
    ds = yt.load(fname)
    p = yt.SlicePlot(ds, 'z', field, center=[0, 0, zcut])
    p.set_width((width, width))
    p.set_cmap(field, 'custom')
    p.set_axes_unit('cm')
    fig = plt.figure(figsize=(10,10))
    grid = AxesGrid(fig, 111, nrows_ncols=(1, 1),
                    axes_pad=0.0, label_mode="L",
                    share_all=False, cbar_location="right",
                    cbar_mode="each", cbar_size="5%", cbar_pad="0%")
    im, cax = plotFRB(grid[0], grid.cbar_axes[0], 
                      np.transpose(p.frb[field]), 
                      lims, top=False)
    im.axes.tick_params(axis='x', top=True)
    im.axes.annotate("Time: {:.5f} s".format(float(ds.current_time)),
                     xy=(0.75, 0.05), xycoords='axes fraction', 
                     textcoords='axes fraction')
    im.axes.annotate("Z: {:.4e} cm".format(zcut),
                     xy=(0.10, 0.05), xycoords='axes fraction', 
                     textcoords='axes fraction')
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'slices{}'.format(num)
        savn = 'slice_{:+13.0f}.png'.format(zcut)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print("Wrote: {}".format(savp))


def bifold(fname, mhead=True, topfield='density', tlims=[1e-1, 4e7], width=1.1e9,
           botfield='temperature', blims=[3e7, 2e9], lintop=False,  show=False,
           linbot=False, name=''):
    """Plots 2 properties as a joined sphere."""
    ds = yt.load(fname)
    p = yt.SlicePlot(ds, 'z', [topfield, botfield])
    p.set_width((width, 2*width))
    p.set_center((width*0.5, 0.0))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    fig = plt.figure(figsize=(12,9))
    #grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
    grid = AxesGrid(fig, 121,
                    nrows_ncols = (2, 1),
                    axes_pad = 0.0, label_mode = "L",
                    share_all = False, cbar_location="right",
                    cbar_mode="each", cbar_size="5%", cbar_pad="0%")
    plotFRB(grid[0], grid.cbar_axes[0], 
                  np.flipud(np.transpose(p.frb[topfield])), 
                  tlims, linear=lintop)
    plotFRB(grid[1], grid.cbar_axes[1], 
                  np.fliplr(np.transpose(p.frb[botfield])), 
                  blims, 
                  top=False, linear=linbot)
    if mhead:
        # we've transposed the matrix, so flip the _match positions
        x_match = ds.parameters['y_match']
        y_match = ds.parameters['x_match']
        grid[0].axes.plot(x_match, y_match, 'o', color='black')
        grid[1].axes.plot(x_match, y_match, 'o', color='black')
    grid[1].axes.annotate("Time: {:.5f} s".format(float(ds.current_time)),
                          xy=(0.75, -0.25), xycoords='axes fraction', 
                          textcoords='axes fraction')
#     fig.tight_layout()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        if name:
            tag = name
        else:
            tag = 'bifold'
        savn = '{}_{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp, bbox_inches='tight')
        plt.close(fig)
        print("Wrote: {}".format(savp))


def plotFRB(gridAx, cbgAx, imArr, lims, top=True, linear=False):
    unit = imArr.units
    field = imArr.info['field']
    yl, yr = imArr.info['ylim']
    xl, xr = imArr.info['xlim']
    vmin, vmax = lims
    if linear:
        norm = NoNorm(vmin=vmin, vmax=vmax)
    else:
        norm = LogNorm(vmin=vmin, vmax=vmax)
    if top:
        im = gridAx.imshow(imArr, 
                           cmap=cmap, aspect='auto', 
                           extent=[yl, yr, xl, xr], # horizontal mi/ma vertical mi/ma
                           norm=norm)
        yticks = im.axes.yaxis.get_major_ticks()
        yticks[0].label.set_visible(False)
        im.axes.tick_params(axis='x', bottom=False)
    else:
        im = gridAx.imshow(imArr, 
                           cmap=cmap, aspect='auto', 
                           extent=[yr, yl, xr, xl], norm=norm)
        im.axes.tick_params(axis='x', top=False)
    
    # get order of mag of highest y value and rescale the ticks
    #print im.get_extent()
    scaling = int(np.log10(abs(xr.value)))
    #scaling = int(np.log10(abs(im.get_extent()[3].value)))
    im.axes.yaxis.set_major_formatter(customFormatter(scaling))
    # yt flips axes, plus imarray is transposed, this must be an hdf5 saving issue
    im.axes.set_ylabel(r'x ($10^{{{}}}$ cm)'.format(scaling))
    
    # get order of mag of highest x value and rescale the ticks
    #scaling = int(np.log10(abs(im.get_extent()[1].value)))
    scaling = int(np.log10(abs(yr.value)))
    im.axes.xaxis.set_major_formatter(customFormatter(scaling))
    # yt flips axes, plus imarray is transposed, this must be an hdf5 saving issue
    im.axes.set_xlabel(r'y ($10^{{{}}}$ cm)'.format(scaling))
    
    scaling = int(np.log10(max(im.get_clim())))
    cax = cbgAx.colorbar(im, format=customFormatter(scaling, width=0, prec=2))
    cax.set_label_text('{} ($10^{{{}}} {}$)'.format(field.capitalize(),
                                                scaling,  
                                                str(unit).replace('**', '^')),
                        fontsize=11.0)
    yticks = cax.ax.get_yticklabels()
    yticks[-1].set_visible(False)
    yticks[0].set_visible(False)
    return im, cax


def mainProps(fname, mhead=True, grids=False, show=False,
              fields=['density', 'pressure', 'temperature'], linear=False, 
              mins=[0.1, 1e+16, 3e7], maxs=[4e7, 6e+24, 2e9]):
    """Plots the list of fields specified in yt.
    
    Args:
        fname (str): filename to plot.
        mhead (bool): mark the position of the matchhead.
        grids (bool): overplot the grid structure.
        show (bool): return figure (true) or save to file (false)
        fields (list of str): list of named fields to plot.
        linear (bool): set linear or log scale(false).
        mins (list of float): minima of scale for each field.
        maxs (list of float): maxima of scale for each field.
    
    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 10))
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (1, size),
                    axes_pad = 1.2, label_mode = "L",
                    share_all = True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    p = yt.SlicePlot(ds, 'z', list(fields))
    #plR = 7e8
    #p.set_width((plR, 2*plR))
    #p.set_center((plR*0.5, plR*0.7))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    header = '{:17.6f} s'
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color':'black', 's': 30})
    if grids:
        p.annotate_grids()
    pvars = zip(fields, mins, maxs)
    for i, (f, mi, mx) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(float(ds.current_time)))
        p.set_cmap(f, 'custom')
        p.set_zlim(f, mi, mx)
        if linear:
            p.set_log(field, False)
        plot = p.plots[f]
       # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'ytprops'
        savn = '{}{}.png'.format(tag, num)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print("Wrote: {}".format(savp))


def colortest(fname, name, cmap=palettable.cmocean.diverging.Balance_20.mpl_colors,
              mhead=True, grids=False, show=False,
              fields=['density'], linear=False, 
              mins=[0.1], maxs=[4e7]):
    """Plots the list of fields specified in yt.
    
    Args:
        fname (str): filename to plot.
        cmap (list of tuple): specify mpl colorlist(yt interpolates by default).
        mhead (bool): mark the position of the matchhead.
        grids (bool): overplot the grid structure.
        show (bool): return figure (true) or save to file (false)
        fields (list of str): list of named fields to plot.
        linear (bool): set linear or log scale(false).
        mins (list of float): minima of scale for each field.
        maxs (list of float): maxima of scale for each field.
    
    """
    ds = yt.load(fname)
    size = len(fields)
    fig = plt.figure(figsize=(5*size, 10))
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (1, size),
                    axes_pad = 1.2, label_mode = "L",
                    share_all = True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    p = yt.SlicePlot(ds, 'z', list(fields))
    p.set_width((1.2e9, 2*1.2e9))
    p.set_center((1.2e9*0.5, 0.0))
    #plr = 6e8
    #p.set_width((plr, plr))
    #p.set_center((plr*0.5, plr*4.1))
    p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    header = '{} {:10.3f} s'
    if mhead:
        x_match = ds.parameters['x_match']
        y_match = ds.parameters['y_match']
        p.annotate_marker((x_match, y_match), coord_system='plot',
                          plot_args={'color':'black', 's': 30})
    if grids:
        p.annotate_grids()

    _ytcmap = cmap
    cols = [tuple([list(x),1]) for x in _ytcmap]
    cols.append(([0.0, 0.0, 0.0], 0))
    setcmap = yt.make_colormap(cols, name='custom')
    pvars = zip(fields, mins, maxs)
    for i, (f, mi, mx) in enumerate(pvars):
        if not i:
            p.annotate_title(header.format(name, float(ds.current_time)))
        p.set_cmap(f, 'custom')
        p.set_zlim(f, mi, mx)
        if linear:
            p.set_log(field, False)
        plot = p.plots[f]
       # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    if show:
        return fig
    else:
        num = ds.parameter_filename[-5:]
        otpf, _ = os.path.split(ds.fullpath)
        tag = 'colors'
        savn = '{}{}.png'.format(tag, name)
        savf = os.path.join(otpf, "png")
        savp = os.path.join(otpf, "png", tag, savn)
        # build filetree and show or save the figure
        if not os.path.exists(savf):
            os.mkdir(savf)
            os.mkdir(os.path.join(savf, tag))
        elif not os.path.exists(os.path.join(savf, tag)):
            os.mkdir(os.path.join(savf, tag))
        plt.savefig(savp)
        plt.close(fig)
        print("Wrote: {}".format(savp))

