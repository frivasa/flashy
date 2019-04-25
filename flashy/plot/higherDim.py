# TDL separate datahaul elements from plotting routines
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
cmap = palettable.cmocean.sequential.Amp_20.mpl_colormap
# import cmocean fails for ipyparallel. 
# so hack it through palettable
_ytcmap = palettable.cmocean.sequential.Gray_20.mpl_colors
# _ytcmap = palettable.cmocean.diverging.Curl_19_r.mpl_colors
# _ytcmap = palettable.cmocean.diverging.Curl_19_r.mpl_colors
# _ytcmap = palettable.cmocean.diverging.Delta_20_r.mpl_colors
_ytcmap = palettable.cmocean.diverging.Balance_20.mpl_colors
# _ytcmap = palettable.cmocean.sequential.Amp_20.mpl_colors
cols = [tuple([list(x),1]) for x in _ytcmap]
cols.append(([0.0, 0.0, 0.0], 0))  # black initial color
# cols.append(([1.0, 1.0, 1.0], 0))  # white initial color
setcmap = yt.make_colormap(cols, name='custom')

# yt plotting functions
def planeSlice(fname, field, lims, zcut=0.0, linear=False, 
               batch=False, width=1.1e9, mark=[]):
    """Makes a slice at zcut
    
    Args:
        name(str): filename to extract data from.
        field(str): field to plot.
        lims(float list): field range.
        zcut(float): 3rd dimension cut depth.
        linear(bool): logscale toggle.
        batch(bool): toggle return figure or write to file.
        width(float): total side size of plot
        mark(float list): specify a position to mark in the plot.
    
    Returns:
        mpl.figure or None
    
    """
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
                      np.transpose(p.frb[field]),  # transpose aligns axes names to data.
                      lims)
    im.axes.tick_params(axis='x', top=True)
    im.axes.annotate("Time: {:.5f} s".format(float(ds.current_time)),
                     xy=(0.75, 0.05), xycoords='axes fraction', 
                     textcoords='axes fraction')
    im.axes.annotate("Z: {:.4e} cm".format(zcut),
                     xy=(0.10, 0.05), xycoords='axes fraction', 
                     textcoords='axes fraction')
    if mark:
        xm, ym = mark
        im.axes.annotate("o",
                         xy=mark[::-1], xycoords='data', 
                         textcoords='data')
    if batch:
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
        return None
    else:
        return fig
    

def bifold(fname, mhead=True, topfield='density', tlims=[1e-1, 4e7], width=1.1e9,
           botfield='temperature', blims=[3e7, 2e9], lintop=False,  batch=False,
           linbot=False, name=''):
    """Plots 2 properties for a hemisphere 2d cut as a joined sphere."""
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
    # fig.tight_layout()
    if batch:
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
    else:
        return fig


def plotFRB(gridAx, cbgAx, imArr, lims, top=True, linear=False):
    """draw frb to a given axesgrid.
    
    Args:
        gridAx(mpl.axes): axes to plot into.
        cbgAx(mpl.cbar): colorbar to set.
        imArr(np.array): data to plot.
        lims(float list): colorbar range.
        top(bool): tick removal for bifold plots.
        linear(bool): logscale toggle.
        
    Returns:
        mpl.axes, mpl.cbar
    
    """
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


def mainProps(fname, mhead=True, grids=False, batch=False, frame=1e9, center=(0.0, 0.0),
              fields=['density', 'pressure', 'temperature'], linear=False, 
              mins=[1.0, 1e+18, 1e7], maxs=[6e7, 3e+25, 8e9], mark=[], cmap=''):
    """Plots the list of fields specified in yt.
    
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
    grid = AxesGrid(fig, (0.075,0.075,0.85,0.85),
                    nrows_ncols = (1, size),
                    axes_pad = 1.2, label_mode = "L",
                    share_all = True, cbar_location="right",
                    cbar_mode="each", cbar_size="10%", cbar_pad="0%")
    p = yt.SlicePlot(ds, 'z', list(fields))
    if sum(center)==0.0:
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
                          plot_args={'color':'black', 's': 30})
    if mark:
        xm, ym = mark
        p.annotate_marker((xm, ym), coord_system='plot', marker='o',
                          plot_args={'color':'white', 's': 30, 'facecolors':"None" }) 
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
            p.set_log(field, False)
        plot = p.plots[f]
       # plot.figure = fig
        plot.axes = grid[i].axes
        plot.cax = grid.cbar_axes[i]
    p._setup_plots()
    if not batch:
        return fig
    else:
        filetag = 'ytprops'
        num = ds.parameter_filename[-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(os.path.dirname(ds.fullpath), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))


def flashSpeeds(fname, thresh=1e-6, filetag='speeds', batch=False, 
                 byM=True, direction=[], plotall=False):
    """Plot species distribution in 'speed' space for a set wedge in the domain.
    ATM the wedge is fixed at a polar angle of 0-45 degrees at 90 degrees azimuth. 
    TDL: variable wedge angle.
    
    Args:
        fname(str): path of file.
        thresh(float): threshold for species fraction.
        filetag(str): prefix for batch mode. 
        batch(bool): skips returning figure, saving it to a structured directory instead.
        byM(bool): plot by mass instead of radius.
    
    Returns:
        (mpl figure) or (None)
    
    """
    ds = yt.load(fname)
    box = ds.box(ds.arr([-3e9, 0.0, -1e8], "code_length"), ds.domain_right_edge)
    # left edge is extended in -z to include the central cell
    cylinder = ds.disk(center=ds.arr([0.0, -3.0e9, 3.0e9], "code_length"),
                       normal=ds.arr([0.0, -1.0/1.414213562, 1.0/1.414213562], "unitary"),
                       radius=(1e11, "code_length"),
    #                    height=(3e9*1.414213562, "code_length"))
                       height=(3e9*1.424213562, "code_length"))  # this is slightly larger (0.7%) to include edge cells
    wedge = ds.intersection([box, cylinder])
    
    _, species = getFields(ds.field_list)
    fields = ['velx', 'vely', 'velz', 'cell_mass'] + species
    
    # ask yt for data, as always this takes forever.
    rawd = []
    for f in fields:
        rawd.append(wedge[f].value)
    
    f, ax = plt.subplots()
    for i in range(len(species)):
        spmass = [x[1][i] for x in celltpls]
        joined = sorted(zip(speedrange, spmass))
        speeds, masses = zip(*joined)
        start = 0
        weights = []
        for c in counts:
            massf = sum(masses[start:start+c])
            weights.append(1.0/c)
            start+=c
        weights = np.array(weights)
        bins = bins[1:]
        mpln, mplbins, patches = ax.hist(bins, bins=len(bins), weights=weights, 
                                         histtype='step', log=True, label=species[i])
    # plot mass vs speed
    ax.set_ylim([thresh, 2])
    ax.yaxis.set_minor_formatter(StrMethodFormatter(''))
    ax.xaxis.set_major_formatter(customFormatter(10))
    ax.xaxis.set_minor_formatter(StrMethodFormatter(''))
    ax.set_xlabel('Speed ($10^{10}$ cm/s)')
    ax.legend()
    fig.subplots_adjust(hspace=0.4)
    
    if not batch:
        return fig
    else:
        num = ds.parameter_filename[-5:]  # checkpoint number 'flash_hdf5_chk_0001'
        dest = os.path.join(os.path.dirname(ds.fullpath), filetag)
        name = os.path.join(dest, '{}{}.png'.format(filetag, num))
        os.makedirs(dest, exist_ok=True)  # bless you, p3
        plt.savefig(name, format='png')
        plt.close(fig)
        print("Wrote: {}".format(name))

### Deprecated ###
def get_frb(fname, field, zcut=0.0, width=3e9, turn=False):
    """TODO: swtich to new method, this one is dead
    returns frb with required field."""
    ds = yt.load(fname)
#     p = yt.SlicePlot(ds, 'z', [field])
    p = yt.SlicePlot(ds, 'z', field, center=[0, 0, zcut])
    
    p.set_width((width, width))
    p.set_cmap(field, 'custom')
    p.set_axes_unit('cm')
    
    #p.set_width((width, 2*width))
    #p.set_center((width*0.5, 0.0))
    #p.set_origin(("center", "left", "domain"))
    p.set_axes_unit('cm')
    if turn:
        return np.transpose(p.frb[field])
    else:
        return p.frb[field]